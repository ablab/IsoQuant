############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import namedtuple

from src.isoform_assignment import *
from src.gene_info import *
from src.long_read_profiles import *

logger = logging.getLogger('IsoQuant')


class PolyAVerifier:
    def __init__(self, gene_info, params):
        self.params = params
        self.gene_info = gene_info

    # check consistency with polyA
    def verify_read_ends(self, combined_read_profile, isoform_id, matching_events):
        logger.debug("+ Validating polyA/T sites")
        if isoform_id is None:
            logger.debug("+ No sites found, ciao")
            return matching_events

        logger.debug("+ Checking isoform %s" % isoform_id)
        read_exons = combined_read_profile.read_exon_profile.read_features
        polya_info = combined_read_profile.polya_info
        isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]

        if self.gene_info.isoform_strands[isoform_id] == '+':
            if polya_info.internal_polya_pos != -1:
                # check for internal intron polyA priming
                matching_events, is_internal = self.check_internal_polya(polya_info.internal_polya_pos, matching_events)
                if is_internal:
                    return matching_events

            if polya_info.external_polya_pos != -1:
                matching_events, corrected_read_end = self.verify_polya(isoform_exons, read_exons,
                                                                        polya_info, matching_events)
                combined_read_profile.corrected_read_end = corrected_read_end

        elif self.gene_info.isoform_strands[isoform_id] == '-':
            if polya_info.internal_polyt_pos != -1:
                # check for internal intron polyA priming
                matching_events, is_internal = self.check_internal_polyt(polya_info.internal_polyt_pos, matching_events)
                if is_internal:
                    return matching_events

            if polya_info.external_polyt_pos != -1:
                matching_events, corrected_read_start = self.verify_polyt(isoform_exons, read_exons,
                                                                          polya_info, matching_events)
                combined_read_profile.corrected_read_start = corrected_read_start

        if not matching_events:
            logger.warning("Empty event list after polyA verification")
            matching_events = [make_event(MatchEventSubtype.none)]
        return matching_events

    def verify_polya(self, isoform_exons, read_exons, polya_info, matching_events):
        isoform_end = isoform_exons[-1][1]
        event_to_remove = -1
        fake_terminal_exon_count = 0
        # checking events
        for i, event in enumerate(matching_events):
            if event.event_type in [MatchEventSubtype.major_exon_elongation_right,
                                    MatchEventSubtype.exon_elongation_right]:
                # substitute major elongation with APA site
                event_to_remove = i
            elif event.event_type == MatchEventSubtype.fake_terminal_exon_right:
                fake_terminal_exon_count += 1

        #assert fake_terminal_exon_count == 0 or event_to_remove == -1
        if event_to_remove != -1:
            del matching_events[event_to_remove]

        corrected_read_end = self.truncate_fake_terminal_exons(read_exons, fake_terminal_exon_count,
                                                               polya_info.external_polya_pos)
        matching_events, corrected_read_end = \
            self.check_reference_terminal_exons(isoform_exons, corrected_read_end, matching_events)

        # accidentally aligned polyA tails must not create APA site
        dist_external_to_polya = abs(isoform_end - corrected_read_end)
        dist_to_internal_polya = abs(isoform_end - polya_info.internal_polya_pos)
        dist_to_polya = min(dist_external_to_polya, dist_to_internal_polya)
        logger.debug("+ Distance to polyA is %d" % dist_to_polya)
        if dist_to_polya > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(make_event(MatchEventSubtype.alternative_polya_site, event_length=dist_to_polya))
        elif dist_to_polya > self.params.delta:
            logger.debug("+ Seems like minor exon elongation site")
            matching_events.append(make_event(MatchEventSubtype.exon_elongation_right, event_length=dist_to_polya))
        return matching_events, corrected_read_end

    def verify_polyt(self, isoform_exons, read_exons, polya_info, matching_events):
        isoform_start = isoform_exons[0][0]
        event_to_remove = -1
        fake_terminal_exon_count = 0
        for i, event in enumerate(matching_events):
            if event.event_type in [MatchEventSubtype.major_exon_elongation_left,
                                    MatchEventSubtype.exon_elongation_left]:
                # substitute elongation with APA site if major
                event_to_remove = i
            elif event.event_type == MatchEventSubtype.fake_terminal_exon_left:
                fake_terminal_exon_count += 1

        #assert fake_terminal_exon_count == 0 or event_to_remove == -1
        if event_to_remove != -1:
            del matching_events[event_to_remove]

        corrected_read_start = self.truncate_fake_starting_exons(read_exons, fake_terminal_exon_count,
                                                                 polya_info.external_polyt_pos)
        matching_events, corrected_read_start = \
            self.check_reference_starting_exons(isoform_exons, corrected_read_start, matching_events)

        # accidentally aligned polyT head must not create APA site
        dist_external_to_polya = abs(isoform_start - corrected_read_start)
        dist_to_internal_polya = abs(isoform_start - polya_info.internal_polyt_pos)
        dist_to_polya = min(dist_external_to_polya, dist_to_internal_polya)
        logger.debug("+ Distance to polyA is %d" % dist_to_polya)
        if dist_to_polya > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(make_event(MatchEventSubtype.alternative_polya_site, event_length=dist_to_polya))
        elif dist_to_polya > self.params.delta:
            logger.debug("+ Seems like minor exon elongation site")
            matching_events.append(make_event(MatchEventSubtype.exon_elongation_left, event_length=dist_to_polya))
        return matching_events, corrected_read_start

    # check if polyA found within intron and inside mapped part of the read
    def check_internal_polya(self, internal_polya_pos, matching_events):
        if internal_polya_pos == -1:
            return matching_events, False

        for i, event in enumerate(matching_events):
            if event.event_type == MatchEventSubtype.incomplete_intron_retention_right:
                matching_events.append(make_event(MatchEventSubtype.fake_polya_site,
                                                  isoform_position=event.isoform_position,
                                                  event_length=internal_polya_pos))
                return matching_events, True

        return matching_events, False

    # check if polyT found within intron and inside mapped part of the read
    def check_internal_polyt(self, internal_polyt_pos, matching_events):
        if internal_polyt_pos == -1:
            return matching_events, False

        fake_terminal_exon_count = 0
        for i, event in enumerate(matching_events):
            if event.event_type == MatchEventSubtype.incomplete_intron_retention_left:
                # polyT found within intron and inside mapped part of the read
                matching_events.append(make_event(MatchEventSubtype.fake_polya_site,
                                                  isoform_position=event.isoform_position,
                                                  event_length=internal_polyt_pos))
                return matching_events, True

        return matching_events, False

    # correct polyA position when fake terminal exons are present or isoform has short terminal exons
    def truncate_fake_terminal_exons(self, read_exons, fake_terminal_exon_count, polya_pos):
        if fake_terminal_exon_count == 0:
            return polya_pos

        assert fake_terminal_exon_count < len(read_exons)
        # correcting fake terminal exons
        dist_to_polya = 0
        for i in range(fake_terminal_exon_count):
            exon = read_exons[-i-1]
            if exon[0] > polya_pos:
                continue
            elif dist_to_polya == 0:
                # no exons counted yet
                dist_to_polya += polya_pos-exon[0]
            else:
                dist_to_polya += interval_len(exon)
        last_good_exon = read_exons[-fake_terminal_exon_count-1]
        return last_good_exon[1] + dist_to_polya

    # correct polyT position when fake terminal exons are present or isoform has short terminal exons
    def truncate_fake_starting_exons(self, read_exons, fake_terminal_exon_count, polyt_pos):
        if fake_terminal_exon_count == 0:
            return polyt_pos

        assert fake_terminal_exon_count < len(read_exons)
        # correcting fake terminal exons
        dist_to_polya = 0
        for i in range(fake_terminal_exon_count):
            exon = read_exons[i]
            if exon[1] < polyt_pos:
                continue
            elif dist_to_polya == 0:
                # no exons counted yet
                dist_to_polya += exon[1] - polyt_pos
            else:
                dist_to_polya += interval_len(exon)
        first_good_exon = read_exons[fake_terminal_exon_count]
        return first_good_exon[0] - dist_to_polya

    # check isoform exons beyond polyA
    def check_reference_terminal_exons(self, isoform_exons, polya_pos, matching_events):
        isoform_terminal_exon_count = 0
        for i in range(len(isoform_exons)):
            if isoform_exons[-i-1][0] < polya_pos:
                break
            isoform_terminal_exon_count += 1

        if isoform_terminal_exon_count == 0:
            return matching_events, polya_pos

        corrected_read_end = polya_pos
        isoform_terminal_exon_length = intervals_total_length(isoform_exons[-isoform_terminal_exon_count:])
        dist_to_polya = abs(isoform_exons[-isoform_terminal_exon_count-1][1] - polya_pos)
        if isoform_terminal_exon_length <= self.params.max_missed_exon_len and \
                dist_to_polya <= self.params.max_missed_exon_len:
            logger.debug("+ Looks like missed terminal exons, shifting polyA")
            corrected_read_end = isoform_exons[-1][1]
            for i in range(isoform_terminal_exon_count):
                matching_events.append(make_event(MatchEventSubtype.exon_misallignment, len(isoform_exons) - 2 - i))

        return matching_events, corrected_read_end

    # check isoform exons prior polyT
    def check_reference_starting_exons(self, isoform_exons, polyt_pos, matching_events):
        # check isoform exons beyond polyA
        isoform_terminal_exon_count = 0
        for i in range(len(isoform_exons)):
            if isoform_exons[i][1] > polyt_pos:
                break
            isoform_terminal_exon_count += 1

        if isoform_terminal_exon_count == 0:
            return matching_events, polyt_pos

        corrected_polyt = polyt_pos
        isoform_terminal_exon_length = intervals_total_length(isoform_exons[:isoform_terminal_exon_count])
        dist_to_polya = abs(isoform_exons[isoform_terminal_exon_count][0] - corrected_polyt)
        if isoform_terminal_exon_length <= self.params.max_missed_exon_len and \
                dist_to_polya <= self.params.max_missed_exon_len:
            logger.debug("+ Looks like missed terminal exons, shifting polyT")
            corrected_polyt = isoform_exons[0][0]
            for i in range(isoform_terminal_exon_count):
                matching_events.append(make_event(MatchEventSubtype.exon_misallignment, i))

        return matching_events, corrected_polyt

