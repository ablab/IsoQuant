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
        if (combined_read_profile.polya_pos == -1 and combined_read_profile.polyt_pos == -1) or isoform_id is None:
            logger.debug("+ No sites found, ciao")
            return matching_events

        logger.debug("+ Checking isoform %s" % isoform_id)
        read_exons = combined_read_profile.read_exon_profile.read_features
        isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]

        if self.gene_info.isoform_strands[isoform_id] == '+' and combined_read_profile.polya_pos != -1:
            matching_events, corrected_polya_pos = self.verify_polya(isoform_exons, read_exons,
                                                                     combined_read_profile.polya_pos, matching_events)
            combined_read_profile.polya_pos = corrected_polya_pos

        elif self.gene_info.isoform_strands[isoform_id] == '-' and combined_read_profile.polyt_pos != -1:
            matching_events, corrected_polyt_pos = self.verify_polyt(isoform_exons, read_exons,
                                                                     combined_read_profile.polyt_pos, matching_events)
            combined_read_profile.polyt_pos = corrected_polyt_pos

        if not matching_events:
            logger.warning("Empty event list after polyA verification")
            matching_events = [make_event(MatchEventSubtype.none)]
        return matching_events

    def verify_polya(self, isoform_exons, read_exons, polya_pos, matching_events):
        isoform_end = isoform_exons[-1][1]
        read_end = read_exons[-1][1]

        # FIXME: remove this set once we check these events never appear together
        events_to_remove = set()
        fake_terminal_exon_count = 0
        # checking events
        for i, event in enumerate(matching_events):
            if event.event_type == MatchEventSubtype.incomplete_intron_retention_right and \
                    read_end > polya_pos + self.params.polya_window:
                # TODO: check for polyA once again
                # polyA found within intron and inside mapped part of the read
                matching_events.append(make_event(MatchEventSubtype.fake_polya_site,
                                                  isoform_position=event.isoform_position,
                                                  event_length=read_end - polya_pos))
                return matching_events, polya_pos

            if event.event_type == MatchEventSubtype.major_exon_elongation_right:
                # substitute major elongation with APA site
                events_to_remove.add(i)
            elif event.event_type == MatchEventSubtype.fake_terminal_exon_right:
                fake_terminal_exon_count += 1

        assert fake_terminal_exon_count == 0 or len(events_to_remove) == 0
        if len(events_to_remove) > 1:
            logger.warning("Too many extension events on the right side: " + str(events_to_remove))
        if events_to_remove:
            del matching_events[events_to_remove.pop()]

        polya_pos = self.correct_polya_coord(read_exons, fake_terminal_exon_count, polya_pos)
        matching_events, polya_pos = self.check_reference_terminal_exons(isoform_exons, polya_pos, matching_events)

        dist_to_polya = abs(isoform_end - polya_pos)
        logger.debug("+ Distance to polyA is %d" % dist_to_polya)
        if dist_to_polya > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(make_event(MatchEventSubtype.alternative_polya_site, event_length=dist_to_polya))
        return matching_events, polya_pos

    def verify_polyt(self, isoform_exons, read_exons, polyt_pos, matching_events):
        isoform_start = isoform_exons[0][0]
        read_start = read_exons[0][0]

        events_to_remove = set()
        fake_terminal_exon_count = 0
        for i, event in enumerate(matching_events):
            if event.event_type == MatchEventSubtype.incomplete_intron_retention_left and \
                    read_start < polyt_pos - self.params.polya_window:
                # polyT found within intron and inside mapped part of the read
                matching_events.append(make_event(MatchEventSubtype.fake_polya_site,
                                                  isoform_position=event.isoform_position,
                                                  event_length=polyt_pos - read_start))
                return matching_events, polyt_pos

            if event.event_type == MatchEventSubtype.major_exon_elongation_left:
                # substitute major elongation with APA site
                events_to_remove.add(i)
            elif event.event_type == MatchEventSubtype.fake_terminal_exon_left:
                fake_terminal_exon_count += 1

        assert fake_terminal_exon_count == 0 or len(events_to_remove) == 0
        if len(events_to_remove) > 1:
            logger.warning("Too many extension events on the left side: " + str(events_to_remove))
        if events_to_remove:
            del matching_events[events_to_remove.pop()]

        polyt_pos = self.correct_polyt_coord(read_exons, fake_terminal_exon_count, polyt_pos)
        matching_events, polyt_pos = self.check_reference_starting_exons(isoform_exons, polyt_pos, matching_events)

        dist_to_polya = abs(isoform_start - polyt_pos)
        logger.debug("+ Distance to polyA is %d" % dist_to_polya)
        if dist_to_polya > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(make_event(MatchEventSubtype.alternative_polya_site, event_length=dist_to_polya))
        return matching_events, polyt_pos

    # correct polyA position when fake terminal exons are present or isoform has short terminal exons
    def correct_polya_coord(self, read_exons, fake_terminal_exon_count, polya_pos):
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

    # check isoform exons beyond polyA
    def check_reference_terminal_exons(self, isoform_exons, polya_pos, matching_events):
        isoform_terminal_exon_count = 0
        for i in range(len(isoform_exons)):
            if isoform_exons[-i-1][0] < polya_pos:
                break
            isoform_terminal_exon_count += 1

        corrected_polya = polya_pos
        if isoform_terminal_exon_count > 0:
            isoform_terminal_exon_length = intervals_total_length(isoform_exons[-isoform_terminal_exon_count:])
            dist_to_polya = abs(isoform_exons[-isoform_terminal_exon_count-1][1] - polya_pos)
            if isoform_terminal_exon_length <= self.params.max_missed_exon_len and \
                    dist_to_polya <= self.params.max_missed_exon_len:
                logger.debug("+ Looks like missed terminal exons, shifting polyA")
                corrected_polya = isoform_exons[-1][1]
                matching_events.append(make_event(MatchEventSubtype.exon_misallignment, len(isoform_exons) - 2))

        return matching_events, corrected_polya

    # correct polyT position when fake terminal exons are present or isoform has short terminal exons
    def correct_polyt_coord(self, read_exons, fake_terminal_exon_count, polyt_pos):
        if fake_terminal_exon_count == 0:
            corrected_polyt = polyt_pos

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

    # check isoform exons prior polyT
    def check_reference_starting_exons(self, isoform_exons, polyt_pos, matching_events):
        # check isoform exons beyond polyA
        isoform_terminal_exon_count = 0
        for i in range(len(isoform_exons)):
            if isoform_exons[i][1] > polyt_pos:
                break
            isoform_terminal_exon_count += 1

        corrected_polyt = polyt_pos
        if isoform_terminal_exon_count > 0:
            isoform_terminal_exon_length = intervals_total_length(isoform_exons[:isoform_terminal_exon_count])
            dist_to_polya = abs(isoform_exons[isoform_terminal_exon_count][0] - corrected_polyt)
            if isoform_terminal_exon_length <= self.params.max_missed_exon_len and \
                    dist_to_polya <= self.params.max_missed_exon_len:
                logger.debug("+ Looks like missed terminal exons, shifting polyT")
                corrected_polyt = isoform_exons[0][0]
                matching_events.append(make_event(MatchEventSubtype.exon_misallignment, 0))

        return matching_events, corrected_polyt

