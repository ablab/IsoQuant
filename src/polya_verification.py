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

        # assert fake_terminal_exon_count == 0 or event_to_remove == -1
        if event_to_remove != -1:
            del matching_events[event_to_remove]

        new_events, polya_pos = self.check_if_close(isoform_end, polya_info.external_polya_pos,
                                                    polya_info.internal_polya_pos, matching_events)
        if new_events is not None:
            logger.debug("Orginal polyA seems good")
            return new_events, polya_pos

        logger.debug("Orginal polyA seems distant, checking for fake terminal exons")
        external_polya_pos, internal_polya_pos = \
            self.correct_polya_positions(read_exons, fake_terminal_exon_count, polya_info)
        logger.debug("+ Corrected external %d, corrected internal %d" % (external_polya_pos, internal_polya_pos))

        polya_pos = internal_polya_pos if internal_polya_pos != -1 else external_polya_pos
        matching_events, corrected_polya_pos = \
            self.check_reference_terminal_exons(isoform_exons, polya_pos, matching_events)

        new_events, polya_pos = self.check_if_close(isoform_end, external_polya_pos, corrected_polya_pos,
                                                    matching_events)
        if new_events is not None:
            logger.debug("Corrected polyA seems good")
            return new_events, polya_pos

        dist_to_polya = abs(polya_info.external_polya_pos - isoform_end)
        logger.debug("+ Distance to polyA is %d" % dist_to_polya)
        if dist_to_polya > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(
                make_event(MatchEventSubtype.alternative_polya_site, event_length=polya_info.external_polya_pos))
        else:
            logger.debug("+ Seems like correct polyA, odd case")
            matching_events.append(
                make_event(MatchEventSubtype.correct_polya_site, event_length=polya_info.external_polya_pos))
        return matching_events, polya_info.external_polya_pos

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

        # assert fake_terminal_exon_count == 0 or event_to_remove == -1
        if event_to_remove != -1:
            del matching_events[event_to_remove]

        new_events, polyt_pos = self.check_if_close(isoform_start, polya_info.external_polyt_pos,
                                                    polya_info.internal_polyt_pos, matching_events)
        if new_events is not None:
            logger.debug("Orginal polyT seems good")
            return new_events, polyt_pos

        logger.debug("Orginal polyT seems distant, checking for fake terminal exons")
        external_polyt_pos, internal_polyt_pos = \
            self.correct_polyt_positions(read_exons, fake_terminal_exon_count, polya_info)
        logger.debug("+ Corrected external %d, corrected internal %d" % (external_polyt_pos, internal_polyt_pos))

        polyt_pos = internal_polyt_pos if internal_polyt_pos != -1 else external_polyt_pos
        matching_events, corrected_polyt_pos = \
            self.check_reference_starting_exons(isoform_exons, polyt_pos, matching_events)

        new_events, polyt_pos = self.check_if_close(isoform_start, external_polyt_pos, corrected_polyt_pos,
                                                    matching_events)
        if new_events is not None:
            logger.debug("Corrected polyT seems good")
            return new_events, polyt_pos

        dist_to_polyt = abs(polya_info.external_polyt_pos - isoform_start)
        logger.debug("+ Distance to polyT is %d" % dist_to_polyt)
        if dist_to_polyt > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(
                make_event(MatchEventSubtype.alternative_polyt_site, event_length=polya_info.external_polyt_pos))
        else:
            logger.debug("+ Seems like correct polyT, odd case")
            matching_events.append(
                make_event(MatchEventSubtype.correct_polya_site, event_length=polya_info.external_polyt_pos))
        return matching_events, polya_info.external_polyt_pos

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

    # correct polyA position when fake terminal exons are present
    def correct_polya_positions(self, read_exons, fake_terminal_exon_count, polya_info):
        assert fake_terminal_exon_count < len(read_exons)

        external_polya_pos = polya_info.external_polya_pos
        internal_polya_pos = polya_info.internal_polya_pos

        if internal_polya_pos == -1:
            return self.shift_polya(read_exons, fake_terminal_exon_count, external_polya_pos), -1

        # find exons beyond internal polyA
        polya_exon_count = 0
        for i in range(len(read_exons)):
            exon = read_exons[-i - 1]
            if exon[1] <= internal_polya_pos:
                break
            if exon[1] - internal_polya_pos > internal_polya_pos - exon[0]:
                # more than half of the exon is polyA
                polya_exon_count += 1

        extra_exons = max(fake_terminal_exon_count, polya_exon_count)
        logger.debug("Fake exons %d, polya exons %d" % (fake_terminal_exon_count, polya_exon_count))
        # correcting fake terminal exons
        return self.shift_polya(read_exons, extra_exons, external_polya_pos), \
               self.shift_polya(read_exons, extra_exons, internal_polya_pos)

    # recalculate polya site considering terminal exons are fake
    def shift_polya(self, read_exons, exon_count, polya_pos):
        if exon_count == 0 or exon_count == len(read_exons):
            return polya_pos

        dist_to_polya = 0
        for i in range(exon_count):
            exon = read_exons[-i - 1]
            if exon[0] > polya_pos:
                continue
            elif dist_to_polya == 0:
                # no exons counted yet
                dist_to_polya += polya_pos - exon[0]
            else:
                dist_to_polya += interval_len(exon)
        return read_exons[-exon_count - 1][1] + dist_to_polya

    # correct polyT position when fake terminal exons are present or isoform has short terminal exons
    def correct_polyt_positions(self, read_exons, fake_terminal_exon_count, polya_info):
        assert fake_terminal_exon_count < len(read_exons)

        external_polyt_pos = polya_info.external_polyt_pos
        internal_polyt_pos = polya_info.internal_polyt_pos

        if internal_polyt_pos == -1:
            return self.shift_polyt(read_exons, fake_terminal_exon_count, external_polyt_pos), -1

        # find exons beyond internal polyA
        polya_exon_count = 0
        for i in range(len(read_exons)):
            exon = read_exons[i]
            if exon[0] >= internal_polyt_pos:
                break
            if exon[1] - internal_polyt_pos < internal_polyt_pos - exon[0]:
                # more than half of the exon is polyT
                polya_exon_count += 1

        extra_exons = max(fake_terminal_exon_count, polya_exon_count)
        logger.debug("Fake exons %d, polya exons %d" % (fake_terminal_exon_count, polya_exon_count))
        # correcting fake terminal exons
        return self.shift_polya(read_exons, extra_exons, external_polyt_pos), \
               self.shift_polya(read_exons, extra_exons, internal_polyt_pos)

    # recalculate polya site considering terminal exons are fake
    def shift_polyt(self, read_exons, exon_count, polyt_pos):
        if exon_count == 0 or exon_count == len(read_exons):
            return polyt_pos

        dist_to_polya = 0
        for i in range(exon_count):
            exon = read_exons[i]
            if exon[1] < polyt_pos:
                continue
            elif dist_to_polya == 0:
                # no exons counted yet
                dist_to_polya += exon[1] - polyt_pos
            else:
                dist_to_polya += interval_len(exon)
        return read_exons[exon_count][0] - dist_to_polya

    # check isoform exons beyond polyA
    def check_reference_terminal_exons(self, isoform_exons, polya_pos, matching_events):
        terminal_exon_count = 0
        while terminal_exon_count < len(isoform_exons) and isoform_exons[-terminal_exon_count - 1][0] >= polya_pos:
            terminal_exon_count += 1
        if terminal_exon_count == 0 or terminal_exon_count == len(isoform_exons):
            return matching_events, polya_pos

        corrected_read_end = polya_pos
        isoform_terminal_exon_length = intervals_total_length(isoform_exons[-terminal_exon_count:])
        dist_to_polya = abs(isoform_exons[-terminal_exon_count - 1][1] - polya_pos)
        if isoform_terminal_exon_length <= self.params.max_fake_terminal_exon_len and \
                dist_to_polya <= self.params.max_fake_terminal_exon_len:
            logger.debug("+ Looks like missed terminal exons, shifting polyA")
            corrected_read_end = isoform_exons[-1][1]
            for i in range(terminal_exon_count):
                matching_events.append(make_event(MatchEventSubtype.exon_misallignment, len(isoform_exons) - 2 - i))

        return matching_events, corrected_read_end

    # check isoform exons prior polyT
    def check_reference_starting_exons(self, isoform_exons, polyt_pos, matching_events):
        # check isoform exons beyond polyA
        terminal_exon_count = 0
        while terminal_exon_count < len(isoform_exons) and isoform_exons[terminal_exon_count][1] <= polyt_pos:
            terminal_exon_count += 1
        if terminal_exon_count == 0 or terminal_exon_count == len(isoform_exons):
            return matching_events, polyt_pos

        corrected_read_start = polyt_pos
        isoform_terminal_exon_length = intervals_total_length(isoform_exons[:terminal_exon_count])
        dist_to_polya = abs(isoform_exons[terminal_exon_count][0] - polyt_pos)
        if isoform_terminal_exon_length <= self.params.max_fake_terminal_exon_len and \
                dist_to_polya <= self.params.max_fake_terminal_exon_len:
            logger.debug("+ Looks like missed terminal exons, shifting polyT")
            corrected_read_start = isoform_exons[0][0]
            for i in range(terminal_exon_count):
                matching_events.append(make_event(MatchEventSubtype.exon_misallignment, i))

        return matching_events, corrected_read_start

    def check_if_close(self, isoform_end, external_polya_pos, internal_polya_pos, matching_events):
        dist_to_external_polya = abs(isoform_end - external_polya_pos)
        dist_to_internal_polya = abs(isoform_end - internal_polya_pos)
        if min(dist_to_internal_polya, dist_to_external_polya) <= self.params.apa_delta:
            # polyA position is close, we are satisfied
            if dist_to_internal_polya < dist_to_external_polya:
                logger.debug("Internal polyAT is good %d, distance %d" %
                             (internal_polya_pos, dist_to_internal_polya))
                matching_events.append(make_event(MatchEventSubtype.correct_polya_site,
                                                  event_length=internal_polya_pos))
                return matching_events, internal_polya_pos
            else:
                logger.debug("External polyAT is good %d, distance %d" %
                             (external_polya_pos, dist_to_external_polya))
                matching_events.append(make_event(MatchEventSubtype.correct_polya_site,
                                                  event_length=external_polya_pos))
                return matching_events, external_polya_pos
        return None, -1
