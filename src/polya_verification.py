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
        read_exons = combined_read_profile.read_split_exon_profile.read_features
        read_introns = combined_read_profile.read_intron_profile.read_features
        polya_info = combined_read_profile.polya_info
        isoform_exons = self.gene_info.all_isoforms_exons[isoform_id]

        if self.gene_info.isoform_strands[isoform_id] == '+':
            if polya_info.internal_polya_pos != -1:
                # check for internal intron polyA priming
                matching_events, is_internal = self.check_internal_polya(polya_info.internal_polya_pos, matching_events)
                if is_internal:
                    return matching_events

            if polya_info.external_polya_pos != -1 or polya_info.internal_polya_pos != -1:
                matching_events, corrected_read_end = self.verify_polya(isoform_exons, read_exons,
                                                                        polya_info, matching_events)
                if corrected_read_end == -1:
                    return matching_events
                matching_events = self.remove_extra_intron_events_right(corrected_read_end, matching_events,
                                                                        read_introns)

        elif self.gene_info.isoform_strands[isoform_id] == '-':
            if polya_info.internal_polyt_pos != -1:
                # check for internal intron polyA priming
                matching_events, is_internal = self.check_internal_polyt(polya_info.internal_polyt_pos, matching_events)
                if is_internal:
                    return matching_events

            if polya_info.external_polyt_pos != -1 or polya_info.internal_polyt_pos != -1:
                matching_events, corrected_read_start = self.verify_polyt(isoform_exons, read_exons,
                                                                          polya_info, matching_events)
                if corrected_read_start == -1:
                    return matching_events
                matching_events = self.remove_extra_intron_events_left(corrected_read_start, matching_events,
                                                                        read_introns)

        if not matching_events:
            logger.warning("Empty event list after polyA verification")
            matching_events = [MatchEvent(MatchEventSubtype.none)]
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
                                                    polya_info.internal_polya_pos, matching_events,
                                                    MatchEventSubtype.correct_polya_site_right)
        if new_events is not None:
            logger.debug("Orginal polyA seems good")
            return new_events, polya_pos

        logger.debug("Orginal polyA seems distant, checking for fake terminal exons")
        matching_events, external_polya_pos, internal_polya_pos = \
            self.correct_polya_positions(read_exons, fake_terminal_exon_count, polya_info, matching_events)
        logger.debug("+ Corrected external %d, corrected internal %d" % (external_polya_pos, internal_polya_pos))

        matching_events, external_polya_pos, internal_polya_pos = \
            self.check_reference_terminal_exons(isoform_exons, external_polya_pos, internal_polya_pos, matching_events)

        new_events, polya_pos = self.check_if_close(isoform_end, external_polya_pos,
                                                    internal_polya_pos, matching_events,
                                                    MatchEventSubtype.correct_polya_site_right)
        if new_events is not None:
            logger.debug("Corrected polyA seems good")
            return new_events, polya_pos

        if polya_info.external_polya_pos == -1:
            return matching_events, -1

        polya_pos = polya_info.external_polya_pos
        dist_to_polya = abs(polya_pos - isoform_end)
        logger.debug("+ Distance to polyA is %d" % dist_to_polya)
        if dist_to_polya > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(
                MatchEvent(MatchEventSubtype.alternative_polya_site_right, event_info=polya_pos))
        else:
            logger.debug("+ Seems like correct polyA, odd case")
            matching_events.append(
                MatchEvent(MatchEventSubtype.correct_polya_site_right, event_info=polya_pos))
        return matching_events, polya_pos

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
                                                    polya_info.internal_polyt_pos, matching_events,
                                                    MatchEventSubtype.correct_polya_site_left)
        if new_events is not None:
            logger.debug("Orginal polyT seems good")
            return new_events, polyt_pos

        logger.debug("Orginal polyT seems distant, checking for fake terminal exons")
        matching_events, external_polyt_pos, internal_polyt_pos = \
            self.correct_polyt_positions(read_exons, fake_terminal_exon_count, polya_info, matching_events)
        logger.debug("+ Corrected external %d, corrected internal %d" % (external_polyt_pos, internal_polyt_pos))

        matching_events, external_polyt_pos, internal_polyt_pos = \
            self.check_reference_starting_exons(isoform_exons, external_polyt_pos, internal_polyt_pos, matching_events)

        new_events, polyt_pos = self.check_if_close(isoform_start, external_polyt_pos, internal_polyt_pos,
                                                    matching_events, MatchEventSubtype.correct_polya_site_left)
        if new_events is not None:
            logger.debug("Corrected polyT seems good")
            return new_events, polyt_pos

        if polya_info.external_polyt_pos == -1:
            return matching_events, -1

        polyt_pos = polya_info.external_polyt_pos
        dist_to_polyt = abs(polyt_pos - isoform_start)
        logger.debug("+ Distance to polyT is %d" % dist_to_polyt)
        if dist_to_polyt > self.params.apa_delta:
            logger.debug("+ Seems like APA site")
            matching_events.append(
                MatchEvent(MatchEventSubtype.alternative_polya_site_left, event_info=polyt_pos))
        else:
            logger.debug("+ Seems like correct polyT, odd case")
            matching_events.append(
                MatchEvent(MatchEventSubtype.correct_polya_site_left, event_info=polyt_pos))
        return matching_events, polyt_pos

    # remove extra introns beyond corrected read ends
    def remove_extra_intron_events_right(self, corrected_read_end, matching_events, read_introns):
        if corrected_read_end == -1:
            return matching_events

        events_to_remove = set()
        for i, event in enumerate(matching_events):
            if event.event_type in [MatchEventSubtype.extra_intron_known,
                                    MatchEventSubtype.extra_intron_novel,
                                    MatchEventSubtype.extra_intron_flanking_right,
                                    MatchEventSubtype.terminal_exon_shift_known,
                                    MatchEventSubtype.terminal_exon_shift_novel] and \
                    read_introns[event.read_region[0]][1] >= corrected_read_end:
                events_to_remove.add(i)

        for i in sorted(events_to_remove, reverse=True):
            del matching_events[i]

        return matching_events

    # remove extra introns beyond corrected read ends
    def remove_extra_intron_events_left(self, corrected_read_start, matching_events, read_introns):
        if corrected_read_start == -1:
            return matching_events

        events_to_remove = set()
        for i, event in enumerate(matching_events):
            if event.event_type in [MatchEventSubtype.extra_intron_known,
                                    MatchEventSubtype.extra_intron_novel,
                                    MatchEventSubtype.extra_intron_flanking_left,
                                    MatchEventSubtype.terminal_exon_shift_known,
                                    MatchEventSubtype.terminal_exon_shift_novel] and \
                    read_introns[event.read_region[0]][0] <= corrected_read_start:
                events_to_remove.add(i)

        for i in sorted(events_to_remove, reverse=True):
            del matching_events[i]

        return matching_events

    # check if polyA found within intron and inside mapped part of the read
    def check_internal_polya(self, internal_polya_pos, matching_events):
        if internal_polya_pos == -1:
            return matching_events, False

        for i, event in enumerate(matching_events):
            if event.event_type == MatchEventSubtype.incomplete_intron_retention_right:
                matching_events.append(MatchEvent(MatchEventSubtype.fake_polya_site_right,
                                                  isoform_region=event.isoform_region,
                                                  event_info=internal_polya_pos))
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
                matching_events.append(MatchEvent(MatchEventSubtype.fake_polya_site_left,
                                                  isoform_region=event.isoform_region,
                                                  event_info=internal_polyt_pos))
                return matching_events, True

        return matching_events, False

    # correct polyA position when fake terminal exons are present
    def correct_polya_positions(self, read_exons, fake_terminal_exon_count, polya_info, matching_events):
        assert fake_terminal_exon_count < len(read_exons)

        external_polya_pos = polya_info.external_polya_pos
        internal_polya_pos = polya_info.internal_polya_pos

        if internal_polya_pos == -1 or len(read_exons) == 0:
            return matching_events, self.shift_polya(read_exons, fake_terminal_exon_count, external_polya_pos), -1

        # find exons beyond internal polyA
        polya_exon_count = 0
        for i in range(len(read_exons)):
            exon = read_exons[-i - 1]
            if exon[1] <= internal_polya_pos:
                break
            if exon[1] - internal_polya_pos > 2 * (internal_polya_pos - exon[0]):
                # more than 2/3 of the exon is polyA
                polya_exon_count += 1

        # add events
        for i in range(fake_terminal_exon_count, polya_exon_count):
            intron_pos = len(read_exons) - i - 2
            matching_events.append(MatchEvent(MatchEventSubtype.fake_terminal_exon_right,
                                              read_region=(intron_pos, intron_pos)))

        extra_exons = max(fake_terminal_exon_count, polya_exon_count)
        logger.debug("Fake exons %d, polya exons %d" % (fake_terminal_exon_count, polya_exon_count))
        # correcting fake terminal exons
        return matching_events,\
               self.shift_polya(read_exons, extra_exons, external_polya_pos), \
               self.shift_polya(read_exons, extra_exons, internal_polya_pos)

    # recalculate polya site considering terminal exons are fake
    def shift_polya(self, read_exons, exon_count, polya_pos):
        if exon_count == 0 or exon_count == len(read_exons) or polya_pos == -1:
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
    def correct_polyt_positions(self, read_exons, fake_terminal_exon_count, polya_info, matching_events):
        assert fake_terminal_exon_count < len(read_exons)

        external_polyt_pos = polya_info.external_polyt_pos
        internal_polyt_pos = polya_info.internal_polyt_pos

        if internal_polyt_pos == -1 or len(read_exons) == 0:
            return matching_events, self.shift_polyt(read_exons, fake_terminal_exon_count, external_polyt_pos), -1

        # find exons beyond internal polyA
        polya_exon_count = 0
        for i in range(len(read_exons)):
            exon = read_exons[i]
            if exon[0] >= internal_polyt_pos:
                break
            if 2 * (exon[1] - internal_polyt_pos) < internal_polyt_pos - exon[0]:
                # more than 2/3 of the exon is polyT
                polya_exon_count += 1

        # add events
        for i in range(fake_terminal_exon_count, polya_exon_count):
            matching_events.append(MatchEvent(MatchEventSubtype.fake_terminal_exon_left, read_region=(i, i)))

        extra_exons = max(fake_terminal_exon_count, polya_exon_count)
        logger.debug("Fake exons %d, polya exons %d" % (fake_terminal_exon_count, polya_exon_count))
        # correcting fake terminal exons
        return matching_events,\
               self.shift_polyt(read_exons, extra_exons, external_polyt_pos), \
               self.shift_polyt(read_exons, extra_exons, internal_polyt_pos)

    # recalculate polya site considering terminal exons are fake
    def shift_polyt(self, read_exons, exon_count, polyt_pos):
        if exon_count == 0 or exon_count == len(read_exons) or polyt_pos == -1:
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
    def check_reference_terminal_exons(self, isoform_exons, external_polya_pos, internal_polya_pos, matching_events):
        terminal_exon_count = 0
        polya_pos = internal_polya_pos if internal_polya_pos != -1 else external_polya_pos
        while terminal_exon_count < len(isoform_exons) and isoform_exons[-terminal_exon_count-1][0] >= polya_pos:
            terminal_exon_count += 1
        if terminal_exon_count == 0 or terminal_exon_count == len(isoform_exons):
            logger.debug("+ No reference short exons")
            return matching_events, external_polya_pos, internal_polya_pos

        isoform_terminal_exon_length = intervals_total_length(isoform_exons[-terminal_exon_count:])
        dist_to_external_polya = abs(isoform_exons[-terminal_exon_count - 1][1] - external_polya_pos)
        dist_to_internal_polya = abs(isoform_exons[-terminal_exon_count - 1][1] - internal_polya_pos)
        dist_to_polya = min(dist_to_external_polya, dist_to_internal_polya)

        if (isoform_terminal_exon_length <= self.params.max_fake_terminal_exon_len and
                dist_to_polya <= self.params.max_fake_terminal_exon_len) or \
                (isoform_terminal_exon_length <= self.params.max_missed_exon_len and
                 abs(isoform_terminal_exon_length - dist_to_polya) <= self.params.delta):
            for i in range(terminal_exon_count):
                isoform_region = (len(isoform_exons) - 2 - i, len(isoform_exons) - 2 - i)
                matching_events.append(MatchEvent(MatchEventSubtype.terminal_exon_misalignment_right, isoform_region))
            corrected_read_end = isoform_exons[-1][1]
            logger.debug("+ Looks like missed terminal exons, shifting polyA: %d" % corrected_read_end)
            return matching_events, corrected_read_end, corrected_read_end

        logger.debug("+ Reference short exons present but not ok")
        return matching_events, external_polya_pos, internal_polya_pos

    # check isoform exons prior polyT
    def check_reference_starting_exons(self, isoform_exons, external_polyt_pos, internal_polyt_pos, matching_events):
        # check isoform exons beyond polyA
        terminal_exon_count = 0
        polyt_pos = internal_polyt_pos if internal_polyt_pos != -1 else external_polyt_pos
        while terminal_exon_count < len(isoform_exons) and isoform_exons[terminal_exon_count][1] <= polyt_pos:
            terminal_exon_count += 1
        if terminal_exon_count == 0 or terminal_exon_count == len(isoform_exons):
            logger.debug("+ No reference short exons")
            return matching_events, external_polyt_pos, internal_polyt_pos

        isoform_terminal_exon_length = intervals_total_length(isoform_exons[:terminal_exon_count])
        dist_to_external_polyt = abs(isoform_exons[terminal_exon_count][0] - external_polyt_pos)
        dist_to_internal_polyt = abs(isoform_exons[terminal_exon_count][0] - internal_polyt_pos)
        dist_to_polyt = min(dist_to_external_polyt, dist_to_internal_polyt)
        # if it looks like we did not align last isoform exons
        if (isoform_terminal_exon_length <= self.params.max_fake_terminal_exon_len and
                dist_to_polyt <= self.params.max_fake_terminal_exon_len) or \
                (isoform_terminal_exon_length <= self.params.max_missed_exon_len and
                 abs(isoform_terminal_exon_length - dist_to_polyt) <= self.params.delta):
            for i in range(terminal_exon_count):
                matching_events.append(MatchEvent(MatchEventSubtype.terminal_exon_misalignment_left, (i, i)))
            corrected_read_start = isoform_exons[0][0]
            logger.debug("+ Looks like missed terminal exons, shifting polyT: %d" % corrected_read_start)
            return matching_events, corrected_read_start, corrected_read_start

        logger.debug("+ Reference short exons present but not ok")
        return matching_events, external_polyt_pos, internal_polyt_pos

    def check_if_close(self, isoform_end, external_polya_pos, internal_polya_pos, matching_events, event_type):
        dist_to_external_polya = abs(isoform_end - external_polya_pos)
        dist_to_internal_polya = abs(isoform_end - internal_polya_pos)
        if min(dist_to_internal_polya, dist_to_external_polya) <= self.params.apa_delta:
            # polyA position is close, we are satisfied
            if dist_to_internal_polya < dist_to_external_polya:
                logger.debug("Internal polyAT is good %d, distance %d" %
                             (internal_polya_pos, dist_to_internal_polya))
                matching_events.append(MatchEvent(event_type, event_info=internal_polya_pos))
                return matching_events, internal_polya_pos
            else:
                logger.debug("External polyAT is good %d, distance %d" %
                             (external_polya_pos, dist_to_external_polya))
                matching_events.append(MatchEvent(event_type, event_info=external_polya_pos))
                return matching_events, external_polya_pos
        return None, -1
