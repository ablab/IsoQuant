
############################################################################
# Copyright (c) 2022-2026 University of Helsinki
# All Rights Reserved
# See file LICENSE for details.
############################################################################

from enum import Enum, unique

from .common import junctions_from_blocks, overlaps


@unique
class Stats(Enum):
    true_negative = 1  #correct before and not changed
    false_negative = 2 #either changed but not correct after or not correct before and not changed
    true_positive = 3  #changed and it is correct after
    false_positive = 4 #changed and it is not correct after but was correct before


class CorrectionStats:
    
    def __init__(self, reference_db, long_reads):
        
        self.reference = reference_db
        self.long_reads = long_reads
        
        self.name_map = self.map_names(self.reference, self.long_reads)
    
    def map_names(self, ref, long_reads):
        transcripts = list(ref.features_of_type('transcript', order_by=('seqid', 'start')))
        names = {}
        for alignment in long_reads.fetch():
            name = alignment.query_name
            name = name.split('_')[0]
            if not name in names.keys():
                for t in transcripts:
                    if name in t["transcript_id"][0]:
                        names[name] = t
                        break
        return names
    
    def get_introns_from_transcript(self, transcript):
        exons = []
        for e in self.reference.children(transcript, order_by='start'):
            if e.featuretype == 'exon':
                    exons.append((e.start, e.end))
        introns= junctions_from_blocks(exons)
        introns = set(introns)
        return introns
        
    def read_stats(self, introns, corrected_introns, alignment):
        name = alignment.query_name
        name = name.split('_')[0]
        transcript = self.name_map[name]
        reference_introns = self.get_introns_from_transcript(transcript)
        n_ref = len(reference_introns)
        n_before = len(introns.intersection(reference_introns))
        n_after = len(corrected_introns.intersection(reference_introns))
        changed = not(introns == corrected_introns)
        correct_before = (n_before == len(introns))
        correct_after = (n_after == len(corrected_introns))
        if changed:
            if correct_after:
                return Stats.true_positive
            elif not correct_before:
                #if len(introns.intersection(reference_introns)) > len(corrected_introns.intersection(reference_introns)):
                #    print("False Positive in Read", alignment.query_name)
                #    print("Before:", sorted(introns))
                #    print("After:", sorted(corrected_introns))
                #    return Stats.false_positive
                #else:
                print("Change but wrong before and after in Read", alignment.query_name, "Correct before:", len(introns.intersection(reference_introns)), "Correct after:", len(corrected_introns.intersection(reference_introns)))
                return Stats.false_negative
            else:
                print("False Positive in Read", alignment.query_name)
                print("Before:", sorted(introns))
                print("After:", sorted(corrected_introns))
                return Stats.false_positive
        else:
            if correct_before:
                return Stats.true_negative
            elif n_before != 0:
                print("No correction but wrong Introns in Read", alignment.query_name)
                return Stats.false_negative
            else:
                return Stats.true_negative
                
    def stats_single(self, before, after, reference_introns):
        if before in reference_introns:
            print("False positive, before:", before, "after:", after)
            return Stats.false_positive
        elif after in reference_introns:
            print("True positive, before:", before, "after:", after)
            return Stats.true_positive
        else:
            print("False negative with change, before:", before, "after:", after)
            return Stats.false_negative
        
    def intron_stats(self, introns, corrected_introns, alignment):
        classification = []
        name = alignment.query_name
        name = name.split('_')[0]
        transcript = self.name_map[name]
        reference_introns = self.get_introns_from_transcript(transcript)
        n_ref = len(reference_introns)
        changed = not(introns == corrected_introns)
        unchanged_introns = [intron for intron in introns if intron in corrected_introns]
        #unchanged_introns = introns.intersection(corrected_introns)
        for intron in unchanged_introns:
            if intron in reference_introns:
                print("True negative", intron)
                classification.append(Stats.true_negative)
            else:
                classification.append(Stats.false_negative)
        #diff_before = list(sorted(introns - unchanged_introns))
        #diff_after = list(sorted(corrected_introns - unchanged_introns))
        diff_before = [intron for intron in introns if not intron in unchanged_introns]
        diff_after = [intron for intron in corrected_introns if not intron in unchanged_introns]
        if len(diff_before) == len(diff_after):
            for i in range(len(diff_before)):
                classification.append(self.stats_single(diff_before[i], diff_after[i], reference_introns))
        elif len(diff_before) < len(diff_after):
            j = 0
            for i in range(len(diff_before)):
                b = diff_before[i]
                a = diff_after[i+j]
                if b[0] == a[0] and b[1] == a[1]-4:
                    classification.append(self.stats_single(b, a, reference_introns))
                elif b[1] == a[1] and a[0] == b[0]-4:
                    classification.append(self.stats_single(b, a, reference_introns))
                else:
                    print(len(diff_before))
                    print(len(diff_after))
                    print(i)
                    print(j)
                    right = diff_after[i+j+1]
                    if overlaps(a, b) and overlaps(b, right):
                        j = j + 1
                        if b in reference_introns:
                            classification.append(Stats.false_positive)
                            print("False positive, before:", b, "after:", a, right)
                        elif a in reference_introns and right in reference_introns:
                            classification.append(Stats.true_positive)
                            print("True positive, before:", b, "after:", a, right)
                        else:
                            classification.append(Stats.false_negative)
                            print("False negative with change, before:", b, "after:", a, right)
                    elif overlaps(a, b):
                        classification.append(self.stats_single(b, a, reference_introns))
                    elif overlaps(b, right):
                        classification.append(self.stats_single(b, right, reference_introns))
        else:
            print("more exons before than after")
            print(diff_before)
            print(diff_after)
        return classification
