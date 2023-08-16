import pysam
import logging

from .common import overlaps, junctions_from_blocks, get_exons
from .transcript_printer import  validate_exons

logger = logging.getLogger('IsoQuant')

class VoidExonCorrector:

    def __init__(self):
        pass


    def correct_read(self, alignment_info):
        return alignment_info.read_exons


class IlluminaExonCorrector:
    
    MAX_SCORE = 1000000000000
    ABSENT_INTRON = (0,0)
    EXON_LENGTH = 50
    SIDE_DIFF = 25
    
    def __init__(self, chromosome, start, end, short_read_file):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        
        self.short_introns, self.counts = self.get_introns(short_read_file, chromosome, start, end)
        
    @classmethod
    def from_data(cls, short_introns):
        corrector = cls.__new__(cls)
        
        corrector.chromosome = 0
        corrector.start = 0
        corrector.end = 0
        corrector.short_introns = short_introns
        corrector.counts = dict()
        return corrector
        
    # get introns from a specific area of the genome
    def get_introns(self, f, chromosome, start, end):
        samfile = pysam.AlignmentFile(f, "rb") 
        intr = samfile.find_introns(samfile.fetch(chromosome, start = start, stop = end))
        samfile.close()
        i_list = set()
        for i in intr.keys():
            if(type(i)!="int"): 
                # gtf files start with 1, bam with 0, we use gtf as standard
                i_list.add((i[0]+1,i[1]))
        # in some cases counts might be necessary, so original is also saved
        return i_list, intr        

    def correct_read(self, alignment_info):
        return self.correct_exons(alignment_info.read_exons)
      
    @staticmethod
    def skipped_score(left, right, old):
        return (old[0] - left[0]) + (right[1] - old[1]) - (right[0] - left[1])
        
    @staticmethod
    def better_skipped(left, right, old, score):
        return self.skipped_score(left, right, old) < score
        
    @staticmethod
    def right_length(left, right, old):
        middle = abs(right[0] - left[1]) <= IlluminaExonCorrector.EXON_LENGTH
        left_side = abs(old[0] - left[0]) <= IlluminaExonCorrector.SIDE_DIFF
        right_side = abs(right[1] - old[1]) <= IlluminaExonCorrector.SIDE_DIFF
        return (middle and left_side and right_side)
        
    #only one side of the new introns is allowed to be the same as the old intron
    #when both the left and the right end are the same it is often a false positive
    @staticmethod
    def one_differs(left, right, old):
        return (not left[0] == old[0] or not right[1] == old[1])


    def correct_exons(self, exons):
        print("exons:", exons)
        introns = junctions_from_blocks(exons)
        corrected_introns = []
        score = IlluminaExonCorrector.MAX_SCORE
        sh = IlluminaExonCorrector.ABSENT_INTRON
        overlapping = []
        appended = False
        for i in introns:
            for s in self.short_introns:
                x = abs(i[0] - s[0]) + abs(i[1] - s[1])
                if overlaps(i, s):
                    #collect all overlapping introns to look for skipped exons
                    #still test if it is the best single match
                    overlapping.append(s)
                    if x < score:
                        score = x
                        sh = s
            #if (i[0] == sh[0] or i[1] == sh[1]) and sh[0] >= exons[0][0] and sh[1] <= exons[-1][1] and self.counts[(sh[0]-1,sh[1])] > 100:
            
            # if the best single match differs by 4 it is usually good to correct
            if ((i[0] == sh[0] and i[1] == sh[1]-4) or (i[1] == sh[1] and sh[0] == i[0]-4)):
            #if ((i[1] == sh[1]-4) or (sh[0] == i[0]-4)) and sh[0] >= exons[0][0] and sh[1] <= exons[-1][1]:
                corrected_introns.append(sh)
                appended = True
            #if nothing has been corrected yet check if there might be a skipped exon
            if len(overlapping) > 1 and not appended:
                score = IlluminaExonCorrector.MAX_SCORE
                left = IlluminaExonCorrector.ABSENT_INTRON
                right = IlluminaExonCorrector.ABSENT_INTRON
                for k in range(0, len(overlapping) -1):
                    x = overlapping[k]
                    for l in range(k, len(overlapping)):
                        y = overlapping[l]
                        # it is important that the exons are in the right order
                        # so there is two cases that might be a skipped exon
                        # It can only be a skipped exon if the gap between the introns has the right length
                        # as well as the differences between the old exon and the new exons
                        # additionally there might be multiple options so it also has to have the best score
                        if x[1] < y[0]: 
                            if self.right_length(x, y, i) and self.one_differs(x, y, i):
                                if self.better_skipped(x, y, i, score):
                                    left = x
                                    right = y
                                    score = self.skipped_score(x, y, i)
                        elif x[0] > y[1]:
                            if self.right_length(y, x, i) and self.one_differs(y, x, i):
                                if self.better_skipped(y, x, i, score):
                                    left = y
                                    right = x
                                    score = self.skipped_score(y, x, i)
                # if the left intron has been changed the requirements have been fulfilled and a correction can happen
                if not left == IlluminaExonCorrector.ABSENT_INTRON:
                    corrected_introns.append(left)
                    corrected_introns.append(right)
                    appended = True
            # if neither of the options lead to a correction the original intron is kept
            if not appended:
               corrected_introns.append(i)
            sh = IlluminaExonCorrector.ABSENT_INTRON
            score = IlluminaExonCorrector.MAX_SCORE
            overlapping = []
            appended = False
        # in case of introns in the wrong order the next steps will not work so debug option here
        if not validate_exons(get_exons((exons[0][0], exons[-1][1]), corrected_introns)):
            logger.debug("old:", introns)
            logger.debug("new:", corrected_introns)
        return get_exons((exons[0][0], exons[-1][1]), corrected_introns)
