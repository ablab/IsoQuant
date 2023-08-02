import pysam

from .common import overlaps, junctions_from_blocks, get_exons
from .transcript_printer import  validate_exons


class VoidExonCorrector:

    def __init__(self):
        pass


    def correct_read(self, alignment_info):
        return alignment_info.read_exons


class IlluminaExonCorrector:
    
    def __init__(self, chromosome, start, end, short_read_file):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        
        self.short_introns, self.counts = self.get_introns(short_read_file, chromosome, start, end)
        
        
    def get_introns(self, f, chromosome, start, end):
        samfile = pysam.AlignmentFile(f, "rb") 
        intr = samfile.find_introns(samfile.fetch(chromosome, start = start, stop = end))
        samfile.close()
        i_list = set()
        for i in intr.keys():
            if(type(i)!="int"): 
                i_list.add((i[0]+1,i[1]))
        return i_list, intr
        

    def correct_read(self, alignment_info):
        return self.correct_exons(alignment_info.read_exons)


    def correct_exons(self, exons):
        introns = junctions_from_blocks(exons)
        corrected_introns = []
        score = 1000000000000
        sh = (0,0)
        overlapping = []
        appended = False
        for i in introns:
            for s in self.short_introns:
                x = abs(i[0] - s[0]) + abs(i[1] - s[1])
                if overlaps(i, s):
                    overlapping.append(s)
                    if x < score:
                        score = x
                        sh = s
            #if (i[0] == sh[0] or i[1] == sh[1]) and sh[0] >= exons[0][0] and sh[1] <= exons[-1][1] and self.counts[(sh[0]-1,sh[1])] > 100:
                
            if ((i[0] == sh[0] and i[1] == sh[1]-4) or (i[1] == sh[1] and sh[0] == i[0]-4)):
            #    corrected_introns.append(sh)
            #if ((i[1] == sh[1]-4) or (sh[0] == i[0]-4)) and sh[0] >= exons[0][0] and sh[1] <= exons[-1][1]:
                corrected_introns.append(sh)
                appended = True
            #else:
            #    corrected_introns.append(i)
            # if len(overlapping) > 1 and not appended:
                # score = 100
                # left = (0,0)
                # right = (0,0)
                # for k in range(0, len(overlapping) -1):
                    # x = overlapping[k]
                    # for l in range(k, len(overlapping)):
                        # y = overlapping[l]
                        # if x[1] < y[0] and (not x[0] == i[0] or not y[1] == i[1]):
                            # if abs(y[0] - x[1]) <= 50 and abs(i[0] - x[0]) <= 25 and abs(y[1] - i[1]) <= 25:
                                # if (i[0] - x[0]) + (y[1] - i[1]) - (y[0] - x[1]) < score:
                                    # left = x
                                    # right = y
                                    # score = (i[0] - x[0]) + (y[1] - i[1]) - (y[0] - x[1])
                        # elif x[0] > y[1] and (not y[0] == i[0] or not x[1] == i[1]):
                            # if abs(x[0] - y[1]) <= 50 and abs(i[0] - y[0]) <= 25 and abs(x[1] - i[1]) <= 25:
                                # if (i[0] - y[0]) + (x[1] - i[1]) - (x[0] - y[1]) < score:
                                    # left = y
                                    # right = x
                                    # score = (i[0] - y[0]) + (x[1] - i[1]) - (x[0] - y[1])
                # if not left == (0,0):
                    # corrected_introns.append(left)
                    # corrected_introns.append(right)
                    # appended = True
                    # #print("Before:", i, "After:", left, right)
            if not appended:
               corrected_introns.append(i)
            sh = (0,0)
            score = 1000000000000
            overlapping = []
            appended = False
        if not validate_exons(get_exons((exons[0][0], exons[-1][1]), corrected_introns)):
            print("old:", introns)
            print("new:", corrected_introns)
        return get_exons((exons[0][0], exons[-1][1]), corrected_introns)
