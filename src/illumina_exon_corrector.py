import pysam

from .common import overlaps
from .short_utils import get_introns

class IlluminaExonCorrector:
	
	def __init__(self, chromosome, start, end, short_read_file):
		self.chromosome = chromosome
		self.start = start
		self.end = end
		
		self.short_introns = get_introns(short_read_file, chromosome, start, end)
		
		
	def correct_read(self, introns):
		corrected_introns = set()
		score = 1000000000000
		sh = (0,0)
		for i in introns:
			for s in self.short_introns:
				x = abs(i[0] - s[0]) + abs(i[1] - s[1])
				if (overlaps(i, s) and abs(i[0] - s[0]) <= 12 and abs(i[1] - s[1]) <= 12):
					if x < score:
						score = x
						sh = s
			if (i[0] == sh[0] and i[1] == sh[1]-4) or (i[1] == sh[1] and sh[0] == i[0]-4):
				corrected_introns.add(sh)
			else:
				corrected_introns.add(i)
			sh = (0,0)
			score = 1000000000000
		return corrected_introns
