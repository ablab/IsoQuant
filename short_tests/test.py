import sys
import re


args = sys.argv

count = 0

f = open(args[1], 'r')
for l in f:
	chunks = l.split(' ')
	x1 = chunks[0]
	x2 = chunks[1]
	y1 = chunks[3]
	y2 = chunks[4]
	x1 = int(re.sub('\D', '', x1))
	x2 = int(re.sub('\D', '', x2))
	y1 = int(re.sub('\D', '', y1))
	y2 = int(re.sub('\D', '', y2))
	if (x1 == y1 and x2 == y2-4):
		#print(x1,x2,y1,y2)
		count = count + 1
	elif (x2 == y2 and y1 == x1-4):
		count = count + 1
		#print(x1,x2,y1,y2)

		
		
print(count)
