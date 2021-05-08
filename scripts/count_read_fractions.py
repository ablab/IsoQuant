import numpy
import sys
from collections import defaultdict

read_info_map = {}
print("Reading read fractions")
for l in open(sys.argv[1]):
    t = l.strip().split()
    read_id = t[0]
    intron_count = t[7].count(",")
    additional_info = t[8]
    if not additional_info.startswith("read_fraction="):
        continue
    read_fraction = float(additional_info[len("read_fraction="):-1])
    read_info_map[read_id] = (intron_count, read_fraction)

print("Reading read assignments")
read_fraction_arrays = defaultdict(list)
for l in open(sys.argv[2]):
    t = l.strip().split()
    read_id = t[0]
    assignment_type = t[5]
    if assignment_type != "inconsistent" or read_id not in read_info_map:
        continue

    read_info = read_info_map[read_id]
    read_fraction_arrays[read_info[0]].append(read_info[1])

print("Calculating histograms")
for intron_count in sorted(read_fraction_arrays.keys()):
    if intron_count <= 3:
        continue
    read_fraction_arrays[3] += read_fraction_arrays[intron_count]

bins = [0.1 * i for i in range(11)]
print("Histogram for bin 0.1")
for intron_count in range(1,4):
    hist = numpy.histogram(read_fraction_arrays[intron_count], bins, density=True)[0]
    print("%d\t%s" % (intron_count, "\t".join(["%.4f" % (x * 10) for x in hist])))

bins = [0.2 * i for i in range(6)]
print("Histogram for bin 0.2")
for intron_count in range(1,4):
    hist = numpy.histogram(read_fraction_arrays[intron_count], bins, density=True)[0]
    print("%d\t%s" % (intron_count, "\t".join(["%.4f" % (x * 20) for x in hist])))


