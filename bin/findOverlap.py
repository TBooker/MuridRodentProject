import sys

lastEnd = 0
for i in open(sys.argv[1]):
	x = i.strip().split()
	start = int(x[1])
	end = int(x[2])
	if start - lastEnd <=0:
		print x
	lastEnd = end