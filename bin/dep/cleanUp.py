import glob


notPresent = []
surplus = []
runs = []
for i in glob.glob('*txt'):
	runs.append(int(i.split('.')[0].split('n')[1]))

for i in range(len(runs) ):
	if i+1 not in runs:
		notPresent.append(i)

previous = 0		
for i in runs:
	if i != previous + 1:
		surplus.append(i)
	previous = i

for i,j in zip(notPresent, surplus):
	print i, j
