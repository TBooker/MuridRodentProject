import glob, subprocess

for i in glob.glob('*neq.trees'):
	x = i.split('.')[0]
	command = ['mv',  x +'.*' , 'done/']
	print ' '.join(command)
