import argparse, random, subprocess, glob, math
import numpy as np
import pandas as pd
import re


def pi(SFS,per_site = True):
	if sum(SFS) ==0: 
		return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	pi = sum([(1.0*i*(N-i)*(SFS[i]))/(binom) for i in xrange(N) if i != 0])
	if per_site == True:
		return pi/sum(SFS)
	else:
		return pi


def xsi(SFS):
	if len(SFS) ==0: return -99
#        print SFS
	N = len(SFS)-1
	total = sum(SFS[1:-1])
	xsi_numerator = sum([float(i*i*SFS[i])/total for i in xrange(N) if i > N/2 and i !=N])
	xsi_denominator = sum([float(i*SFS[i])/total for i in xrange(N) if i > N/2 and i !=N])
	if xsi_denominator == 0:
		return -99
	else:
		return (xsi_numerator/xsi_denominator)/N

def fwh(SFS,per_site = True):
	if len(SFS) ==0: return -99
	N = len(SFS)-1
	binom = (N * (N -1))/2
	fwh = sum([(1.0*i*i*(SFS[i]))/binom for i in xrange(N) if i != 0 or  i != N] )
	if per_site == True:
		return fwh/sum(SFS)
	else:
		return fwh


def theta_W(SFS,per_site=True):
	if len(SFS) ==0: return -99
	N = len(SFS)-1
	S = sum(SFS[1:-1]) ## this slice takes the interior of the SFS, gets S
	harmonic = sum(1.0/d for d in xrange(1, N-1,1))
	if per_site == True:
		return float(S)/(harmonic*sum(SFS))
	else:
		return float(S)/(harmonic)

def tajima(SFS):
	#if len(SFS) <2: return None
	th_pi = pi(SFS,per_site=False)
	N = len(SFS)-1
	S = sum(SFS[1:-1])
	length = sum(SFS)
	a1 = sum([1.0/i for i in range(1, N)])
	a2 = sum([1.0/(i**2) for i in range(1, N)])
	b1 = float(N+1)/(3*(N-1))
	b2 = float(2 * ( (N**2) + N + 3 )) / (9*N*(N-1))## Tajima 1989 Equation 9
	c1 = b1 - 1.0/a1
	c2 = b2 - float(N+2)/(a1 * N) + float(a2)/(a1 ** 2)
	e1 = float(c1) / a1
	e2 = float(c2) / ( (a1**2) + a2 )
	if ((e1 * S )+ ((e2 * S) * (S - 1))) == 0.0:   # Tajima 1989 Equation 38
		return -99
	else: 
		D = (float(th_pi - (float(S)/(a1)))/math.sqrt((e1 * S )+ ((e2 * S) * (S - 1) )))
		return D

def parseLabel(label):
	d = label.split(':')[0]
	if d == 'u':
		mult = -1.
	elif d == 'd':
		mult = 1.
	dist = sum(map(float,label.split(':')[1].split('-')))/2
	return dist * mult

def bootstrap(array):
	temp = array[np.random.randint(array.shape[0], size=array.shape[0]), :]
	return np.sum(temp, axis=0)

def SFSsummary(sfs):
	return np.array([pi(sfs), theta_W(sfs), tajima(sfs), xsi(sfs),fwh(sfs)])

		
def main():
	parser = argparse.ArgumentParser(description="Gets the SFS files from a directory and collates the site frequency spectra from a directory containing many entries with the suffix *sfs.txt.gz")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str, 
			help = "name the directory containing the site frequnecy spectra files")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str, 
			help = "name the output file")
	parser.add_argument("--boots", 
			required = False,
			dest = "boots",
			type = int,
			default = 1000, 
			help = "the number of bootstrap replicates that you want to perform [1000]")

	args = parser.parse_args()
	data = {}
	for i in glob.glob(args.input +'/*sfs.txt.gz'):
		try:
			for k in pd.read_csv(i, compression = 'gzip').values:
				if True in pd.isna(k[2:]):
					continue
				try:
					data[k[0]].append(k[2:])
				except KeyError:
					data[k[0]] = [k[2:]]
		except pd.errors.EmptyDataError:
			continue
			
	output = []

#	print data['d:0-1000']
	for key in sorted(data.keys()):
		data[key] = np.array(data[key])
		full_sfs = np.sum(data[key], axis=0)
		summary = np.empty([100, 5 ])
		for i in range(100):
			x =  bootstrap(data[key])
			summary[i] = SFSsummary(x)

		stats = [pi, theta_W, tajima, xsi,fwh]
		outline = [key, parseLabel(key), full_sfs.sum()]
		for k in range(len(stats)):
			outline += [stats[k](full_sfs)] + list(np.percentile(summary.T[k], [2.5,97.5], interpolation = 'nearest'))
			
		output.append(outline)
	names = ['label','dist', 'sites']
	for i in ['pi', 'theta_W', 'tajima', 'xsi','fwh']:
		names.append(i)
		names.append(i+'_lower')
		names.append(i+'_upper')
		
	pd.DataFrame(output, columns = names).to_csv(args.output,index = False)

main()
