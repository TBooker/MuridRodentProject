import glob, argparse
import numpy as np 
import scipy.integrate
from scipy.stats import gamma

## Produce a "frozen" form of the gamma distribution of interest. This will give the distribution with the parameters of interest and allow us to integrate over it easily
## The parameterisation order is: shape	and scale = mean*shape
def prop_muts_in_range(shape, mean):
	rv = gamma(shape, scale = mean/shape)

	nes_0_1 = scipy.integrate.quad(rv.pdf, 0, 1)[0]
	nes_1_10 = scipy.integrate.quad(rv.pdf, 1, 10)[0]
	nes_10_inf = scipy.integrate.quad(rv.pdf, 10, np.inf)[0]
	
	return [nes_0_1, nes_1_10, nes_10_inf]

def parsepolyDFE(input):
	data = [i.strip() for i in open(input).readlines()]
	lnL = [i for i in data if i.startswith('---- Best joint likelihood')][0].split(' ')[5]
	results = data[data.index('--  Model: C')+1 : data.index('--  Model: C')+5]
	retDict = {}
	
	retDict['lnL'] = lnL
	for i,j in zip(results[2].split(), results[3].split()):
		if i == '--': continue

		if i == 'p_b':
			retDict['pa_est'] = [float(j)]
		elif i == 'S_b':		
			retDict['Sb_est'] = [float(j)]
		else:
		
			retDict[i] = [float(j)]

	retDict['product_est'] = [retDict['pa_est'][0] * retDict['Sb_est'][0]]
	return retDict

def main():
	parser = argparse.ArgumentParser(description="It contains a python implementation of Peter's prop_muts_in_range program. It gives the proportion of mutations in the different categories of selection strength.")
	parser.add_argument("--input", 
			required = True,
			dest = "input",
			type = str, 
			help = "the directory containing the polyDFE output files")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type = str, 
			help = "Name the output file")
	parser.add_argument("--taxa", 
			required = True,
			dest = "taxa",
			type = str, 
			help = "Name the taxa")
	parser.add_argument("--type", 
			required = True,
			dest = "type",
			type = str, 
			help = "Name the site type")

	args = parser.parse_args()
	
	results = []
	count = 0
	for i in glob.glob(args.input + '/*polyDFE.dDFE.output'):
		rep =  i.split('/')[-1].split('.')[0]
		polyDFE_results = parsepolyDFE(i)
		bins = prop_muts_in_range( float(polyDFE_results['b'][0]), abs(float(polyDFE_results['S_d'][0])) ) 
		if rep == 'vanilla':
			point = bins
		else:
			count +=1
			results.append( np.array(bins) )
#		if count == 10: break

	boots = np.array(results)

	lower = np.percentile(boots.T, q = 2.5, axis=1)
	upper = np.percentile(boots.T, q = 97.5, axis=1)

	output = open(args.output, 'w')
	
	output.write( ','.join([args.taxa, args.type, '0-1', str(point[0]), str(lower[0]), str(upper[0]) ]) + '\n' )
	output.write( ','.join([args.taxa, args.type, '1-10',  str(point[1]), str(lower[1]), str(upper[1]) ]) + '\n' )
	output.write( ','.join([args.taxa, args.type, '10-Inf',  str(point[2]), str(lower[2]), str(upper[2]) ]) + '\n' )

	output.close()

main()

#print(rv.rvs(100))
