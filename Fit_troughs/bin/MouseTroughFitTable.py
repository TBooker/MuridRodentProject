from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np, pandas as pd, argparse, random
from tom import brace
from scipy.integrate import quad
from scipy.stats import gamma
import math

def I_component(i, W):
	return (-1.**i) * (1./i) *(1 - 2**(-1.*i)) * math.pow( W - 1, i)
	
def combinedSelSweepRecovery(params, dist, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination distances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	Ya = params['NeSa'] # This is 2Nes
	pa = params['pa']
	

	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth

	Ne = 426200.0
	mut_rate = 5.4e-9
	g = 0
	
	r = dist / (4*Ne)
	s = (1.*Ya) / (2*Ne)
	Va = mut_rate * pa * Ya * 2

	W = 2 * Ne * sites * Va 

	A = W/(W+(B**-2))

	S =  2 * Ne * sites * Va *  (Ya ** (-4.*r/s))

	E = (1./W) * (1./S) 

	I = math.log(2) + sum([ I_component(i, W)  for i in range(1,100)])
	
	model = (B * ( (1. - (B*W*I)) + (B*W*I-A)*(1.-E)))/(1-(A * ( 1.-E)))
	
	return model - data

def combinedSel(params, dist, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination distances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	Ya = params['NeSa'] # This is 2Nes
	pa = params['pa']

	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth


	Ne = 426200.0
	#Ne = 1000 * 0.978
	mut_rate = 5.4e-9
	#mut_rate = 0.0083/4000

	g = 0

	r = dist / (4 * Ne)
	s = (1.*Ya) / (2*Ne)
#	Va = pa
	Va = mut_rate * pa * Ya 


	S =  4 * Ne * sites * Va *  (Ya ** (-4.*r/s))
	
	model = 1./ ((1./B) + B*S)

	return model - data

def combinedTwoSpike(params, dist, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination distances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa_1 = params['NeSa'] # estimated as 4Nes in the model
	pa_1 = params['pa']
	NeSa_2 = params['NeSa_2'] # estimated as 4Nes in the model
	pa_2 = params['pa_2']


	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth

	Ne = 426200.0
#	Ne = 1000
	mut_rate = 5.4e-9
#	mut_rate = 0.0083/4000
	g = 0
	r = dist / (4*Ne)
	S = 0
	for ns, pa in zip([NeSa_1, NeSa_2],[pa_1, pa_2]):

		s = (1.*ns) / (2*Ne)
		Va = mut_rate * pa * ns * 2


		S +=  2 * Ne * sites * Va *  (ns ** (-4.*r/s))
	
	model = 1./ ((1./B) + B*S)


	return model - data


def exponentialPDF(Nes,lamb):
	return (lamb) * np.exp(-1.*lamb * Nes)


def sweep( Nes, Ner4):
	return math.pow(Nes ,-1.*Ner4/(Nes/2)) ## Rho is in terms of 4Ner

def combined(x_Nes, Ner4, lamb, pa, mut_rate): #lamb here is the mean of the exponential distribution
	Va = mut_rate * pa * x_Nes # the rate of sweeps
	return Va * sweep(x_Nes, Ner4) * exponentialPDF(x_Nes, 1/lamb)


def combinedSelExpdist(params, dist, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination distances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as Nes in the model
	pa = params['pa']
	
	ModelSums = []
	lamb = NeSa.value # The mean of the Exp. distribution	

	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth


	Ne = 426200.0
	mut_rate = 5.4e-9
	
	g = 0

	for m in dist:
		#ModelSum += [ quad(combined, 0, 10, args=(Rho * (dist + i ), NeSa)) for i in range(ExonWidth) ]
		ModelSums.append( quad(combined, 1, 1000000, args=(m, lamb, pa, mut_rate))[0] )

	Psc = (2 * Ne  * Fraction * ExonWidth * np.array(ModelSums) )

	model = 1/ ( (1/B) + B*Psc )

	return model - data


def main():
	parser = argparse.ArgumentParser(description="Combine all the sfs files coming out of the sfs_from_slim_update_bootstrap.py script")


	parser.add_argument("-i","--combined", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the file(s) that contains the combined pi data")
	parser.add_argument("-b","--bgs", 
		required = True,
		dest = "bgs",
		type =str, 
		help = "The name of the file containing the BGS information")
	parser.add_argument("--cne", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to model CNEs",
		default = False)
	parser.add_argument("--bgs_null", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to set B (as in bgs) to 1 for the calculations",
		default = False)
	parser.add_argument("--pi_0", 
		required = False,
		type = float,
		help = "Specify the neutral expectation of pi in your data [0.01]",
		default = 0.01)
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		help = "The name of the output files DONT GIVE AN EXTENSION",
		default = False)


	args = parser.parse_args()
	count = 0 

	data = pd.read_csv(args.input).sort_values('dist').dropna(axis=0, how='any')
	bgs = pd.read_csv(args.bgs).sort_values('dist').dropna(axis=0, how='any')


	bgs['B'] = bgs['Bsmooth']

	if args.cne:
		pass
	else:
		data = data[abs(data['dist']) < 2500]
		bgs = bgs[bgs['dist'] < 2500]
		
  	data = data[data['dist'] >1]
 	bgs = bgs[bgs['dist'] >1]
	
	dist = np.array(data['dist']) # mean Rho = 0.009

	if args.cne:
		combined = data['pi']/args.pi_0
	else:
		combined =  np.array(data['pi'])/args.pi_0

	B = np.array(bgs['B'])
	
	B[B > 1.] = 1. # limits all BGS values to a maximum of 1

	if args.bgs_null:
		B[B > 0.] = 1. # sets all values to 1
	


	best_models = {}
	
	for m in ['s','e', '2']:
		

		
		best_AIC = 1e6
		reps = 0		

		while reps < 10:
			params = Parameters()

			params.add('pa', value = random.random(), min = 0, max = 1.0) # Give the minimiser random seeds
			params.add('NeSa', value = random.randint(1,1000), min = 2, max = 1e8)

			if m == 'e':
			   functionToMinimize = combinedSelExpdist
			elif m =='2':
			   functionToMinimize = combinedTwoSpike
			   params.add('pa_2', value = random.random(), min = 0, max = 1) # Give the minimiser random seeds
			   params.add('NeSa_2', value = random.randint(1,1000), min = 2, max = 1e8)
			# 	#	params.add('NeSa', value = 200, vary = False) # Give the minimiser random seeds
			# 		params.add('NeSa_2', value = 9.17*2, vary = False) # Fix the starting values
			# 		params.add('pa_2', value = 0.0098, vary = False) # Fix the starting values

			elif m == 's':
			   functionToMinimize = combinedSel
   
			#  	params.add('NeSa', value = 100*2, vary = False) # Give the minimiser random seeds
			#	params.add('pa', value = 0.00016, vary = False) # Give the minimiser random seeds
			# #
			elif m == 'sr': # Sweep recovery model
			   functionToMinimize = combinedSelSweepRecovery
   
			   params.add('NeSa', value = 200*2, vary = False) # Give the minimiser random seeds
			   params.add('pa', value = 0.00016, vary = False) # Give the minimiser random seeds

			minner = Minimizer(functionToMinimize, params, fcn_args=(dist, combined, B, args.cne))

			result = minner.minimize()
			print 'Model',m,'rep',reps
			report_fit(result)

			if m == 's' or m == 'e':
				if result.covar is None: continue
				reps +=1
#				report_fit(result)
			
				rep_result = {
				'NeSa_1' :result.params['NeSa'].value,
				'pa_1' :result.params['pa'].value,
## LmFit gives the variance/covariance matrix, so we use that to get standard errors
				'NeSa_1_sd' :np.sqrt(result.covar[1][1]),
				'pa_1_sd' :np.sqrt(result.covar[0][0]),

				'NeSa_2' :-99,
				'pa_2' :-99,
				'NeSa_2_sd' :-99,
				'pa_2_sd' :-99,
				'AIC' : result.aic
				}

			elif m == '2':				
				if result.covar is None: continue
				reps +=1
				rep_result = {
				'NeSa_1' :result.params['NeSa'].value,
				'pa_1' :result.params['pa'].value,
## LmFit gives the variance/covariance matrix, so we use that to get standard errors
				'NeSa_1_sd' :np.sqrt(result.covar[1][1]),
				'pa_1_sd' :np.sqrt(result.covar[0][0]),

				'NeSa_2' :result.params['NeSa_2'].value,
				'pa_2' :result.params['pa_2'].value,
## LmFit gives the variance/covariance matrix, so we use that to get standard errors
				'NeSa_2_sd' :np.sqrt(result.covar[3][3]),
				'pa_2_sd' :np.sqrt(result.covar[2][2]),
				'AIC' : result.aic
				}		
			final = combined + result.residual

			if rep_result['AIC'] < best_AIC: 
				best_AIC = rep_result['AIC']
				best_models[m] = rep_result
				df = pd.DataFrame([dist,combined,final,B] , index = ['distance','pi','fitted', 'BGS']).transpose()
				name = args.output  + '.' + str(args.pi_0) + '.' + m + '.csv'
				df.to_csv(name, index = False)
	
	final = []
	for i in ['s', '2', 'e']:
	 	best_models[i]['model'] = i
		final.append(best_models[i])
		
	finalDF =  pd.DataFrame(final)
	print(finalDF)
	finalDF.to_csv(args.output + '.'+ str(args.pi_0) + '.ModelParams.csv', index = False)

if '__name__':
	main()

