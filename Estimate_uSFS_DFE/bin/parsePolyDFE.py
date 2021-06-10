import pandas as pd, argparse, glob
from scipy.stats import chi2
from tom import brace
import numpy as np 
from collections import OrderedDict

def parsepolyDFE(input):
	data = [i.strip() for i in open(input).readlines()]
	lnL = [i for i in data if i.startswith('---- Best joint likelihood')][0].split(' ')[5]
	results = data[data.index('--  Model: C')+1 : data.index('--  Model: C')+5]
	retDict = {}
	retDict['lnL'] = float(lnL)
 	for i,j in zip(results[2].split(), results[3].split()):
 		if i == '--': continue
 
 		if i == 'p_b':
 			retDict['pa_est'] = [float(j)]
 		elif i == 'S_b':		
 			retDict['Sb_est'] = [float(j)]
 		else:
# 		
 			retDict[i] = [float(j)]
 	retDict['alpha_div'] = [float([i for i in data if i.startswith('---- alpha_div')][0].split(' ')[-1])]
	retDict['alpha_dfe'] = [float([i for i in data if i.startswith('---- alpha_dfe')][0].split(' ')[-1])]

	return retDict
		
		
def getBootRanges(df, vanilla, taxa):
	
	dfs = []
	for i in ['Sb_est','pa_est','b','S_d','alpha_dfe','alpha_div']:
		temp = OrderedDict() 
		
		temp['taxa'] = taxa
		temp['stat'] = i
		temp['est'] = vanilla[i]
		temp['lower'] = np.percentile(df[i], 2.5)
		temp['median'] = np.median(df[i],)
		temp['upper'] = np.percentile(df[i], 97.5)
		dfs.append(pd.DataFrame(temp, columns=temp.keys()))
#)
	new_df = pd.concat(dfs)
#
# 	temp['Sb'] = vanilla['Sb_est']
# 	temp['Sb_lower'] = np.percentile(df['Sb_est'], 2.5)
# 	temp['Sb_median'] = np.median(df['Sb_est'],)
# 	temp['Sb_upper'] = np.percentile(df['Sb_est'], 97.5)
# 	
# 	temp['pa'] = vanilla['pa_est']
# 	temp['pa_lower'] = np.percentile(df['pa_est'], 2.5) 
# 	temp['pa_median'] = np.median(df['pa_est'],)
# 	temp['pa_upper'] = np.percentile(df['pa_est'], 97.5)
# 
# 	temp['b'] = vanilla['b']
# 	temp['b_lower'] = np.percentile(df['b'], 2.5)
# 	temp['b_median'] = np.median(df['b'])
# 	temp['b_upper'] = np.percentile(df['b'], 97.5)
# 
# 	temp['Sd'] = vanilla['S_d']
# 	temp['Sd_lower'] = np.percentile(df['S_d'], 2.5)
# 	temp['Sd_median'] = np.median(df['S_d'])
# 	temp['Sd_upper'] = np.percentile(df['S_d'], 97.5) 
# 
# 	temp['alpha_dfe'] = vanilla['alpha_dfe']
# 	temp['alpha_dfe_lower'] = np.percentile(df['alpha_dfe'], 2.5)
# 	temp['alpha_dfe_median'] = np.median(df['alpha_dfe'])
# 	temp['alpha_dfe_upper'] = np.percentile(df['alpha_dfe'], 97.5) 
# 
# 	temp['alpha_div'] = vanilla['alpha_div']
# 	temp['alpha_div_lower'] = np.percentile(df['alpha_div'], 2.5)
# 	temp['alpha_div_median'] = np.median(df['alpha_div'])
# 	temp['alpha_div_upper'] = np.percentile(df['alpha_div'], 97.5) 

	return new_df


def main():
	parser = argparse.ArgumentParser(description="This script makes a CSV File with the results of polyDFE in a nice table")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "Give the name of a directory contatining the multiple directories whose names end with '_boots'")
		
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the Dataframe you'll write")

	parser.add_argument("-t","--taxa", 
		required = True,
		dest = "taxa",
		type =str, 
		help = "The name of the taxa this data comes from")
	
	args = parser.parse_args()

		
	allBoots = []
	for i in glob.glob(args.input+'/*.output'):
		#	print i
			rep = i.split('/')[-1].split('.')[2]
			

			temp = parsepolyDFE(i)
			ID = i.split('/')[-1].split('.')[0]
			temp['ID'] = ID	
			if ID == 'vanilla':
				vanilla = pd.DataFrame.from_dict(temp)
			else:
				allBoots.append( pd.DataFrame.from_dict(temp) )
 	DFE = pd.concat(allBoots).sort_values('ID', ascending=1)
# 

	output = getBootRanges( DFE , vanilla, args.taxa)
#
# 	a = getBootRanges( fullDFE[ fullDFE['Divergence'] == '-'] , Sb, pa)
# 	b = getBootRanges( fullDFE[ fullDFE['Divergence'] == '+'] ,Sb, pa )
# 	c = getBootRanges( dDFE[ dDFE['Divergence'] == '-'] , Sb, pa)
# 	d = getBootRanges( dDFE[ dDFE['Divergence'] == '+'] , Sb, pa)
# 	thisDFE = pd.concat([a,b,c,d])
# 	allDFEs.append( thisDFE )
	
	output.to_csv( args.output , index = False)
	
	
if '__name__':
	main()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	