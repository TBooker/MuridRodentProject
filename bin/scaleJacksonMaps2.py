import pandas as pd 
import numpy as np
from scipy import stats
from scipy.stats.mstats import gmean, hmean


def geo_mean_overflow(iterable):
    a = np.log(iterable)
    return np.exp(a.sum()/len(a))

## This function works with the recombination map encoding used by Liu et al 2014. Check their .csv files to see how others should be formatted
def getRecomIntervalPedigree(interval):
	map = pd.read_csv('/home/booker/work/MouseAdap/recombinationMaps/Liu/avg.map.csv') 
	chrom = int(interval[0].split('r')[1])
	cmap = map[map['Chr'] == chrom].copy().reset_index()
	new = [] 
	length = interval[2] - interval[1]
	prev_cM = 0 
	for index, row in cmap.iterrows(): 
		if prev_cM != row.cM: 
			new.append(row) 
		prev_cM = row.cM    
	newDF = pd.DataFrame(new) 
	newDF['cM_diff'] = newDF['cM'].diff()  
	newDF['bp_diff'] = newDF['Pos'].diff()/1e6
	newDF['cM_Mb'] = newDF['cM_diff']/newDF['bp_diff'] ## cM/Mb values for each interval

	# Grab the chunk of the dataframe that corresponds to the simulated interval
	slicey = newDF[(newDF['Pos'] > interval[1]) & (newDF['Pos']  < interval[2])]
	# If the region you're simulating doesn't contain any recorded rates, use the nearest values
	if slicey.shape[0] ==0:
		newDF['pos_diff'] = newDF['Pos'] - region[2]
		newDF['pos_diff'][newDF['pos_diff'] < 0] = np.nan
		rate = newDF.loc[newDF['pos_diff'].idxmin()].cM_Mb		#print newDF
		return [rate], [length-1]
	else:
		lastRate = newDF.loc[slicey.index.max()+1]
		return list(slicey['cM_Mb'])+[lastRate.cM_Mb], list(slicey['Pos'] - interval[1]) + [length-1] 
	
region = ('chr19',7000000,8000000)
print getRecomIntervalPedigree(region)
