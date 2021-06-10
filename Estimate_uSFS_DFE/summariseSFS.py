import pandas as pd
import argparse, glob
import site_frequency_spectrum as SFS
class summarised_SFS:
	def __init__(self,sfs):
		self.pi = SFS.pi(sfs)
		self.TD = SFS.tajima(sfs)
		self.xsi = SFS.xsi(sfs)
		self.sites = sum(sfs)		
		self.div = sfs[-1] / sum(sfs)


def main():
	parser = argparse.ArgumentParser(description="Give the directory containing data for the")

	parser.add_argument("-t","--txt_dir", 
			required = True,
			dest = "txt_dir",
			type =str, 
			help = "the path to a directory containing the txt files contining site frequency spectra")
	parser.add_argument("--output", 
			required = True,
			dest = "output",
			type =str, 
			help = "The name you want to give to the output file")
	args = parser.parse_args()

	output = []
	for i in glob.glob(args.txt_dir+'/*sfs'):
		
		taxa = i.split('/')[1].split('.')[0]
		print taxa
		temp = pd.read_csv(i, header = None)
		boots = []

		for row in temp.iterrows():
			index, dataFull = row
			data = dataFull.tolist()
#    temp.append(data.tolist())
			if data[0] == 'vanilla':
				point_estimates = summarised_SFS(data[1:])
			else:
		
				boots.append( summarised_SFS(data[1:]))
		pi =  sorted([b.pi*100 for b in boots])		
		TD =  sorted([b.TD for b in boots])
		xsi = sorted([b.xsi for b in boots])
		div = sorted([b.div*100 for b in boots])		
		upper = int(len(boots)*0.975)
		lower = int(len(boots)*0.025)

		output += [[taxa,'sites',point_estimates.sites/1e6,point_estimates.sites/1e6,point_estimates.sites/1e6]]

		output += [[taxa,'pi',point_estimates.pi*100, pi[lower], pi[upper]]]
		output += [[taxa,'tajima',point_estimates.TD, TD[lower], TD[upper]]]
		output += [[taxa,'xsi',point_estimates.xsi, xsi[lower], xsi[upper]]]
		output += [[taxa,'div',point_estimates.div*100, div[lower], div[upper]]]
	out_df = pd.DataFrame(output, columns = ['taxa','statistic','value','lower','upper'])
	out_df.to_csv(args.output, index = False)
			

main()