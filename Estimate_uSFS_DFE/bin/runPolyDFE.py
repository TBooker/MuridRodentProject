import pandas as pd
import argparse, glob
import site_frequency_spectrum as SFS


def polyDFEconfig(sel,neu):
        outString = '1\t1\t'+str(len(sel)-1)+'\n' 
        outString += '\t'.join(map(str, neu[1:-1] + [sum(neu) , neu[-1], sum(neu) ])) + '\n'
        outString += '\t'.join(map(str, sel[1:-1] + [sum(sel) , sel[-1], sum(sel) ]))
        return outString


#curl -F"d=593" -F"d0=930" -F"xlow=0.1" -F"xhigh=0.9" -F"datafile=@polymorphisms.txt" -o "MK_full.html" http://benhaller.com/cgi-bin/R/asymptoticMK_run.html

def main():
	parser = argparse.ArgumentParser(description="Give the directory containing data for the")

	parser.add_argument("-n","--neu_txt", 
			required = True,
			dest = "neu_txt",
			type =str, 
			help = "the path to a directory containing the txt files contining neutral site frequency spectra")
	parser.add_argument("-s","--sel_txt", 
			required = True,
			dest = "sel_txt",
			type =str, 
			help = "the path to a directory containing the txt files contining selected site frequency spectra")
	parser.add_argument("-o","--output_dir", 
			required = True,
			dest = "output_dir",
			type =str, 
			help = "The name of the directory you want to dump the output file to")
	args = parser.parse_args()

	output = []
	
	for i, j in zip(open(args.neu_txt), open(args.sel_txt)):
		
		neu = i.strip().split(',')
		sel = j.strip().split(',')

		neuSFS = map(float, neu[1:])
		selSFS = map(float, sel[1:])

		polyDFE =  polyDFEconfig(selSFS,neuSFS)
		name = args.output_dir + '/' + neu[0]+'.polyDFE.config.txt'
		temp = open(name,'w')		
		temp.write(polyDFE)
		temp.close()
# 		taxa = i.split('/')[1].split('.')[0]
# 		print taxa
# 		temp = pd.read_csv(i, header = None)
# 		boots = []
# 
# 		for row in temp.iterrows():
# 			index, dataFull = row
# 			data = dataFull.tolist()
# #    temp.append(data.tolist())
# 			if data[0] == 'vanilla':
# 				point_estimates = summarised_SFS(data[1:])
# 			else:
# 		
# 				boots.append( summarised_SFS(data[1:]))
# 		pi =  sorted([b.pi*100 for b in boots])		
# 		TD =  sorted([b.TD for b in boots])
# 		xsi = sorted([b.xsi for b in boots])
# 		div = sorted([b.div*100 for b in boots])		
# 		upper = int(len(boots)*0.975)
# 		lower = int(len(boots)*0.025)
# 
# 		output += [[taxa,'sites',point_estimates.sites/1e6,point_estimates.sites/1e6,point_estimates.sites/1e6]]
# 
# 		output += [[taxa,'pi',point_estimates.pi*100, pi[lower], pi[upper]]]
# 		output += [[taxa,'tajima',point_estimates.TD, TD[lower], TD[upper]]]
# 		output += [[taxa,'xsi',point_estimates.xsi, xsi[lower], xsi[upper]]]
# 		output += [[taxa,'div',point_estimates.div*100, div[lower], div[upper]]]
# 	out_df = pd.DataFrame(output, columns = ['taxa','statistic','value','lower','upper'])
# 	out_df.to_csv(args.output, index = False)
# 			

main()