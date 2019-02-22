"""
Make an input file for the mouse map converter with mm10 chromosomes + bp
coordinates. Specify the size of the interval between positions on the command
line
"""

import argparse

import argparse

parser = argparse.ArgumentParser(description="""Make an input file for the mouse map converter with mm10 chromosomes""")
parser.add_argument('gap', help='the size of the gap between each position in the input file')
args = parser.parse_args()

gap = int(args.gap)

length_dict = {}
with open('/home/ben/data/mice_reference/reference_UCSC/mm10.fa.fai', 'rt') as f:
    for line in f:
        l = line.strip().split()
        length_dict[l[0]] = int(l[1])

for i in range(1, 20):
    with open('chr' + str(i) + '_input_mm10_' + str(gap) + 'bp.txt', 'wt') as f:
        for j in range(0, length_dict['chr' + str(i)], gap):
            if j == 0:
                f.write(str(i) + '\t1\n')
            else:
                f.write(str(i) + '\t' + str(j) + '\n')
