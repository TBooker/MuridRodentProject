I can generate files to be input into the mouse map
converter (http://cgd.jax.org/mousemapconverter/), which
is able to convert between mm10 bp positions and cM/Mb
recombination rates from the cox map.

The input files should be in the format:
1 1
1 1000000
1 2000000
etc, where the first column is the name of the chromosome
and the second column is the position in base pairs on that
chromosome. I can use the converter to return the sex-averaged
map for the same positions. 

Use the mm10 .fai file + the python script make_map_converter_input.py
to get one input file per chromosome.

The output files consist of the two columns of the input file and then
two columns of the output results from the conversion - chromosome again 
and cumulative recombination distance in cM/Mb (sex-averaged).

