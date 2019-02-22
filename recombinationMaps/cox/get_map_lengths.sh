# parse the ends of the per-chromosome files from the mouse map converter in
# order to make a table with the total physical and map lengths of the cox
# map for the mm10 assembly

for i in {1..19}
do
  tail -n1 chr${i}_input_mm10_1000bp_output_cMMb.txt | cut -f2,5 >> a
  echo "chr${i}" >> b
done

paste b a > lengths_per_chromosome.txt
rm a b
