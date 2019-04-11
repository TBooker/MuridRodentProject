#Download the GTF file for the mouse genic annotations

rsync -avz --progress rsync://ftp.ensembl.org/ensembl/pub/release-93/gtf/mus_musculus//Mus_musculus.GRCm38.93.gtf.gz ./

# Extract the relevant elements from the GTF
cat Mus_musculus.GRCm38.93.gtf | awk '$3 == "CDS"' >AllCDS.gtf
cat Mus_musculus.GRCm38.93.gtf | awk '$3 == "five_prime_utr"' >All5pUTR.gtf
cat Mus_musculus.GRCm38.93.gtf | awk '$3 == "three_prime_utr"' >All3pUTR.gtf

# If you're lucky, like me, and have an awesome person like Rory Craig 
# working with you, you'll have a list of the phyogenetically conserved 
# non-coding regions in the genome, all nicely sorted and merged. If
# not, I'm very sorry.
 
# Those are in the file CNEs.bed, check Rory's pipeline for info on this

# First, I need to convert each of the GTFs above into BED files
python bin/gtf2bed.py -i All5pUTR.gtf.gz -o  All5pUTR.bed --element UTR
python bin/gtf2bed.py -i All3pUTR.gtf.gz -o  All3pUTR.bed --element UTR
python bin/gtf2bed.py -i AllCDS.gtf.gz -o  AllCDS.bed --element CDS

# Sort each file:

sort -k1,1 -k2,2n All5pUTR.bed > All5pUTR.s.bed
sort -k1,1 -k2,2n All3pUTR.bed > All3pUTR.s.bed
sort -k1,1 -k2,2n AllCDS.bed > AllCDS.s.bed

# Now I need to merge overlapping elements within each type.
# I'm assuming that 3' and 5' UTRs have the same DFE so can just merge them into one file 

bedtools merge -i All5pUTR.s.bed > All5pUTR.m.bed
bedtools merge -i AllCDS.s.bed > AllCDS.m.bed
bedtools merge -i All3pUTR.s.bed > All3pUTR.m.bed

cat All5pUTR.m.bed All3pUTR.m.bed > AllUTR.bed
sort -k1,1 -k2,2n AllUTR.bed > AllUTR.s.bed
bedtools merge -i AllUTR.s.bed > AllUTR.m.bed

#bedtools merge -i CNEs.bed > CNEs.m.bed


# Now I'll subtract overlapping elements of different types
# I give priority to the different elements in the following order:

# CDS > UTR > CNE

# This means if a UTR and a CDS overlap, the CDS is retained, and so on

bedtools subtract -a AllUTR.m.bed -b AllCDS.m.bed > AllUTR.f.bed
bedtools subtract -a CNEs.stripped.bed -b AllUTR.f.bed > CNEs.f1.bed
bedtools subtract -a CNEs.f1.bed -b AllCDS.m.bed > CNEs.f2.bed

# Bedtools merge removes the 4th column, so let's add those back in:
sed -e "s/$/\tCDS/" AllCDS.m.bed > AllCDS.c.bed
sed -e "s/$/\tUTR/" AllUTRs.f.bed > AllUTRs.c.bed
sed -e "s/$/\tCNE/" CNEs.f2.bed > CNEs.c.bed

# Now let's combine all the elements into one file

cat AllCDS.c.bed AllUTRs.c.bed CNEs.c.bed > combinedElements.bed

# Sorting and merging 'combinedElements.bed' should give a file that has the exact same number of entries!

sort -k1,1 -k2,2n combinedElements.bed > combinedElements.s.bed

#There are a BUNCH of bookended elements in  this file, so if you 
# merge then things will get blended that you don't want.

# Now clean up the intermediate files

rm *.m.bed
rm *.c.bed
rm *.f1.bed
rm *.f2.bed
rm *.f.bed
rm All3pUTR.bed All5pUTR.bed AllCDS.bed All5pUTR.s.bed All3pUTR.s.bed  AllCDS.s.bed AllUTR.bed AllUTR.s.bed combinedElements.bed


# The only thing left is to BGZIP and TABIX the resulting file

~/bin/htslib-1.9/tabix combinedElements.s.bed 
~/bin/htslib-1.9/tabix combinedElements.s.bed.gz 
