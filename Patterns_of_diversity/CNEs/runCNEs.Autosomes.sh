
mkdir Mmc
for i in $(seq 1 19)
do
parallel -j 60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmc/chr$i.Mmc.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmc/Mmc.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --min_individuals 10 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmc/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done
# #/localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/
# #/localdisk/data/troughs/VCF/Mmc/
exit 0

mkdir Mmd_ger
for i in $(seq 1 19)
do
mkdir Mmc
for i in $(seq 1 19)
do
parallel -j 60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_ger/chr$i.Mmd.Ger.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ger/Mmd_ger.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmd_ger/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_ger/chr$i.Mmd.Ger.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ger/Mmd_ger.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmd_ger/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done
# 
# mkdir Mmd_ger
# for i in $(seq 1 19)
# do
# parallel "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_ger/chr$i.Mmd.Ger.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ger/Mmd_ger.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 8    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed  --output Mmd_ger/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done
# #exit 0

mkdir Mmd_fra
for i in $(seq 1 19)
do
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_fra/chr$i.Mmd.Fra.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_fra/Mmd_fra.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmd_fra/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done
# mkdir Mmd_fra
# for i in $(seq 1 19)
# do
# parallel "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_fra/chr$i.Mmd.Fra.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_fra/Mmd_fra.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 8    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed  --output Mmd_fra/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done

mkdir Mmd_ira
for i in $(seq 1 19)
do
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_ira/chr$i.Mmd.Ira.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ira/Mmd_ira.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmd_ira/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done
# 
# mkdir Mmd_ira
# for i in $(seq 1 19)
# do
# parallel "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmd_ira/chr$i.Mmd.Ira.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ira/Mmd_ira.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 8    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed  --output Mmd_ira/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done


mkdir Mmm_afg
for i in $(seq 1 19)
do
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmm_afg/chr$i.Mmm.Afg.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_afg/Mmm_afg.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 6 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmm_afg/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done

# mkdir Mmm_afg
# for i in $(seq 1 19)
# do
# parallel  "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmm_afg/chr$i.Mmm.Afg.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_afg/Mmm_afg.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 6    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmm_afg/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done

mkdir Mmm_cze
for i in $(seq 1 19)
do
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmm_cze/chr$i.Mmm.Cze.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_cze/Mmm_cze.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmm_cze/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done

# mkdir Mmm_cze
# for i in $(seq 1 19)
# do
# parallel  "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmm_cze/chr$i.Mmm.Cze.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_cze/Mmm_cze.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 8    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed  --output Mmm_cze/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done


mkdir Mmm_kaz
for i in $(seq 1 19)
do
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmm_kaz/chr$i.Mmm.Kaz.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_kaz/Mmm_kaz.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Mmm_kaz/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done

# mkdir Mmm_kaz
# for i in $(seq 1 19)
# do
# parallel "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Mmm_kaz/chr$i.Mmm.Kaz.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_kaz/Mmm_kaz.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 8    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed  --output Mmm_kaz/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done

mkdir Ms
for i in $(seq 1 19)
do
parallel -j60 "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Ms/chr$i.Ms.g.vcf.gz  --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Ms/Ms.alleles.mm10.coords.CpGprone.fa --outgroup_1 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_Rn6_alleles.fa --min_individuals 8 --outgroup_cpg_1 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Rn.alleles.mm10.coords.CpGprone.fa --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed --output Ms/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/)
done
# 
# mkdir Ms
# for i in $(seq 1 19)
# do
# parallel "python ../bin/vcf2sfs.py --vcf /localdisk/data/troughs/VCF/Ms/chr$i.Ms.g.vcf.gz --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Ms/Ms.alleles.mm10.coords.CpGprone.fa  --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa   --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_caroli_alleles.fasta --min_individuals 8    --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa  --outgroup_cpg_2 /localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mc.alleles.mm10.coords.CpGprone.fa     --bed /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/CNEs/w.100_10000/bed.s/{}/chr$i.bed  --output Ms/{}.chr$i.txt" ::: $(ls /localdisk/home/s0784966/MuridRodentTroughs/mm10_analysis_windows/exons/w.1000/bed.s/)
# done
