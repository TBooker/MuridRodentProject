#Mmc
mkdir Mmc
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmc/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmc/chr{}.Mmc.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmc/Mmc.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 10" ::: $(seq 1 19)
cat  Mmc/chr*.out > Mmc/autosomes.out
gzip Mmc/chr*.out 

# Mmd_fra 
mkdir Mmd_fra
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmd_fra/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmd_fra/chr{}.Mmd.Fra.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_fra/Mmd_fra.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa \
--outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa \
 --min_individuals 8" ::: $(seq 1 19)
cat  Mmd_fra/chr*.out > Mmd_fra/autosomes.out
gzip Mmd_fra/chr*.out 

# Mmd_ger 
mkdir Mmd_ger
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmd_ger/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmd_ger/chr{}.Mmd.Ger.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ger/Mmd_ger.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 8" ::: $(seq 1 19)
cat  Mmd_ger/chr*.out > Mmd_ger/autosomes.out
gzip Mmd_ger/chr*.out 

# Mmd_ira 
mkdir Mmd_ira
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmd_ira/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmd_ira/chr{}.Mmd.Ira.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmd_ira/Mmd_ira.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 8" ::: $(seq 1 19)
cat  Mmd_ira/chr*.out > Mmd_ira/autosomes.out
gzip Mmd_ira/chr*.out 

# Mmm_afg 
mkdir Mmm_afg
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmm_afg/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmm_afg/chr{}.Mmm.Afg.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_afg/Mmm_afg.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 6" ::: $(seq 1 19)
cat  Mmm_afg/chr*.out > Mmm_afg/autosomes.out
gzip Mmm_afg/chr*.out 

# Mmm_cze 
mkdir Mmm_cze
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmm_cze/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmm_cze/chr{}.Mmm.Cze.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_cze/Mmm_cze.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 8" ::: $(seq 1 19)
cat  Mmm_cze/chr*.out > Mmm_cze/autosomes.out
gzip Mmm_cze/chr*.out 

# Mmm_kaz 
mkdir Mmm_kaz
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Mmm_kaz/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Mmm_kaz/chr{}.Mmm.Kaz.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mmm_kaz/Mmm_kaz.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 8" ::: $(seq 1 19)
cat  Mmm_kaz/chr*.out > Mmm_kaz/autosomes.out
gzip Mmm_kaz/chr*.out 

# Ms
mkdir Ms
parallel "python /localdisk/home/s0784966/MuridRodentTroughs/dataAnalysis/bin/bed2PK.py \
 --bed /localdisk/home/s0784966/MuridRodentTroughs/annotations/UTRs/3_prime/chr{}.bed \
 --output Ms/chr{}.out \
 --vcf /localdisk/data/troughs/VCF/Ms/chr{}.Ms.g.vcf.gz \
 --cpg /localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Ms/Ms.alleles.mm10.coords.CpGprone.fa \
 --outgroup_1 /localdisk/data/troughs/mice_synthetic_references/Mf/mm10_Mf.fa \
 --outgroup_2 /localdisk/data/troughs/mice_liftover_references/mm10_pahari_alleles.fasta --outgroup_cpg_1 \
/localdisk/data/troughs/mice_synthetic_references_cpgprone_info/Mf/Mf.alleles.mm10.coords.CpGprone.fa --outgroup_cpg_2 \
/localdisk/data/troughs/mice_liftover_references_cpgprone_info/Mp.alleles.mm10.coords.CpGprone.fa --min_individuals 8" ::: $(seq 1 19)
cat  Ms/chr*.out > Ms/autosomes.out
gzip Ms/chr*.out 

