#!/bin/bash
#SBATCH --job-name=FCPU_RS
#SBATCH --array=33-164
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=36 #4 #6 #36
#SBATCH --mem-per-cpu=5G #40G #120G #80G #45G #25G #60G #20 #30G
#SBATCH --time=168:00:00
#SBATCH --error=./ErrOut/job.%A_%a.err
#SBATCH --output=./ErrOut/job.%A_%a.out
#SBATCH --partition=#jclarke

module load anaconda
conda activate r-resamGWAS
module load plink/1.90

i=$SLURM_ARRAY_TASK_ID  #make a folder by name results and ErrOut along with plink folder

awk -v i="$i" -F "," '{print $2 "\t" $i}' /lustre/work/schnablelab/ravimural/Maize_GWAS/Phenotype_162.csv | sed '/NA/d' > pheno_${i}.txt # print genotype i-th phenotype and remove NA #Renk2021 and Zhou2021 data Included

sed 1d pheno_${i}.txt | awk '{print $1}' > genoToKeep1_${i}.txt # genotypes to keep from oryginal file

paste genoToKeep1_${i}.txt genoToKeep1_${i}.txt > genoToKeep_${i}.txt 

#ReadGenoFile
plink --bfile /lustre/work/schnablelab/ravimural/Maize_GWAS/plink/geno1051 --make-founders --keep genoToKeep_${i}.txt --maf 0.05 --make-bed --out /scratch/geno # remove genotypes with missing pheno and create bed file with maf >5%


### Kinship and PC  

#################plink --bfile /scratch/geno --out /scratch/plink --indep-pairwise 500 100 0.2 # remove SNPs in strong local LD

#################plink --bfile /scratch/geno --extract /scratch/plink.prune.in --make-bed --out /scratch/pruneGeno # make bed file for kinship and PC calculation 


R CMD BATCH -$SLURM_ARRAY_TASK_ID runrMVPFcpuRM.R runrMVPFcpuRM_${i} # run GWAS


pheno=$(head pheno_${i}.txt -n1 | awk '{print $2}')
echo $pheno Doe!
mv /scratch/$pheno.csv results/

#rm pheno.txt genoToKeep*
#done
