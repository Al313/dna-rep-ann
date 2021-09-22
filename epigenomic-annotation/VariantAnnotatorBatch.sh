#!/bin/bash


wd="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation"
dt=$(date +%y%m%d_%T | sed -e "s/:/-/g")
RSCRIPT=${wd}/variantAnnotator.R


tadBedFileInput="$wd/annotation-beds/chromatin-status-tads.bed.bgz"  #args[2]
binBedFileInput="$wd/annotation-beds/chromatin-status-bins.bed.bgz"  #args[3]
repTimingFileInput="$wd/annotation-beds/repliseq_encode_mean_binned.bed.bgz"  #args[4]
repTimingFileInputBoxtel="$wd/annotation-beds/boxtel_paper_all_RepliSeq_median.bed.bgz" #args[5]
cpgInputFile="$wd/annotation-beds/cpgIslandExtUnmasked.bed.bgz"  #args[6]
repOrientationInputFile="$wd/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz"  #args[7]

job_info="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/slurm_out"


for i in {1..7049}; do



job_dir="${job_info}/jobs/${dt}.genestat"; mkdir -p $job_dir
err_dir="${job_info}/errs/${dt}.genestat"; mkdir -p $err_dir
outjob_dir="${job_info}/outs/${dt}.genestat"; mkdir -p $outjob_dir


job_script=${job_dir}/${i}.jobscript.txt

        if [[ ! -p ${job_script} ]]; then
        touch ${job_script};
        fi

echo "#!/bin/bash

#SBATCH --job-name=${i}
#SBATCH --output=${outjob_dir}/${i}.output.txt
#SBATCH --error=${err_dir}/${i}.error.txt
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=a.movasati@uu.nl

guixr load-profile /hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/.guix-profile-berner-proj --<<EOF

Rscript $RSCRIPT ${i} ${tadBedFileInput} ${binBedFileInput} ${repTimingFileInput} ${repTimingFileInputBoxtel} ${cpgInputFile} ${repOrientationInputFile} 

EOF
" > ${job_script}

sbatch ${job_script}

done




