#!/bin/bash

manifest_path=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/external-files/manifest_HMF_PCAWG.gene_ann.txt.gz

wd="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis"
dt=$(date +%y%m%d_%T | sed -e "s/:/-/g")


RSCRIPT=${wd}/dna-rep-ann/detGeneStatuses-bash.R


counter=0
zcat $manifest_path | tail -n +2 | while read cohort sample dir germ_vcf som_vcf cnv_tsv; do

#if [ ${counter} -lt 2 ]; then
counter=$((${counter}+1))
	
echo -e "\n######### [$counter] $sample #########"

	job_dir="${wd}/slurm_out/jobs/${dt}.genestat/"; mkdir -p $job_dir
	err_dir="${wd}/slurm_out/errs/${dt}.genestat/"; mkdir -p $err_dir
	outjob_dir="${wd}/slurm_out/outs/${dt}.genestat/"; mkdir -p $outjob_dir

	job_name=$sample
	job_script=${job_dir}/${job_name}.jobscript.txt

	if [[ ! -p ${job_script} ]]; then
        touch ${job_script};
	fi

echo "#!/bin/bash

#SBATCH --job-name=${job_name}
#SBATCH --output=${err_dir}/${job_name}.output.txt
#SBATCH --error=${outjob_dir}/${job_name}.error.txt
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=a.movasati@uu.nl

guixr load-profile /hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/.guix-profile-berner-proj --<<EOF

Rscript $RSCRIPT $dir/$germ_vcf $dir/$som_vcf $dir/$cnv_tsv $sample

EOF
" > ${job_script}


sbatch ${job_script}

#else break
#fi

done
