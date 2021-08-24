#!/bin/bash



wd="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis"
dt=$(date +%y%m%d_%T | sed -e "s/:/-/g")
RSCRIPT=${wd}/dna-rep-ann/combineTimingMs.R

start_pos=(1 200 400 600 800 1000 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000 4200 4400 4600 4800 5000 5200 5400 5600 5800 6000 6200 6400 6600 6800 7000)

end_pos=(199 399 599 799 999 1199 1399 1599 1799 1999 2199 2399 2599 2799 2999 3199 3399 3599 3799 3999 4199 4399 4599 4799 4999 5199 5399 5599 5799 5999 6199 6399 6599 6799 6999 7049) 

for i in {0..35}; do


job_dir="${wd}/slurm_out/jobs/${dt}.genestat/"; mkdir -p $job_dir
err_dir="${wd}/slurm_out/errs/${dt}.genestat/"; mkdir -p $err_dir
outjob_dir="${wd}/slurm_out/outs/${dt}.genestat/"; mkdir -p $outjob_dir


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

Rscript $RSCRIPT ${start_pos[${i}]}  ${end_pos[${i}]}

EOF
" > ${job_script}

sbatch ${job_script}

done



