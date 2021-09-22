#!/bin/bash



wd="/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/epigenomic-annotation"
vcf_dir_hmf="/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics"
tadBedFileInput="$wd/annotation-beds/chromatin-status-tads.bed.bgz"  #args[3]
binBedFileInput="$wd/annotation-beds/chromatin-status-bins.bed.bgz"  #args[4]
repTimingFileInput="$wd/annotation-beds/repliseq_encode_mean_binned.bed.bgz"  #args[5]
repTimingFileInputBoxtel="$wd/annotation-beds/boxtel_paper_all_RepliSeq_median.bed.bgz" #args[6]
cpgInputFile="$wd/annotation-beds/cpgIslandExtUnmasked.bed.bgz"  #args[7]
repOrientationInputFile="$wd/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz"  #args[8]

out_dir="$wd/HMF/annotated-vcfs"   #args[9]
job_dir="$wd/logs"

rowNo=$(cat ${wd}/metadata_whitelisted.tsv | sed -e '1d' | wc -l)
primaryTumorLocs=$(cat ${wd}/metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 7)
sampleNames=$(cat ${wd}/metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 3)
sampleIds=$(cat ${wd}/metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 2)


# marking the time for jobs
dt=$(date +%y%m%d_%T)



#samples.txt contains sample IDs obtained from spargling-genomics database. It currently contains only three sample IDs!

for i in $(cat metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 7 | sort | uniq); do mkdir -p bin/${i}; done #This is to make the cancer-specific folder structure!

counter=0
for i in $(eval echo "{1..$rowNo}"); do

primaryTumorLoc=$(echo ${primaryTumorLocs} | cut -d " " -f ${i})
sampleName=$(echo ${sampleNames} | cut -d " " -f ${i})
sampleId=$(echo ${sampleIds} | cut -d " " -f ${i})


path_to_vcf="${vcf_dir_hmf}/${sampleName}/purple/${sampleId}.purple.somatic.vcf.gz"

#if [ "${primaryTumorLoc}" == "Lung" ]; then


if [ $counter -lt 3 ]; then

if [ -f ${path_to_vcf} ]; then



#if [[ -f $out_dir/${primaryTumorLoc}/${sampleId}/annotated-${sampleId}.txt ]]; then
        
#echo "${sampleId} is already annotated!"

#else

if [[ -d $out_dir/cancer-type/${primaryTumorLoc}/${sampleId} ]]; then
rm -rf $out_dir/cancer-type/${primaryTumorLoc}/${sampleId};
echo "Directory removed!";
fi;

mkdir ${out_dir}/cancer-type/${primaryTumorLoc}/${sampleId}

write_dir="${out_dir}/cancer-type/${primaryTumorLoc}"


# marking the time

if [[ ! -d ${job_dir}/jobs/${primaryTumorLoc}/${dt} ]]; then
mkdir -p ${job_dir}/jobs/${primaryTumorLoc}/${dt};
mkdir -p ${job_dir}/errs/${primaryTumorLoc}/${dt};
mkdir -p ${job_dir}/outs/${primaryTumorLoc}/${dt};

echo ${dt} >> ${wd}/sample-ids/cancer-type/${primaryTumorLoc}/not-annotated.txt
echo ${dt} >> ${wd}/sample-ids/cancer-type/${primaryTumorLoc}/annotated.txt

if [ $(cat annotated.txt | grep "${dt}" | wc -l) -lt 1 ]; then
echo ${dt} >> ${wd}/sample-ids/cancer-type/not-annotated.txt
echo ${dt} >> ${wd}/sample-ids/cancer-type/annotated.txt
fi

fi


job_script=$job_dir/jobs/${primaryTumorLoc}/${dt}/${sampleId}.sh


if [[ ! -p $job_script ]]; then
        touch $job_script;
fi


echo "#!/bin/bash

#SBATCH --job-name=${primaryTumorLoc}_$(basename $job_script)
#SBATCH --output=$job_dir/outs/${primaryTumorLoc}/${dt}/${sampleId}.txt
#SBATCH --error=$job_dir/errs/${primaryTumorLoc}/${dt}/${sampleId}.txt
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=a.movasati@uu.nl

guixr load-profile $wd/.guix-profile-berner-proj --<<EOF

Rscript $wd/variantAnnotator.R $sampleId $path_to_vcf $tadBedFileInput $binBedFileInput $repTimingFileInput $repTimingFileInputBoxtel $cpgInputFile $repOrientationInputFile $write_dir

if [[ ! -f ${out_dir}/cancer-type/${primaryTumorLoc}/${sampleId}/annotated-${sampleId}.txt ]]; then
echo ${sampleId} >> ${wd}/sample-ids/cancer-type/${primaryTumorLoc}/not-annotated.txt
echo ${sampleId} >> ${wd}/sample-ids/cancer-type/not-annotated.txt
else
echo ${sampleId} >> ${wd}/sample-ids/cancer-type/${primaryTumorLoc}/annotated.txt
echo ${sampleId} >> ${wd}/sample-ids/cancer-type/annotated.txt
fi

EOF" > $job_script

if [[ ! -f ${job_script}.done ]]; then
	sbatch $job_script
else
	echo "Path for storing SBATCH outputs not available. Please construct the correct directory structure before running this program!"
fi


#fi #for samples that already are annotated if statement


else

echo '${sampleId} does not exist!'

fi #for path availability if statement

counter=$((${counter}+1))


else

break

fi #for counter if statement


#fi #for cancer type if statement

done








