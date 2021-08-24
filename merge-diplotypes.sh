#!/bin/bash
  
mergeTables (){
        dir=$1
        sub_path=$2
        out_txt=$3
        skip_write_header=$4

        counter=0
        for i in $dir/*; do
                counter=$((counter+1))

                in_txt=$i/$sub_path
                sample_name=$(basename $i)
                echo -ne "Processing [$counter]: $sample_name\r"

                ## Write header using first file
                if [[ $counter -eq 1 && $skip_write_header -ne 1 ]]; then
                        zcat $in_txt | head -n 1 |
                        awk '{print "sample""\t"$0}' |
                        gzip -c > $out_txt
                fi

                zcat $in_txt | tail -n +2 |
                awk -v sample_name="$sample_name" '{print sample_name"\t"$0}' |
                gzip -c >> $out_txt

        done

        echo -e '\n'
}

out_txt=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/all_diplotypes.txt.gz
mergeTables \
        /hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/diplotype-output \
        gene_diplotypes.txt.gz \
        $out_txt

