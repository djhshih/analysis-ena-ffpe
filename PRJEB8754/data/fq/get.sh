#!/bin/bash

urls=$(sed 1d ../../annot/sample-info_stage0.tsv | cut -f 7)

generate_script() {
	for url in ${urls}; do
		fastq1=$(echo ${url} | cut -d ';' -f 1)
		echo wget -c $fastq1
		fastq2=$(echo ${url} | cut -d ';' -f 2)
		echo wget -c $fastq2
	done
}

generate_script > download.sh
