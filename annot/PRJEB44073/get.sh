#!bin/env/bash

wget -O sample-info_stage0.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB44073&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,instrument_model,read_count,study_title,experiment_alias,fastq_md5,fastq_ftp,sample_alias&format=tsv&download=true&limit=0"