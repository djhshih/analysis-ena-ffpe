#!/usr/bin/env bash

wget -O annotation_stage0.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP065941&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,instrument_model,read_count,study_title,experiment_alias,run_alias,fastq_md5,fastq_ftp,submitted_ftp,sra_ftp,sample_alias,sample_title,bam_ftp&format=tsv&download=true&limit=0"