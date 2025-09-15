#!/bin/env/bash

wget "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP044740&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,instrument_model,study_title,fastq_md5,fastq_ftp,sra_ftp,sample_alias,sample_title&format=tsv&download=true&limit=0" -O sample-info_stage0.tsv