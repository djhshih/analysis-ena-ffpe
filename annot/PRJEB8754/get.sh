#!/bin/bash
# Download original sample information table

curl -o sample-info_stage0.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB8754&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,sra_ftp,bam_ftp&format=tsv&download=true&limit=0"
