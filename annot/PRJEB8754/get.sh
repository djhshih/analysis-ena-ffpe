#!/bin/bash
# Download original sample information table

curl -o sample-info_unprocessed.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB8754&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,instrument_model,library_name,fastq_ftp,sample_alias&format=tsv&download=true&limit=0"
