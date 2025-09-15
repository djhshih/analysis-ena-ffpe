#!/bin/bash

Rscript MicroSEC_parallel.R ../PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered/microsec inputs/microsec.sample_info.tsv Y
Rscript MicroSEC_parallel.R ../SRP044740/vcf_filtered_pass_orientation/microsec inputs/microsec.sample_info.tsv Y

