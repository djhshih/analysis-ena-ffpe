#!/bin/bash
ls inputs/* > dlazy_samples.txt
djobs dlazy_samples.txt cromwell run ../../../wdl/bam_variant_mutect2_no_docker.wdl -i

pdlazy job -j 4
