#!/bin/env/bash

set -euo pipefail

ls ../../data/vcf_filtermutectcalls_obf/*/*.vcf > dlazy_samples.txt

djobs dlazy_samples.txt bash filter_alignment_artifacts.sh

pdlazy job -j 4
