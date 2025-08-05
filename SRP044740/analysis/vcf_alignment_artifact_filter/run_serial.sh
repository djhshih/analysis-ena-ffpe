#!/bin/env/bash

set -euo pipefail

ls ../../data/vcf_filtermutectcalls_obf/*/*.vcf > dlazy_samples.txt

djobs dlazy_samples.txt filter_alignment_artifacts.sh

dlazy job
