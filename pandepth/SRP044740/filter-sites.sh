#!/usr/bin/env bash

set -euo pipefail

N_SAMPLES=40

zcat SRP044740.SiteDepth.gz | \
	awk -v n="$N_SAMPLES" '$3 >= 15*n {print $0}' \
	> SRP044740.SiteDepth.dp600.tsv
