#!/usr/bin/env bash
set -euox pipefail

bash make-pileup.sh
bash make-bed.sh
