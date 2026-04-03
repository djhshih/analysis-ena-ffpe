#!usr/in/env bash

set -euox pipefail

python non-exome-exclusion.py
python blacklist-exclusion.py
python germline-exclusion.py
python micr-exclusion.py
