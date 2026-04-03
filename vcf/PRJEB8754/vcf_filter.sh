#!/usr/bin/env bash

set -euox pipefail

# bash filter_pass-orient-pos-sb.sh
bash filter_dp-vaf.sh
bash filter_blacklist.sh


