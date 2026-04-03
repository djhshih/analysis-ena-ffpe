#!/usr/bin/env bash

set -euox pipefail

bash filter_pass-orientation.sh
bash filter_exome_target.sh
bash filter_blacklist.sh
bash filter_germline.sh

