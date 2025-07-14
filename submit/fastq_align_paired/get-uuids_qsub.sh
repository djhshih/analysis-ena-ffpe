#!/bin/bash

# Get workflow UUIDs for each sample

printf 'sample_id\tworkflow_uuid\n' > uuids.tsv

grep -Eo 'WorkflowActor-([^ ]+)' log/*.stdout | 
	sed -E 's|.*/([^.]+)\.stdout:WorkflowActor-(.*)+|\1\t\2|' \
	>> uuids.tsv
