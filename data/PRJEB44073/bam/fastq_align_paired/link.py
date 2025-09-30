#!/usr/bin/env python3
import os
import glob
# import shutil

bam_paths = glob.glob("cromwell-executions/fastq_align_paired/*/call-bam_sort_coord/execution/*.bam")
bai_paths = glob.glob("cromwell-executions/fastq_align_paired/*/call-bam_sort_coord/execution/*.bai")
bam_bai_paths = bam_paths + bai_paths

target_dir = "../"
os.makedirs(target_dir, exist_ok=True)

for path in bam_bai_paths:
    base = os.path.basename(path)
    
    sample_name = base.split(".")[0]
    sample_dir = f"{target_dir}/{sample_name}"
    os.makedirs(sample_dir, exist_ok=True)
    target_path = f"{sample_dir}/{base}"
    
    if os.path.exists(target_path):
        print(f"'{base}' already exists at '{sample_dir}'\n\tSKIPPING...")
        continue
    
    os.link(path, target_path)
    print(f"Created hard link for: {base}, \n\tat {sample_dir}\n")

print("Done.")
