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
    
    target_path = f"{target_dir}/{base}"
    
    if os.path.exists(target_path):
        print(f"'{base}' already exists at '{target_dir}'\n\tSKIPPING...")
        continue
    
    os.link(path, target_path)
    print(f"Created hard link for: {base}, \n\tat {target_dir}\n")

print("Done.")
