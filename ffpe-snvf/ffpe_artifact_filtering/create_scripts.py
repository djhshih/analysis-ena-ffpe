import os
import glob
from tkinter.ttk import Separator
import polars as pl

## PRJEB8754
vcf_paths = glob.glob("../../vcf/PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered/*/*.vcf")
filtered_outdir = "../../ffpe-snvf/PRJEB8754/vcf_pass-orient-pos-sb_ad_filtered"
bam_dir = "../../data/PRJEB8754/bam"

models = ["mobsnvf", "vafsnvf", "sobdetector"]

outdir = "script"
os.makedirs(outdir, exist_ok=True)

for vcf_path in vcf_paths:
    if "Frozen" in vcf_path:
        continue
    
    vcf_filename = os.path.basename(vcf_path)
    sample_name = vcf_filename.removesuffix(".vcf")
    
    for model in models:
        content = f"#!/bin/bash\nbash ../{model}.sh.template '../../../data/bam/{sample_name}.bam' '../{vcf_path}' '../{filtered_outdir}'\n"
        filename = f"{outdir}/{model}_{sample_name}.sh"
        
        with open(filename, "w") as f:
            f.write(content)
            
        print(f"Created {model} filtering script for {vcf_path}")


## SRP044740
vcf_paths = sorted(glob.glob("../../vcf/SRP044740/vcf_filtered_pass_orientation/*/*.vcf"))
filtered_outdir = "../../ffpe-snvf/SRP044740/filtered_pass_orientation"
bam_dir = "../../data/SRP044740/bam"

models = ["mobsnvf", "vafsnvf", "sobdetector"]

outdir = "script"
os.makedirs(outdir, exist_ok=True)

for vcf_path in vcf_paths:
    if "FROZ" in vcf_path:
        continue
    
    vcf_filename = os.path.basename(vcf_path)
    sample_name = vcf_filename.removesuffix(".vcf")
    
    for model in models:
        content = f"#!/bin/bash\nbash ../{model}.sh.template '../{bam_dir}/{sample_name}.bam' '../{vcf_path}' '../{filtered_outdir}'\n"
        filename = f"{outdir}/{model}_{sample_name}.sh"
        
        with open(filename, "w") as f:
            f.write(content)
            
        print(f"Created {model} filtering script for {vcf_path}")


## SRP065941
vcf_paths = sorted(glob.glob("../../vcf/SRP065941/vcf_filtered_pass_orientation/*/*.vcf"))
filtered_outdir = "../../ffpe-snvf/SRP065941/vcf_filtered_pass_orientation"
bam_dir = "../../data/SRP065941/bam"

models = ["mobsnvf", "vafsnvf", "sobdetector"]

outdir = "script"
os.makedirs(outdir, exist_ok=True)

for vcf_path in vcf_paths:
    if "frozen" in vcf_path:
        continue
    
    vcf_filename = os.path.basename(vcf_path)
    sample_name = vcf_filename.removesuffix(".vcf")
    
    for model in models:
        content = f"#!/bin/bash\nbash ../{model}.sh.template '../{bam_dir}/{sample_name}/{sample_name}.bam' '../{vcf_path}' '../{filtered_outdir}'\n"
        filename = f"{outdir}/{model}_{sample_name}.sh"
        
        with open(filename, "w") as f:
            f.write(content)
            
        print(f"Created {model} filtering script for {vcf_path}")
        
        
## PRJEB44073
vcf_paths = sorted(glob.glob("../../vcf/PRJEB44073/vcf_filtered_pass_orientation/*/*.vcf"))
filtered_outdir = "../../ffpe-snvf/PRJEB44073/vcf_filtered_pass_orientation"
bam_dir = "../../data/PRJEB44073/bam"

annot = pl.read_csv("../../annot/PRJEB44073/sample-info_stage2.tsv", separator = "\t")
ffpe_samples = annot.filter(pl.col("preservation") == "FFPE").get_column("sample_alias").to_list()

vcf_paths = [path for path in vcf_paths if any(sample in path for sample in ffpe_samples)]

models = ["mobsnvf", "vafsnvf", "sobdetector"]

outdir = "script"
os.makedirs(outdir, exist_ok=True)

for vcf_path in vcf_paths:
    
    vcf_filename = os.path.basename(vcf_path)
    sample_name = vcf_filename.removesuffix(".vcf")
    
    for model in models:
        content = f"#!/bin/bash\nbash ../{model}.sh.template '../{bam_dir}/{sample_name}/{sample_name}.bam' '../{vcf_path}' '../{filtered_outdir}'\n"
        filename = f"{outdir}/{model}_{sample_name}.sh"
        
        with open(filename, "w") as f:
            f.write(content)
            
        print(f"Created {model} filtering script for {vcf_path}")
        
