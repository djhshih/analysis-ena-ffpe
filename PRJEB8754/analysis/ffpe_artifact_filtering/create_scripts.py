import os
import glob

vcf_paths = glob.glob("../../data/vcf_pass-n-orientation_ad_filtered/*/*.vcf")
filtered_outdir = "../../ffpe-snvf/vcf_pass-n-orientation_ad_filtered"

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

