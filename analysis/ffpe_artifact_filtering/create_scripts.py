import os
import glob

vcf_paths = glob.glob("../../data/vcf/*/*.vcf")

models = ["mobsnvf", "vafsnvf", "sobdetector"]

outdir = "script"
os.makedirs(outdir, exist_ok=True)

for path in vcf_paths:
	vcf_filename = os.path.basename(path)
	sample_name = vcf_filename.removesuffix(".vcf")
	
	for model in models:
     
		content = f"#!/bin/bash\nbash ../{model}.sh.template '{sample_name}.bam' '{vcf_filename}'\n"
		filename = f"{outdir}/{model}_{sample_name}.sh"
		
		with open(filename, "w") as f:
			f.write(content)

