import os
import glob

outdir = "rg_corrected"
os.makedirs(outdir, exist_ok=True)

paths = glob.glob("rg/*rg.txt")

for path in paths:

	sample_name = os.path.basename(path).split(".")[0]

	filename = os.path.basename(path)
	outpath = f"{outdir}/{filename}"

	with open(path, 'r') as file:
		rg = file.read()
		rg_fixed = rg.replace("ID:", f"ID:{sample_name}").replace("PU:", f"PU:{sample_name}")
	
	with open(outpath, 'w') as outfile:
		outfile.write(rg_fixed)
  
	print(f"Corrected RG for : {path}")
 
print(f"\nOutputs saved to {outdir}/\nDone.")



