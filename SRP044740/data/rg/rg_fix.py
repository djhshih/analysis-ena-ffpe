import os
import glob

outdir = "rg_corrected"
os.makedirs(outdir, exist_ok=True)

paths = glob.glob("*rg.txt")

for path in paths:

	filename = os.path.basename(path)
	outpath = f"{outdir}/{filename}"

	with open(path, 'r') as file:
		rg = file.read()
		rg_fixed = rg.replace("ID:", f"ID:{paths[0].split("_")[0]}").replace("PU:", f"PU:{path.split("_")[0]}")
	
	with open(outpath, 'w') as outfile:
		outfile.write(rg_fixed)
  
	print(f"Corrected RG for : {path}")
 
print(f"\nOutputs saved to {outdir}/\nDone.")



