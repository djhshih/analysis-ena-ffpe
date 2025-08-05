import polars as pl
import os
import glob

sd_paths = glob.glob("../../ffpe-snvf/vcf_pass-orient-pos-sb_ad_filtered/sobdetector/*/*.sobdetector.snv")

for path in sd_paths:
    df = pl.read_csv(path, separator="\t", infer_schema_length=10000)
    artifacts = df.filter(pl.col("artiStatus") == "artifact")
    snv = df.filter(pl.col("artiStatus") == "snv")
    
    dirname = os.path.dirname(path)
    basename = os.path.basename(path).strip(".snv")
    
    artifacts.write_csv(f"{dirname}/{basename}.artifacts.snv", separator="\t")
    snv.write_csv(f"{dirname}/{basename}.filtered.snv", separator="\t")
    


