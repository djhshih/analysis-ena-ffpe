from pathlib import Path
from statistics import median
import pysam
import polars as pl
import subprocess

repo_root = Path("../..").resolve()


## Functions
def calculate_median_insert_size(bam_path, max_reads=2_000_000, threads=1):
	"""Compute median insert size from proper pairs.
	
	Samples up to max_reads proper pairs from the start of the BAM.
	Raises ValueError if no valid pairs found.
	"""
	sizes = []
	n_zero_tlen = 0
	
	with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
		for read in bam:
			if len(sizes) >= max_reads:
				break
			if (
				read.is_proper_pair          # both mates aligned sensibly together
				and not read.is_secondary    # not an alternate mapping
				and not read.is_supplementary # not a split-read fragment
				and not read.is_duplicate    # not a PCR/optical duplicate
				and not read.is_unmapped     # actually aligned somewhere
			):
				if read.template_length == 0:
					n_zero_tlen += 1
					continue
				sizes.append(abs(read.template_length))
	
	if n_zero_tlen:
		print(
			f"Skipped {n_zero_tlen} proper-pair reads with TLEN=0 in {bam_path}"
		)
	
	if not sizes:
		raise ValueError(f"No valid proper pairs found in {bam_path}")
	
	return int(median(sizes))

def calculate_median_coverage(
	bam_path,
	bed_path,
	threads=2,
	min_mapq=20,
	tmp_outdir="mosdepth",
	skip_if_exists=True,
):
	"""
	Run mosdepth on a BAM restricted to regions in `bed_path` and return the
	median per-region coverage from the total distribution.

	Parameters
	----------
	bam_path : Path
		Input BAM file.
	bed_path : Path
		BED file of regions to evaluate.
	threads : int
		Threads for mosdepth (`--threads`).
	min_mapq : int
		Minimum mapping quality (`--mapq`).
	tmp_outdir : str | Path
		Directory to write mosdepth outputs into.
	skip_if_exists : bool
		If True (default), skip running mosdepth when the expected
		`*.mosdepth.region.dist.txt` already exists. Set False to force re-run.

	Returns
	-------
	int
		Median region coverage (depth at which cumulative proportion >= 0.5).
	"""
	bam_path = Path(bam_path)
	bed_path = Path(bed_path)
	tmp_outdir = Path(tmp_outdir)
	tmp_outdir.mkdir(exist_ok=True, parents=True)

	sample_name = bam_path.name.split(".")[0]
	prefix = tmp_outdir / sample_name
	dist_file = tmp_outdir / f"{sample_name}.mosdepth.region.dist.txt"

	# ---- Run mosdepth (or skip if output already exists) ----
	if skip_if_exists and dist_file.exists():
		print(f"\t\t[mosdepth] Found existing {dist_file}, skipping run.")
	else:
		cmd = [
			"mosdepth",
			"--by", str(bed_path),
			"--no-per-base",
			"--threads", str(threads),
			"--mapq", str(min_mapq),
			str(prefix),
			str(bam_path),
		]
		print(f"[mosdepth] Running: {' '.join(cmd)}")
		try:
			subprocess.run(cmd, check=True)
		except FileNotFoundError as e:
			raise RuntimeError(
				"\t\tmosdepth executable not found in PATH. Install it or activate the right env."
			) from e
		except subprocess.CalledProcessError as e:
			raise RuntimeError(
				f"\t\tmosdepth failed with exit code {e.returncode} for {bam_path}"
			) from e

	if not dist_file.exists():
		raise FileNotFoundError(
			f"\t\tExpected mosdepth output not found: {dist_file}"
		)

	# ---- Parse distribution and compute median ----
	mosdepth_dist = pl.read_csv(
		dist_file,
		separator="\t",
		has_header=False,
		new_columns=["chrom", "depth", "proportion"],
	)
	median = (
		mosdepth_dist
		.filter((pl.col("chrom") == "total") & (pl.col("proportion") >= 0.5))
		.sort("depth", descending=True)
		.select("depth")
		.head(1)
		.item()
	)
	return median if median > 1 else 2

def sample_filter(vcf_paths: list, dataset: str) -> list:
	if dataset == "PRJEB44073":
		samples = [f"{i:04d}" for i in range(17, 33)]  # ['0017', '0018', ..., '0032']
		filtered = [
			path for path in vcf_paths
			if any(sample in Path(path).name for sample in samples)
		]
		return filtered
	else:
		filtered = [
			path for path in vcf_paths
			if "FFPE" in Path(path).name.upper()
		]
		return filtered

def get_bam_path(dataset: str, sample_name: str, fail_if_not_exit: bool = False, repo_root: Path = repo_root) -> Path:

	bam_root = Path(repo_root) / "data" / dataset / "bam"
	bam_candidates = [
		bam_root / sample_name / f"{sample_name}.bam",
		bam_root / f"{sample_name}.bam",
	]

	bam_path = next((p for p in bam_candidates if p.exists()), None)

	if bam_path is None:
		if fail_if_not_exit:
			raise FileNotFoundError(
				f"BAM for sample {sample_name!r} not found at any of: "
				+ ", ".join(str(p) for p in bam_candidates)
			)
		else:
			return None
	else:
		return bam_path

def get_bed_path(dataset: str, sample_name: str, repo_root: Path = repo_root) -> Path:
	exome_bed   = repo_root / "data/regions/sureselect-all-exon-v5_hg38_regions_clean.bed"
	exome_bed_padded = repo_root / "data/regions/sureselect-all-exon-v5_hg38_regions_200bp-pad.bed"
	srp044740_bed = repo_root / "pandepth/SRP044740/SRP044740_pileup-avg-dp15.bed"
	srp065941_bed = repo_root / "pandepth/SRP065941/SRP065941_pileup-avg-dp15.bed"

	if dataset == "PRJEB44073":
		return exome_bed, exome_bed_padded

	if dataset == "SRP044740":
		return srp044740_bed, srp044740_bed

	if dataset == "SRP065941":
		if sample_name.startswith(("T1", "T4")):
			return srp065941_bed, srp065941_bed
		if sample_name.startswith(("T2", "T3")):
			return exome_bed, exome_bed_padded
		raise ValueError(
			f"Sample {sample_name!r} in dataset {dataset!r} does not match any known prefix (T1-T4)."
		)

	raise ValueError(f"Unknown dataset: {dataset!r}")

## Create Execution Scripts
### Wrapper function
def generate_ffperase_scripts(
	dataset: str,
	variant_set: str,
	vcf_root: Path = repo_root / "vcf",
	vcf_ext: str = "vcf.gz",
	ref_path: Path = repo_root / "data/ref/Homo_sapiens_assembly38.fasta",
	script_outdir: Path = repo_root / "ffpe-snvf/ffpe_artifact_filtering/script_ffperase",
	ffperase_launcher: Path = repo_root / "common-ffpe-snvf/templates/ffpe-snvf/ffperase.sh",
	repo_root: Path = repo_root,
) -> list[Path]:
	"""
	Generate ffperase execution scripts for every (filtered) VCF in a dataset.

	Returns the list of script paths written.
	"""

	print(f"Processing Dataset: {dataset} | Variant Set: {variant_set}")
	script_outdir.mkdir(exist_ok=True, parents=True)

	# ---- Discover VCFs ----
	vcf_dir = vcf_root / dataset / variant_set
	vcf_paths = sorted(p.resolve() for p in vcf_dir.glob(f"*/*.{vcf_ext}"))
	vcf_paths = sample_filter(vcf_paths, dataset)

	if not vcf_paths:
		print(f"[warn] No VCFs found under {vcf_dir} after filtering.")

	written_scripts: list[Path] = []

	# ---- Generate per-sample scripts ----
	for i, vcf_path in enumerate(vcf_paths, start=1):
		sample_name = vcf_path.parent.name
		bam_path = get_bam_path(dataset, sample_name)
		mos_bed_path, ffperase_bed_path = get_bed_path(dataset, sample_name)

		print(f"\t{i}. Creating scripts for sample: {sample_name}")
		if not bam_path:
			print(f"\t\tBAM does not exist for sample: {sample_name}")
			continue

		ins_size = calculate_median_insert_size(bam_path)
		coverage = calculate_median_coverage(bam_path, mos_bed_path)

		result_outdir = repo_root / "ffpe-snvf" / dataset / variant_set / "ffperase" / sample_name
		result_outdir.mkdir(exist_ok=True, parents=True)

		script_file = script_outdir / f"{sample_name}_ffperase.sh"
		contents = f"""#!/usr/bin/env bash
set -euo pipefail

bash {ffperase_launcher} \\
\t--vcf           {vcf_path} \\
\t--bam           {bam_path} \\
\t--reference     {ref_path} \\
\t--bed           {ffperase_bed_path} \\
\t--coverage      {coverage} \\
\t--median-insert {ins_size} \\
\t--sample-name   {sample_name} \\
\t--step          full \\
\t--mutation-type snvs \\
\t--outdir        {result_outdir}
"""
		script_file.write_text(contents)
		written_scripts.append(script_file)

		print(f"\t\tWrote {script_file}  (median_coverage={coverage}, median_insert_size={ins_size})")


### Generate for each dataset
generate_ffperase_scripts(
	dataset="SRP065941",
	variant_set="filtered_pass-orientation",
)

generate_ffperase_scripts(
	dataset="PRJEB44073",
	variant_set="filtered_pass-orientation",
)

generate_ffperase_scripts(
	dataset="SRP044740",
	variant_set="filtered_pass-orientation",
)
