#!/usr/bin/env python3
import os
import sys
import glob
import pysam

module_path = os.path.abspath("../../python")
sys.path.append(module_path)

import mutation_signatures as ms


def eval_mut_sig(vcf_path, reference_genome, outdir_root) -> None:
	sample_name = os.path.basename(vcf_path).removesuffix(".gz").removesuffix(".vcf")

	outdir = f"{outdir_root}/mutation_signatures/{sample_name}"
	os.makedirs(outdir, exist_ok=True)

	variants_96c = ms.variants_mut_profile(vcf_path, sample_name, reference_genome)
	mutaion_counts_96c = ms.create_96c_mutation_counts(variants_96c)

	print(f"{sample_name}:")

	variants_96c.write_csv(f"{outdir}/{sample_name}_96c_mutation_profiles.tsv", separator="\t")
	print(f"\tSaved mutation profiles to: {outdir}/{sample_name}_96c_mutations.tsv")

	mutaion_counts_96c.write_csv(f"{outdir}/{sample_name}_96c_mutation_count.tsv", separator="\t")
	print(f"\tSaved 96 channnel mutation counts to: {outdir}/{sample_name}_96c_mutation_count.tsv")

	ms.SBS96_plot(mutaion_counts_96c[:, 2].to_numpy(), label=sample_name, height=3, width=10, s=8, xticks_label=False, file=f"{outdir}/{sample_name}_96c_mutation_signature.pdf")
	print(f"\tSaved 96 channnel mutation counts to: {outdir}/{sample_name}_96c_mutation_signature.pdf")


def main() -> None:
	vcf_paths = sorted(glob.glob("../../data/vcf/*/*.vcf"))
	reference_genome_path= "../../data/ref/Homo_sapiens_assembly38.fasta"
	reference_genome = pysam.FastaFile(reference_genome_path)
	outdir_root = "../../evaluations/alt_flags_no_pon"

	for vcf_path in vcf_paths:
		eval_mut_sig(vcf_path, reference_genome, outdir_root)
  

if __name__ == "__main__":
    main()


