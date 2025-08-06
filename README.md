# Analysis of paired FFPE - fresh frozen samples

## Description
Here we use data from ENA PRJEB8754 to compare the performance of various FFPE artifact filtering tools on paired FFPE and Fresh Frozen samples in the dataset.

## Dependencies
- python
- R
- SOBDetector (included)
- htspan ([djhshih/htspan](https://github.com/djhshih/htspan))
- samtools
- bcftools
- cromwell
- GATK
- gsutil
- rgsam ([djhshih/rgsam](https://github.com/djhshih/rgsam))
- dlazy ([djhshih/dlazy](https://github.com/djhshih/dlazy/tree/main))
- openjdk

### R packages
- tidyverse
- io
- precrec
- jsonlite
- argparser
- glue
- patchwork
- grid
- hrbrthemes
- viridis
- MicroSEC ([MANO-B/MicroSEC](https://github.com/MANO-B/MicroSEC))
- Rsamtools
- BiocGenerics
- Biostrings
- doParallel

### Python Libraries
- polars
- pysam
- pandas
- numpy
- matplotlib
- seaborn
- elementpath

## Usage

1. Clone the repository:
    ```bash
    git clone https://github.com/djhshih/analysis-ena-ffpe.git
    ```

2. Install all the dependencies listed above. 

3. Move to the `annot` directory and run `annot.py`:
    ```bash
    python annot.py
    ```
    This will generate several annotation files which will be used later in our analysis.

4. Move to the `data` directory and run:

    ```bash
    bash get_gatk_data.sh
    cd ref
    bash get.sh
    ```
    This will download the reference genome and additional data from broadinstitute.org to be used later in the analysis.

5. Move to the `data/fq` directory and run `fastq_ftp_download.sh`:
    ```bash
    bash fastq_ftp_download.sh
    ```
    This will download the FASTQ files from ENA.

6. Move to the `data/rg` directory and run `get-fq-head.sh` followed by `rgsam.sh`:
    ```bash
    bash get-fq-head.sh
    bash rgsam.sh
    ```
    This will generate the read group information for the FASTQ files.

7. Move to the `analysis/fastq` aligned pair directory and run the following:
    ```bash
    python prepare.py
    bash align_fastqs.sh
    python link.py
    ```
    This will align the FASTQ files to the reference genome and generate BAM files.

8. Move to the `analysis/bam_variant_mutect2` directory and run:
    ```bash
    python prepare.py
    bash call_variants_mutect2.sh
    python link.py
    ```
    This will call variants using Mutect2 on the aligned BAM files.

9. Move to the `analysis/ffpe_artifact_filtering` directory and run:
    ```bash
    python create-scripts.py
    bash ffpe-snvf.sh
    python make-microsec-inputs.py
    bash ffpe-snvf_microsec.sh
    ```
    This will prepare the outputs from each ffpe artifact filtering tool to be used for performance evaluation.
    Outputs will be saved in the `ffpe-snvf`
10. Move to the `analysis/evauations` directory and run:
    ```bash
    python eval-mutation-signatures.py
    python make-summary-statistics.py
    Rscript make-eval-plots.R
    ```
    This will generated performance evaluaation plots, tables and summary statistics in the `evaluations` directory.

## Data source

### PRJEB8754
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0127146

https://www.ebi.ac.uk/ena/browser/view/PRJEB8754

### SRP044740
Study Title: NA

exomes of 13 FFPE breast tumor samples and 13 corresponding frozen samples

https://www.ebi.ac.uk/ena/browser/view/SRP044740


### PRJEB44073
Study Title: The Mutational Concordance of Fixed Formalin Paraffin Embedded and Fresh Frozen Gastro-oesophageal Tumours Using Whole Exome Sequencing
https://www.mdpi.com/2077-0383/10/2/215

https://www.ebi.ac.uk/ena/browser/view/PRJEB44073

WXS of 16 matched fresh frozen and FFPE gastro-oesophageal tumours


## Variant Calling Decisions

https://www.biostars.org/p/448808/

https://gatk.broadinstitute.org/hc/en-us/community/posts/360059696811-Mutect2-not-calling-a-4-bp-deletion-in-BRCA1-with-50-AF

https://gatkforums.broadinstitute.org/gatk/discussion/24507/mutect2-repeatedly-not-detecting-somatic-variant-idh2-r172k-with-solid-read-support-and-5-af

https://gatk.broadinstitute.org/hc/en-us/community/posts/360057582511-HaplotypeCaller-data-generated-from-amplicon-sequencing

https://gatk.broadinstitute.org/hc/en-us/community/posts/12450796994459-Asking-for-advice-on-Mutect2-calling-in-somatic-but-amplicon-data

https://www.reddit.com/r/bioinformatics/comments/a71z2f/running_mutect2_with_dontusesoftclippedbases/


TO Do: Add more information about what changes are made to the variant calling workflow. What is skipped. What is added. etc.

## Issues


### Mutect2 and Amplicon Sequencing
There are too many reads in amplicon sequencing so mutect 2 needed to be run with special flags suitable for amplicon sequencing.

Extra flags used for Mutect2: 

```
--disable-read-filter NotDuplicateReadFilter --downsampling-stride 50 --linked-de-bruijn-graph --max-reads-per-alignment-start 0 --annotations-to-exclude StrandBiasBySample --annotations-to-exclude ReadPosRankSumTest

```

Reasoning for these flags:
- `--disable-read-filter NotDuplicateReadFilter`: This is used to avoid filtering out reads that are marked as duplicates, which is common in amplicon sequencing.
- `--downsampling-stride 50`: This is used to reduce the number of reads processed by Mutect2, which is necessary due to the high coverage in amplicon sequencing. 
- `--linked-de-bruijn-graph`: This is used to improve the sensitivity of variant calling in amplicon sequencing.
- `--max-reads-per-alignment-start 0`: This is used to avoid downsampling reads at the start of an alignment, which can be problematic in amplicon sequencing.
- `--annotations-to-exclude StrandBiasBySample`: This is used to exclude the strand bias annotation, which is not relevant for amplicon sequencing as it is often seen that different PCR primer pairs have differet efficiencies leading to once strand being amplified faster than the other leading to strand bias.
- `--annotations-to-exclude ReadPosRankSumTest`: This is used exclude penalizing variants which appear reads near the start or the end of the read in too many reads. This is not relevant for amplicon sequencing as the reads are of fixed length and the variants are expected to appear at the same position in most of the reads.


### FFPE vs Fresh Frozen Variant Intersection
From the VCFs we see that a lot of variants are called in FFPE but not in Fresh Frozen and vice versa. Some of this can be attributed to tumor heterogeniety and artifacts. But opening up the BAM files of an FFPE and matched Fresh Frozen sample in IGV, we can see that a variant called in FFPE is still present in the Fresh Frozen sample or vice versa albeit with either lower number of reads supporting it or with stronger strand bias.

Example: In the Pat01 sample from PRJEB8754, we can see that for Chr1:114,716,199 a G>T mutation was called in the FF VCF but was not Present in the FFPE VCF. However, upon inspection in IGV we can see that the variant is present in the FFPE sample too despite in exremely low proportion.

![alt text](image.png)

In another example from the same patient, we can see that in CHR2:29220907 a A>G mutation was called in both FFPE and but was not present in the FF VCF. Upon isnpection in IGV we can see that the variant is present in the FF sample too and also with a higher propertion of reads supporting the alternative allele, albeit with a stronger strand bias.

![alt text](image-1.png)


Some variants appear in both FF and FFPE but get called in one due to a slightly higher frequency. We have to find a way to adjust the variant calling so that the variants are more reliable for benchmarking.

