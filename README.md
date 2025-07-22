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

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0127146

https://www.ebi.ac.uk/ena/browser/view/PRJEB8754

