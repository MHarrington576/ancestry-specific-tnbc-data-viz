# For CU Data Visualization 2026: ancestry-specific-tnbc-data-viz
***

> This program uses the [U.S. Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/) to access open-source data from The Cancer Genome Atlas (TCGA) project.
> This program uses the Xena TCGA Hub to access [open-source breast cancer clinical data](https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix) from TCGA.

## Usage
1. Install the GDC Client in a local Linux environment
2. To perform an analysis of mutation rates from MAF files:
    * Run make-gdc-ancestry-maf-manifest.sh
    * Create subdirectories for each ancestry group
    * In the corresponding subdirectory, download the MAF files from the GDC
    * Run maf-mutation-freq-analysis.R
3. To performa an analysis of clinical data:
    * Have Xena clinical matrix downloaded in the directory
    * Run extract-clinical-data.py (see comments on user-defined inputs)
    * Run subtype-stats-calculator.py
