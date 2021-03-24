# Altered transcriptome-proteome coupling indicates aberrant proteostasis in Parkinson’s disease 

This repository holds all code and data to reproduce results, figures and and supplementary data presented in "Altered transcriptome-proteome coupling indicates aberrant proteostasis in Parkinson’s disease".
The paper is currently available on [medrxiv](https://www.medrxiv.org/content/10.1101/2021.03.18.21253875v1) (not peer reviewed).


## Prerequisites

* [R >= 3.6.0](https://www.r-project.org/)
* All R packages used in the analysis are listed at the top of each Rmd file
* Many of these are bioconductor packages so it helps to have BiocManager installed:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(<"packagename">)
```
* The installation of ermineR requires JAVA and setting $JAVA_HOME variable: [ermineR](https://github.com/PavlidisLab/ermineR)  

## How to..  

### Pathway enrichment analysis

* The paper discusses results of pathway enrichment analysis that were performed on ranked gene lists.
* The ranking is based on the Pearson correlation coefficient (gene-wise, across-samples and within groups (HA, YG, PD). 
* The ranking takes into account the difference in correlation between groups (see paper methods).
* The genelist ranking (based on the correlation coefficients given in `./results/rds/gene_cl.rds` and `./results/rds/gene_cl_pd.rds`) are generated in `./pea.Rmd`.
* The pathway enrichment analysis is performed with [ermineR](https://github.com/PavlidisLab/ermineR) in `./pea.Rmd`. 
* Requirements for this analysis are the `.rds` files in `./results/rds/` and the GO annotations for ermineR in `./referenceData/go_daily-termdb.rdf-xml`. R packages needed are listed at the beginning of the `./pea.Rmd` file.
* The knitted, self_contained html of `./pea.Rmd` is also available (`./pea.html`, download and open in browser) and has all pathway enrichment results listed.

### Reproduce Paper Figures

* The main figures (excluding figure 1, which is a schematic figure) can be reproduced with `./analysis.Rmd`.
* The code loads `.rds` files need for the analysis form `./rdsData/`. These are generated in `./preprocess.Rmd` as described below.
* The figures are available in `./result_figs/Draft/` and `./result_figs/Supp/`.
* The `analysis.Rmd` also creates and saves the `gene_cl.rds` and `gene_cl_pd.rds` files to `./results/rds/` (needed in `./pea.Rmd`).
* Would not advice to knit this file. Hasnt been edited appropriately or tested. The purpose is mainly to document code.

### Rerun analysis from scratch

* The rawest data available are transcript-level counts in `./salmonOut` and protein intensities in `./rdsData/rawData.rds`.
* With these and the metadata in `./rdsData/info.rds` `./rdsData/info_rna.rds` and `./rdsData/info_rna_pd.rds` all results can be reproduced.
* The file `./preprocess.Rmd`, loads the required data and:
  * aggregates transcript-level counts from `./salmonOut` to gene-level using [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) and the annotation file `./referenceData/gencode.v32.annotation.gtf` (needs to be unzipped).
  * filters out zeros and applies batch correction to the proteomics data in `./rdsData/rawData.rds`
  * saves all generated `.rds` files to `./rdsData`.
* Utility functions (also the batch correction function) are loaded from `./functions.R`.
* Next steps would be to go through `./analysis.Rmd` and then `./pea.Rmd`
* Would not advice to knit `./analysis.Rmd`. Hasnt been edited appropriately or tested. The purpose is mainly to document code.

## Content

* `./preprocess.Rmd` Code to preprocess data (Transcript counts and raw protein intensities)
* `./analysis.Rmd` Code to reproduce paper figures 
* `./pea.Rmd`, `./pea.html` Code to run pathway enrichment analysis and knitted html to browse through results
* `./referenceData/` Holds reference data needed for the analysis. Remember to unzip compressed files (before running `./preprocess.Rmd`)
* `./salmonOut/` Output from [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) for each sample in the analysis. Details on parameters of the salmon command see methods section of the paper. Sample_ids (Salmon dir. names) correspond to RNAseq_id (see also supplementary File 1, or './rdsData/info.rds`.
* `./rdsData/` 
  * `info.rds` Holds all sample metadata. (The proteomics data uses reporter.intensity.id as sample identifier, this can be mapped to RNAseq_id for the RNASeq data).
  * `info_rna.rds` and `info_rna_pd.rds` Holds sample metadata for RNA data (including RIN) for the groups HA, YG and PD respectively. 
  * `cpmAll.rds` and `cpmAll_pd.rds` Holds gene-level counts in counts per million (CPM), after filtering. These files are created in `./preprocess.Rmd`
  * `rawData.rds` Holds protein intensities before filtering and before normalization. This file is used in `./preprocess.Rmd` to generate the batch corrected matrix which is then used in the analysis.
* `./result_figs` All paper figures (some might have been edited with graphics software afterwards and thus dont resemble exactly the same figure as in the paper)
* `./results` Holds the rds files for the pathway enrichment analysis (i.e. correlation coefficients for the groups (YG, HA and PD).
* Each `.Rmd` file loads functions from `./functions.R`. 

## Version info 

### R sessionInfo() output  

```
> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8          LC_NUMERIC=C                  LC_TIME=en_GB.UTF-8           LC_COLLATE=en_GB.UTF-8        LC_MONETARY=en_GB.UTF-8
 [6] LC_MESSAGES=en_GB.UTF-8       LC_PAPER=en_GB.UTF-8          LC_NAME=en_GB.UTF-8           LC_ADDRESS=en_GB.UTF-8        LC_TELEPHONE=en_GB.UTF-8
[11] LC_MEASUREMENT=en_GB.UTF-8    LC_IDENTIFICATION=en_GB.UTF-8

attached base packages:
 [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] mixOmics_6.10.9             lattice_0.20-40             MASS_7.3-51.5               boot_1.3-24                 cowplot_1.1.0               ggplotify_0.0.5
 [7] ggrepel_0.8.2               pheatmap_1.0.12             stringr_1.4.0               knitr_1.29                  readr_1.3.1                 tidyr_1.1.1
[13] gridExtra_2.3               readxl_1.3.1                ggfortify_0.4.10            preprocessCore_1.48.0       DESeq2_1.26.0               SummarizedExperiment_1.16.1
[19] DelayedArray_0.12.3         BiocParallel_1.20.1         matrixStats_0.56.0          ggpubr_0.4.0                tibble_3.0.3                ggbiplot_0.55
[25] scales_1.1.1                plyr_1.8.6                  ensembldb_2.10.2            AnnotationFilter_1.10.0     GenomicFeatures_1.38.2      AnnotationDbi_1.48.0
[31] Biobase_2.46.0              GenomicRanges_1.38.0        GenomeInfoDb_1.22.1         IRanges_2.20.2              S4Vectors_0.24.4            BiocGenerics_0.32.0
[37] RColorBrewer_1.1-2          dplyr_1.0.2                 colorspace_1.4-1            xlsx_0.6.4.2                kableExtra_1.1.0            colormap_0.1.4
[43] dendextend_1.14.0           igraph_1.2.5                ermineR_1.0.1.9000          patchwork_1.0.1             ComplexHeatmap_2.2.0        viridis_0.5.1
[49] viridisLite_0.3.0           ggplot2_3.3.2               reshape2_1.4.4              magrittr_1.5                nvimcom_0.9-102

loaded via a namespace (and not attached):
  [1] backports_1.1.8          circlize_0.4.10          Hmisc_4.4-1              BiocFileCache_1.10.2     lazyeval_0.2.2           splines_3.6.3
  [7] digest_0.6.25            htmltools_0.5.0          checkmate_2.0.0          memoise_1.1.0            cluster_2.1.0            openxlsx_4.1.5
 [13] annotate_1.64.0          Biostrings_2.54.0        rARPACK_0.11-0           askpass_1.1              prettyunits_1.1.1        jpeg_0.1-8.1
 [19] blob_1.2.1               rvest_0.3.6              rappdirs_0.3.1           haven_2.3.1              xfun_0.16                crayon_1.3.4
 [25] RCurl_1.98-1.2           jsonlite_1.7.0           genefilter_1.68.0        survival_3.1-8           glue_1.4.1               gtable_0.3.0
 [31] zlibbioc_1.32.0          XVector_0.26.0           webshot_0.5.2            GetoptLong_1.0.2         V8_3.2.0                 car_3.0-9
 [37] shape_1.4.4              abind_1.4-5              DBI_1.1.0                rstatix_0.6.0            Rcpp_1.0.5               xtable_1.8-4
 [43] htmlTable_2.0.1          progress_1.2.2           clue_0.3-57              gridGraphics_0.5-0       foreign_0.8-75           bit_4.0.4
 [49] Formula_1.2-3            htmlwidgets_1.5.1        httr_1.4.2               ellipsis_0.3.1           pkgconfig_2.0.3          gemmaAPI_2.0.3.9000
 [55] XML_3.99-0               rJava_0.9-13             nnet_7.3-13              dbplyr_1.4.4             locfit_1.5-9.4           tidyselect_1.1.0
 [61] rlang_0.4.7              munsell_0.5.0            cellranger_1.1.0         tools_3.6.3              generics_0.0.2           RSQLite_2.2.0
 [67] broom_0.7.0              evaluate_0.14            bit64_4.0.2              zip_2.1.0                purrr_0.3.4              xml2_1.3.2
 [73] biomaRt_2.42.1           compiler_3.6.3           rstudioapi_0.11          curl_4.3                 png_0.1-7                ggsignif_0.6.0
 [79] geneplotter_1.64.0       stringi_1.4.6            RSpectra_0.16-0          forcats_0.5.0            ProtGenerics_1.18.0      Matrix_1.2-18
 [85] vctrs_0.3.2              pillar_1.4.6             lifecycle_0.2.0          BiocManager_1.30.10      GlobalOptions_0.1.2      corpcor_1.6.9
 [91] data.table_1.13.0        bitops_1.0-6             rtracklayer_1.46.0       R6_2.4.1                 latticeExtra_0.6-29      rio_0.5.16
 [97] assertthat_0.2.1         xlsxjars_0.6.1           openssl_1.4.2            rjson_0.2.20             withr_2.2.0              GenomicAlignments_1.22.1
[103] Rsamtools_2.2.3          GenomeInfoDbData_1.2.2   hms_0.5.3                rpart_4.1-15             rvcheck_0.1.8            rmarkdown_2.3
[109] carData_3.0-4            base64enc_0.1-3          ellipse_0.4.2
```
 
