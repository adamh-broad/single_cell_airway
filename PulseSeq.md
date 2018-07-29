Analysis of 'pulse-seq' lineage trace of 65,000 airway epithelial cells
================

Load required R packages
------------------------

### Can be installed using 'install.package'

``` r
library(NMF)
```

    ## Loading required package: pkgmaker

    ## Loading required package: registry

    ## Loading required package: rngtools

    ## Loading required package: cluster

    ## NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 7/8

    ##   To enable shared memory capabilities, try: install.extras('
    ## NMF
    ## ')

``` r
library(RColorBrewer)
library(statmod)
library(ggplot2)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

``` r
library(Seurat)
```

    ## Loading required package: Matrix

``` r
library(Matrix)
### Load all the required functions for this analysis
source("Fxns.R")
```

    ## 2018-07-29 13:14:30 INFO: Loading default colors

### Load UMI count data from GEO and clustering and dimensionality reduction (output of 'Pre\_process' analysis)

``` r
## Downloading UMI count data
#download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_atlas_UMIcounts.txt.gz", destfile="GSE92332_atlas_UMIcounts.txt.gz")
## Reading UMI count data from file
## Note that processed UMI data already has ~500 contaminating immune (dendritic) cells, and all low quality (<1000 genes detected) cells removed
ps_umis = readRDS("PulseSeq_UMI_counts.rds")
info(sprintf("Data dimensions: %s" , paste(dim(ps_umis), collapse = "x")))
```

    ## 2018-07-29 13:14:36 INFO: Data dimensions: 27998x66265

``` r
ps = CreateSeuratObject(raw.data=ps_umis,min.cells = 0, min.genes = 0)
ps = NormalizeData(object = ps,normalization.method = "LogNormalize", scale.factor = 10000,display.progress = TRUE)
ps@dr = readRDS("dr.rds")
ps@meta.data= read.delim("meta_data_regressed.txt")
rownames(ps@meta.data) = ps@cell.names
rownames(ps@dr$tsne@cell.embeddings) = ps@cell.names
```

Plot the t-SNE, with cells split by time-point (Fig 3b)
-------------------------------------------------------

### Run PCA, t-SNE

``` r
d = FetchData(ps, c("timepoint", "color", "res.2_named", "tSNE_1", "tSNE_2"))
g = ggplot(d, aes(x=tSNE_1, y=tSNE_2)) + geom_point(data=subset(d, color=="Tom"), size=0.45, stroke=0, color="grey92") + geom_point(data=subset(d, color!="Tom"), size=0.5, stroke=0, color=default.cols(3)[3]) + facet_wrap(~timepoint)
print(g)
```

<img src="PulseSeq_figs/PulseSeq-tsne-1.png" style="display: block; margin: auto;" />
