# A revised airway epithelial hierarchy includes CFTR-expressing ionocytes

## Single-cell analysis code
This repo contains <a href="https://github.com/adamh-broad/single_cell_airway/blob/master/PulseSeq.md">quick start</a> R code for analysis of 66,265 airway epithelial cells. The easiest way to run it is to clone the repo (download ZIP button at top) and open the file 'Analysis.Rmd' in <a href="https://www.rstudio.com/">RStudio</a>. Each code chunk can then be run by clicking the run button (top right corner of each chunk).

There are several R (3.3 or later is required) packages needed to run the code, each can be installed using the 'install.packages' command. For example, to install the NMF package used to render a heatmap, run the R command 'install.packages("NMF")'. The 'sva' package used for batch correction must be installed from Bioconductor:

```{r }
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
``` 

## Abstract
We combine single-cell RNA-seq and in vivo lineage tracing to study the cellular composition and hierarchy of the murine tracheal epithelium. We identify a new rare cell type, the pulmonary ionocyte; functional variations in club cells based on their proximodistal location; a distinct cell type that resides in high turnover squamous epithelial structures that we named “hillocks”; and disease-relevant subsets of tuft and goblet cells. With a new method, Pulse-Seq, we show that tuft, neuroendocrine, and ionocyte cells are continually and directly replenished by basal progenitor cells. Remarkably, the cystic fibrosis gene, CFTR, is predominantly expressed in both mouse and human pulmonary ionocytes. Genetic perturbation of murine ionocytes causes a loss of Cftr expression and disrupts airway fluid and mucus physiology, which are also altered in cystic fibrosis. By associating cell type-specific expression programs with key disease genes, we establish a new cellular narrative for airways disease. 

Experimental workflow            |  GFP(Foxi1)+ ionocytes
:-------------------------:|:-------------------------:
![](https://github.com/adamh-broad/single_cell_airway/blob/master/fig1a.jpg)  |  ![](https://github.com/adamh-broad/single_cell_airway/blob/master/fox1_gfp.jpg)

## Related Resources
* <a href="https://www.nature.com/"> Associated paper </a>
* <a href="https://portals.broadinstitute.org/single_cell/study/airway-epithelium">Single Cell Portal (Broad Institute)</a>
* <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354">GEO Database link</a>

For questions or issues email:
ahaber -at- broadinstitute.org
