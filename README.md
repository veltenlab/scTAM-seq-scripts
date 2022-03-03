# A collection of scripts for analyzing data obtained from single-cell Targeted Analysis of the Methylome (scTAM-seq)

## [manuscript](manuscript)
The collection of R scripts needed to re-generate the plots that we already generated for the manuscript. The scripts are ordered by the Figure/panel they are needed for.

If you want to re-generate all of the plots of the publication, you'll have to download several datasets and put them into the base folder of this repository (called scTAM-seq-scripts):
- The count matrices from GEO
- The CITE-seq reference dataset from: https://figshare.com/ndownloader/files/28408917
- The reference bulk DNAm data from: 

In addition, you'll have to install the following R-packages
- RnBeads: (```source('https://rnbeads.org/data/install.R')```)
- Seurat: (```BiocManager::install('Seurat')```)
- Signac: (```install.packages('Signac')```)
- monocle3: (```install.packages('monocle3')```)
- SeuratWrappers:(```devtools::install_github('satijalab/seurat-wrappers')```)
- viridis: (```install.packages('viridis')```)
- pheatmap: (```install.packages('pheatmap')```)
- ggsci: (```install.packages('ggsci')```)

## [tapestri](tapestri)
This contains the re-implementation of the Mission Bio tapestri pipeline to generate all of the files needed for downstream processing.

