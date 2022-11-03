# A collection of scripts for analyzing data obtained from single-cell Targeted Analysis of the Methylome (scTAM-seq)

## [manuscript](manuscript)
The collection of R scripts needed to re-generate the plots that we already generated for the manuscript. The scripts are ordered by the Figure/panel they are needed for.

If you want to re-generate all of the plots of the publication, you'll have to download several datasets and put them into the base folder of this repository (called scTAM-seq-scripts):
- The count matrices from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198019). They should be extracted into a folder called *data* in the base directory of this repository.
- The CITE-seq reference dataset from: https://figshare.com/ndownloader/files/28408917
- The reference bulk DNAm data from: https://figshare.com/account/projects/134468/articles/19323623. This data should also be placed into the *data* folder.

In addition, you'll have to install the following R-packages
- RnBeads: (```source('https://rnbeads.org/data/install.R')```)
- Seurat: (```BiocManager::install('Seurat')```)
- Signac: (```install.packages('Signac')```)
- monocle3: (```install.packages('monocle3')```)
- SeuratWrappers:(```devtools::install_github('satijalab/seurat-wrappers')```)
- viridis: (```install.packages('viridis')```)
- ComplexHeatmap: (```BiocManager::install('ComplexHeatmap')```)
- ggsci: (```install.packages('ggsci')```)
- rstan: (```install.packages('rstan')```)

## [tapestri](tapestri)
This contains the re-implementation of the Mission Bio tapestri pipeline to generate all of the files needed for downstream processing.

# Citation

If you use this software in your publication, please cite:

Bianchi, A., Scherer, M., Zaurin, R. et al. scTAM-seq enables targeted high-confidence analysis of DNA methylation in single cells. Genome Biol 23, 229 (2022). [https://doi.org/10.1186/s13059-022-02796-7](https://doi.org/10.1186/s13059-022-02796-7)
