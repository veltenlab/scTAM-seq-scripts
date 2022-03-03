# Installation
We use the downloaded version of the Mission Bio Tapestri platform (v2), which we downloaded from here:

```
wget https://dl.missionbio.io/onprem-miniconda/TapestriPipeline-v2.0.1-Linux-x86_64.sh
```

And install the pipeline using the installation script. You will have to add some parts of the path, where the program has been installed to specific steps of the bash pipelines (variable $TAPESTRI).

```
sh TapestriPipeline-v2.0.1-Linux-x86_64.sh $TAPESTRI
```

The script will automatically generate a conda environment, which will be used in the bash pipeline scripts.

In addition, RnBeads is needed for parts of the Cellfinder script. It can be installed as follows.

```
conda create -p $RNBEADS bioconductor-rnbeads bioconductor-rnbeads.hg19 bioconductor-rnbeads.mm10
```

## Using the conda yml
Alternatively, you can use the yaml files in this repository to generate the corresponding conda environments.

```
conda env create -f tapestri.yml -p $TAPESTRI
conda env create -f rnbeads.yml -p $RNBEADS
```
