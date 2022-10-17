# Installation
We use the downloaded version of the Mission Bio Tapestri platform (v2). This uses code from MissionBio, to which access has to be granted by MissionBio. Please contact the [maintainers](michael.scherer@crg.eu) of this pipeline for further instructions on the installation. Then, install the pipeline using the installation script. You will have to add some parts of the path, where the program has been installed to specific steps of the bash pipelines (variable $TAPESTRI).

```bash
sh TapestriPipeline-v2.0.1-Linux-x86_64.sh $TAPESTRI
```

The script will automatically generate a conda environment, which will be used in the bash pipeline scripts.

In addition, RnBeads is needed for parts of the Cellfinder script. It can be installed as follows.

```bash
conda create -p $RNBEADS bioconductor-rnbeads bioconductor-rnbeads.hg19 bioconductor-rnbeads.mm10
```

## Using the conda yml
Alternatively, you can use the yaml files in this repository to generate the corresponding conda environments.

```bash
conda env create -f tapestri.yml -p $TAPESTRI
conda env create -f rnbeads.yml -p $RNBEADS
```

## Running the pipeline
Depending whether you want to process bone marrow or Bcell data, you'll have to adapt the files run\_tapestri\_[BCells/BM].sh and tapestri\_pipeline\_[BCells/BM].sh, respectively. Please provide the paths to the downloaded files from GEO and an output folder.

The following parameters have to be changed:
-$OUTPUT: should be set to the desired output directory
-$INPUT: should be set to the path, where the input FASTQ files are stored
-$TAPESTRI: should be the name or path of the tapestri conda environment
-$RNBEADS: should be the name or path of the RnBeads conda environment
