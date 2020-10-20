## Pan-CoV PhIP-Seq Manuscript Repository

This repository contains all necessary code for analysis described in <>
The repository contains the code for a [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) pipeline,
dubbed [phip-flow](https://github.com/matsengrp/phip-flow), that produces an xarray dataset 
containing sample and peptide metadata tied to the raw counts matrix like so:

<p align="center">
  <img src="cartoons/Xarray.png" width="350">
</p>

With the colored columns representing coordinate dimensions and the grey squares representing data arrays
organized by respective shared dimensions.

Once this dataset is obtained, the jupyter notebooks containing all analysis may be run using the environment
provided here along with the [phippery](git@github.com:matsengrp/phippery.git) package installed

### Running the alignment pipeline

All raw fastq files for the samples can be obtained upon request and will be provided in a tarball, NGS.tar.
place this file in the alignment pipeline directory and run `tar -xvf NGS.tar` to produce the `NGS` Directory.

The pipeline can be run anywhere with `Nextflow` and `Docker`/`Singularity` installed. 
The scripts inside the directory provide an example of how to run the pipeline on our local cluster, gizmo.

To 

### Running the analysis notebooks


