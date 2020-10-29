## Pan-CoV PhIP-Seq Manuscript Repository

This repository contains all necessary code for analysis described in <>
The repository contains the code for a [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) pipeline,
dubbed [phip-flow](https://github.com/matsengrp/phip-flow), that produces an xarray dataset 

Once this dataset is obtained, the jupyter notebooks containing all analysis may be run using the environment
provided here along with the [phippery](git@github.com:matsengrp/phippery.git) package installed

### Running the alignment pipeline

All raw fastq files for the samples can be obtained upon request and will be provided in a tarball, `NGS.tar`.
Place this file in the `alignment-pipeline/` directory and run `tar -xvf NGS.tar` to produce the `NGS` Directory.

The pipeline can be run anywhere with `Nextflow` and `Docker` (or `Singularity`) installed. 
The scripts inside the directory provide an example of how to run the pipeline on our local Fred Hutch cluster, gizmo.
To run on you're own compute system: 
(1) Create a config file similar to `nextflow.gizmo.config` that describes the job submission parameters and partitions you would like to use.
(2) Create a run script using `run_gizmo_cat.sh` as a template with your own paths for temporary directories and modules you might need to load.
More instruction on creating config files can be found [here](https://www.nextflow.io/docs/latest/config.html#configuration-file)

After the pipeline finished merging all counts from the alignments,
the output of the pipeline is a file, `phip_data/pan-cov-ds.phip`
containing sample and peptide metadata tied to the raw counts matrix like so:

<p align="center">
  <img src="cartoons/Xarray.png" width="350">
</p>

With the colored columns representing coordinate dimensions and the grey squares representing data arrays
organized by respective shared dimensions.

### Running the analysis notebooks

Upon request, we can also provide you with the pre-aligned xarray dataset described above in a
[Pickle](https://docs.python.org/2/library/pickle.html) 
dumped binary file, named `pan-cov-ds.phip`. This is the file generated from the pipeline
described above, as well as the only file you'll need to run all the analysis notebooks.

We suggest using [conda](https://www.anaconda.com/) to create the environment like so:
```
conda env create -f environment.yml && conda activate pan-cov-manuscript
mkdir -p _ignore && cd _ignore
git clone git@github.com:matsengrp/phippery.git
cd phippery && python setup.py install --user && cd ../../
```
Once obtained, place the `pan-cov-ds.phip` into the `_ignore` directory then simply
```
cd analysis-notebooks
jupyter notebook
```
to launch jupyter notebooks.

These notebooks contains all analysis run for the manuscript using custom code from 
[phippery](https://github.com/matsengrp/phippery)

**Pan-CoV-Analysis.ipynb**

This notebook will run all analysis and create a diretory with all plots as well as "hits_df.csv"
which gives a row descibing each one of the peptide hits computed at the given False positive rate.
The parameters up top may easily changed to see how results change at different FPR or metrics
for ranking peptides.

**Peptide-Alignments.ipynb**

Once "hits_df.csv" has been created, you can run the alignments using this notebooks.
It's trivial to change the penalty scores to see how the results change across the library for all hits. 
