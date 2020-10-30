## Pan-CoV PhIP-Seq Manuscript Repository

This repository contains all necessary code for analysis described in Stoddard et. al. "Viral proteome-wide epitope profiling reveals binding signatures of SARS-CoV-2 immune response and cross-reactivity with endemic human coronaviruses".
The repository contains the code for a [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) pipeline,
dubbed [phip-flow](https://github.com/matsengrp/phip-flow), that produces an xarray dataset.

Once this dataset is obtained, the jupyter notebooks containing all analysis may be run using the environment
provided here along with the [phippery](git@github.com:matsengrp/phippery.git) package installed.

### Abstract

A major goal of current SARS-CoV-2 vaccine efforts is to elicit antibody responses that confer protection.
Mapping the epitope targets of the SARS-CoV-2 antibody response is critical for innovative vaccine design, diagnostics, and development of therapeutics. 
Here, we developed a phage display library to map antibody binding sites at high resolution within the complete viral proteomes of all human-infecting coronaviruses in patients with mild or moderate/severe COVID-19.
The dominant immune responses to SARS-CoV-2 were targeted to regions spanning the Spike protein, Nucleocapsid, and ORF1ab.
Some epitopes were seen in the majority of individuals while others were rare. Further, we found variation in the number of epitopes targeted by different individuals.
We also identified a set of cross-reactive sequences that were bound by antibodies in a subset of SARS-CoV-2 negative individuals. 
Finally, we uncovered a subset of enriched epitopes from commonly circulating human coronaviruses with significant homology to highly reactive SARS-CoV-2 sequences.

### Analysis environment

The library design script and analysis notebooks contain all code needed to run analysis 
for the manuscript using custom code from 
[phippery](https://github.com/matsengrp/phippery) along with a few other popular python packages.

We suggest using [conda](https://www.anaconda.com/) to create the environment like so:
```
conda env create -f environment.yml && conda activate pan-cov-manuscript
mkdir -p _ignore && cd _ignore
git clone git@github.com:matsengrp/phippery.git
cd phippery && python setup.py install --user && cd ../../
```

### library design

The pan-human CoV was created using a script that can also be found 
[here](https://github.com/jbloomlab/phipseq_oligodesign) 
and was written by Kate H.D. Crawford in the Bloom lab.
The fasta files needed to create the library are located in the library-design directory. 
Simply run

```
cd library-design
tar -xvf prot_files.tar
```

Next, to create the library, make sure your environment is activated then

```
python phip_seq_oligodesign.py prot_files out_dir protein
```

### Running the alignment pipeline

All raw fastq files for the samples can be obtained upon request and will be provided in a tarball, `NGS.tar`.
Place this file in the `alignment-pipeline/` directory and run `tar -xvf NGS.tar` to produce the `NGS` directory.

The pipeline can be run anywhere with `Nextflow` and `Docker` (or `Singularity`) installed. 
The scripts inside the directory provide an example of how to run the pipeline on our local Fred Hutch cluster, gizmo.
To run on you're own compute system: 
(1) Create a config file similar to `nextflow.gizmo.config` that describes the job submission parameters and partitions you would like to use.
(2) Create a run script using `run_gizmo_cat.sh` as a template with your own paths for temporary directories and modules you might need to load.
More instruction on creating config files can be found [here](https://www.nextflow.io/docs/latest/config.html#configuration-file).

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

Once obtained, place the `pan-cov-ds.phip` into the `_ignore` directory then simply
```
cd analysis-notebooks
jupyter notebook
```
to launch jupyter notebooks.


**Pan-CoV-Analysis.ipynb**

This notebook will run all analysis and create a diretory with all plots as well as `hits_df.csv`
which gives a row describing each one of the peptide hits computed at the given false positive rate (FPR).
The parameters up top may easily changed to see how results change at different FPR or metrics
for ranking peptides.

**Peptide-Alignments.ipynb**

Once `hits_df.csv` has been created, you can run the local alignments using this notebook.
It's trivial to change the penalty scores to see how the results change across the library for all hits. 
