#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow \
    run main.nf \
    -c nextflow.gizmo.config \
    --params_file config.json \
    -with-report output/nextflow_report.html \
    -work-dir /fh/scratch/delete30/matsen_e/jgallowa/temp/work/ \
    -ansi-log false \
    -resume
