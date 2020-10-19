import phippery
#from phippery.normalize import rank_data_table, standardized_enrichment, cpm, size_factors
from phippery.utils import sample_id_coordinate_subset, peptide_id_coordinate_subset
#from phippery.tidy import tidy_ds
#from phippery.modeling import gamma_poisson_model

import xarray as xr
#import numpy as np
#import pandas as pd
import pickle
#import matplotlib.pyplot as plt
#from plotnine import *
#%matplotlib inline
#%config InlineBackend.figure_format = 'retina'

#import os
import sys
#import copy
#from functools import reduce
#from collections import defaultdict

import 

#import warnings
#warnings.filterwarnings('ignore')

# load the raw counts dataset
ds = pickle.load(open(sys.argv[1],"rb"))

# we've decided LIB3 == MEGSUB
batch_col = ds.sample_table.loc[:,"library_batch"]
batch_col[batch_col == "LIB3"] = "MEGSUB"

# change some names
prot_col = ds.peptide_table.loc[:,"Protein"]

name_map = {
    "ORF3a_protein" : "ORF3a",
    "ORF6_protein" : "ORF6", 
    "ORF7a_protein" : "ORF7a",
    "ORF8_protein" : "ORF8",
    "replicase_polyprotein_1ab" : "replicase",
    "envelope" : "Envelope",
    "membrane" : "Membrane",
    "nucleocapsid" : "Nucleocapsid",
    "replicase" : "ORF1ab",
    "spike" : "Spike"
}

for old_name, new_name in name_map.items(): 
    prot_col[prot_col == old_name] = new_name

# throw out HIV positive and temporal measures of the same samples.
non_HIV = sample_id_coordinate_subset(ds, "patient_status", is_not_equal_to="hiv pos control")
non_60 = sample_id_coordinate_subset(ds, "patient_status", is_not_equal_to="conv outpatient 60d")
non_90 = sample_id_coordinate_subset(ds, "patient_status", is_not_equal_to="conv outpatient 90d")
non_hiv_neg = sample_id_coordinate_subset(ds, "patient_status", is_not_equal_to="HIV neg control")
samples_to_keep = set.intersection(set(non_HIV), set(non_60), set(non_90), set(non_hiv_neg))

# we're only keeping a single strain per virus.
strains_to_keep = {
    '229E':'SC0865',
    'HIV1':'BG505',
    'HKU1':'Caen1',
    'MERS,':'KFMC4',
    'NL63,':'ChinaGD01',
    'OC43,':'SC0776',
    'SARS,':'Urbani',
    'SARSCoV2,':'WuhanHu1',
    'batSL,':'CoVZC45'
}

# throw out 
peptides_to_keep = peptide_id_coordinate_subset(ds, where="Strain", is_in=list(strains_to_keep.values()))
ds = ds.loc[dict(sample_id=list(samples_to_keep), peptide_id=peptides_to_keep)]
pickle.dump(ds, open(sys.argv[2], "wb"))
