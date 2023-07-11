#!/usr/bin/env python
# coding: utf-8

# ## Apply an image-based profiling pipeline using pycytominer
# 
# As described fully in [Caicedo et al. 2017](https://doi.org/10.1038/nmeth.4397), an image-based profiling pipeline consists of three core steps:
# 
# 1. Aggregation
# 2. Normalization
# 3. Feature selection
# 
# [Pycytominer](https://github.com/cytomining/pycytominer) is a python package, built on top of pandas, that facilitates all of these steps and more.
# 
# ### Data levels
# 
# The concept of "data levels" is important to understand when implementing an image-based profiling pipeline.
# 
# | Data | Level |
# | :---- | :---- |
# | Images | Level 1 |
# | Single cell profiles (SQLite) | Level 2 |
# | Aggregated profiles with metadata information | Level 3 |
# | Normalized aggregated profiles | Level 4a |
# | Normalized and feature selected profiles | Level 4b |
# | Consensus profiles | Level 5 |

# In[1]:


import pathlib
import pandas as pd

from pycytominer import annotate, normalize, feature_select, consensus
from pycytominer.cyto_utils import cells


# ### Step 0 - Initialize data and options

# In[2]:


data_dir = "data"
plate_id = "218360"
platemap = f"{plate_id}.txt"

sqlite_file = f"sqlite:///{data_dir}/{plate_id}.sqlite"
sqlite_file


# In[3]:


output_dir = pathlib.Path("profiles")
output_dir.mkdir(exist_ok=True)


# In[4]:


# Load platemap file
platemap_file = pathlib.Path(f"{data_dir}/{platemap}")
platemap_df = pd.read_csv(platemap_file, sep="\t")

print(platemap_df.shape)
platemap_df.head(3)


# In[5]:


# Output options
compression_options = {"method": "gzip", "mtime": 1}


# ### Step 1 - Single cell aggregation
# 
# In this step, **level 2 profiles** (single cells) are processed to **level 3 profiles** (well-level profiles).

# In[6]:


# Initialize the single cells class
sc_data = cells.SingleCells(
    sqlite_file,
    strata=["Metadata_Plate", "Metadata_Well"],
    features="infer",
    compartments=["cells", "cytoplasm", "nuclei"],
    merge_cols=["TableNumber", "ImageNumber"],
    aggregation_operation="median",
    load_image_data=True,
    subsample_frac=1,
    subsample_n="all",
    subsampling_random_state="none"
)


# In[7]:


# Count cells
cell_count_df = sc_data.count_cells()

cell_count_df = (
    cell_count_df.merge(
        platemap_df,
        left_on="Metadata_Well",
        right_on="well_position"
    )
)

# Output cell count to file
cell_count_file = pathlib.Path(f"{output_dir}/cell_counts.tsv")
cell_count_df.to_csv(cell_count_file, sep="\t", index=False)

print(cell_count_df.shape)
cell_count_df.head()


# In[8]:


# Perform the aggregation - output well level profiles
output_file = pathlib.Path(f"{output_dir}/{plate_id}.csv.gz")

# Aggregate profiles can output a file, or save the result to a variable
# Here, we output the intermediate result to a file
sc_data.aggregate_profiles(
    output_file=output_file,
    compute_subsample=True,
    compression_options=compression_options,
    float_format=None
)


# In[9]:


# Read in and preview what was output in the previous step
# Note: This is not necessary to perform in automated pipelines, but is a nice check.
aggregated_df = pd.read_csv(output_file)

print(aggregated_df.shape)
aggregated_df.head(3)


# ## Step 2 - Annotate wells using the platemap file
# 
# In this step, **level 3 profiles** (well-level profiles) are annotated with platemap metadata.

# In[10]:


# Annotate profiles
annotate_file = pathlib.Path(f"{output_dir}/{plate_id}_augmented.csv.gz")

annotate(
    profiles=output_file,
    platemap=platemap_df,
    join_on=["Metadata_well_position", "Metadata_Well"],
    output_file=annotate_file,
    compression_options=compression_options,
)


# In[11]:


# Read in and preview what was output in the previous step 
annotated_df = pd.read_csv(annotate_file)

print(annotated_df.shape)
annotated_df.head(3)


# ## Step 3 - Normalize well-level profiles
# 
# In this step, **level 3 profiles** (well-level profiles) are normalized to form **level 4a profiles**.

# In[12]:


# Normalize profiles
normalize_file = pathlib.Path(f"{output_dir}/{plate_id}_normalized.csv.gz")

normalize(
    profiles=annotate_file,
    features="infer",
    meta_features="infer",
    samples="Metadata_treatment == '0.1% DMSO'",
    method="standardize",
    output_file=normalize_file,
    compression_options=compression_options,
)


# In[13]:


# Read in and preview what was output in the previous step 
normalized_df = pd.read_csv(normalize_file)

print(normalized_df.shape)
normalized_df.head(3)


# ## Step 4 - Apply feature selection
# 
# In this step, we apply a series of feature selection steps to **level 4a profiles** (normalized well-level profiles) to form **level 4b profiles**.

# In[14]:


# Apply feature selection
feature_select_file = pathlib.Path(f"{output_dir}/{plate_id}_normalized_feature_select.csv.gz")

feature_select_opts = [
    "variance_threshold",
    "drop_na_columns",
    "correlation_threshold",
    "blocklist",
    "drop_outliers"
]

feature_select(
    profiles=normalize_file,
    features="infer",
    samples="all",
    operation=feature_select_opts,
    output_file=feature_select_file,
    compression_options=compression_options,
)


# In[15]:


# Read in and preview what was output in the previous step 
feature_select_df = pd.read_csv(feature_select_file)

print(feature_select_df.shape)
feature_select_df.head(3)


# ## Step 5 - Form consensus signatures
# 
# In this step, we collapse replicates (**level 4 profiles**) into a single profile.
# This forms **level 5 profiles**.

# In[16]:


# Generate consensus profiles
consensus_file = pathlib.Path(f"{output_dir}/{plate_id}_consensus.csv.gz")

consensus(
    profiles=feature_select_df,
    replicate_columns=["Metadata_clone_number", "Metadata_treatment"],
    operation="modz",
    features="infer",
    output_file=consensus_file,
    modz_args={
        "method": "spearman",
        "min_weight": 0.01,
        "precision": 4,
    },
    compression_options=compression_options,
)


# In[17]:


# Read in and preview what was output in the previous step 
consensus_df = pd.read_csv(consensus_file)

print(consensus_df.shape)
consensus_df.head(3)

