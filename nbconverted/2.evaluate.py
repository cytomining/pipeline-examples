#!/usr/bin/env python
# coding: utf-8

# ## Evaluate profile quality with cytominer-eval
# 
# A necessary step in every image-base profiling pipeline is to evaluate profile quality.
# 
# [cytominer-eval](https://github.com/cytomining/cytominer-eval) is a python package, built on top of pandas with a purpose of quickly evaluating profile quality.
# 
# Currently, there are three metrics enabled:
# 
# 1. Replicate reproducibility - a measurment of how many profiles are above a null distribution threshold
# 2. Enrichment - Fisher Exact test at a specific enrichment percentile testing if replicates are more similar than non-replicates
# 3. Grit - a measurement of phenotype consistency and strength, as compared to a negative control
# 
# Typically, we evaluate quality on a per-plate basis, but this is not a strict requirement.

# In[1]:


import pathlib
import pandas as pd
import plotnine as gg

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features

from cytominer_eval import evaluate


# In[2]:


# Load level 4b profiles derived from 1.profile.ipynb
profile_file = pathlib.Path("profiles/218360_normalized_feature_select.csv.gz")
profile_df = pd.read_csv(profile_file)

print(profile_df.shape)
profile_df.head(3)


# In[3]:


# Distinguish metadata and morphology features
meta_features = infer_cp_features(profile_df, metadata=True)
morphology_features = infer_cp_features(profile_df)

len(morphology_features)


# ## Set a variable that describes which columns to use for replicate information
# 
# In our case, two columns together tell us which profiles are replicates.

# In[4]:


replicate_groups = ["Metadata_clone_number", "Metadata_treatment"]


# ### Calculate replicate reproducibility

# In[5]:


percent_reproducible = evaluate(
    profiles=profile_df,
    features=morphology_features,
    meta_features=meta_features,
    replicate_groups=replicate_groups,
    operation="replicate_reproducibility",
    similarity_metric="pearson",
)

percent_reproducible


# In[6]:


enrichment = evaluate(
    profiles=profile_df,
    features=morphology_features,
    meta_features=meta_features,
    replicate_groups=replicate_groups,
    operation="enrichment",
    similarity_metric="pearson",
    enrichment_percentile=[0.99, 0.95, 0.9, 0.75, 0.5]
)

enrichment


# ### Calculate replicate enrichment

# ## Calculate grit
# 
# In this plate, it doesn't make sense to calculate grit.
# Instead, we use all CRISPR perturbations from the [Cell Health](https://github.com/broadinstitute/cell-health) experiment.
# 
# ### Grit intuition
# 
# Grit can be calculated in the scenario where there is some sort of replicate hierarchy.
# 
# For example, grit can be calculated in the following scenarios:
# 
# | Scenario | Group ID | Replicate ID |
# | :------- | :------- | :----------- |
# | CRISPR perturbation | Gene target | Guide ID |
# | Compound treatment | Mechanism of action | Compound ID |
# | Single cell (scGrit) | Replicate ID | Single cell profile |
# 
# We can interpret grit to mean _phenotype strength_.
# A high grit score indicates that the replicate ID is a) similar to other perturbations in the group; b) different compared to a non-targeting control.
# Refer to https://github.com/broadinstitute/grit-benchmark for more details.

# In[7]:


# Load Cell Health data
commit = "07e4b40c39dd27084be36fbef4d64c5654b2960f"
base_url = f"https://github.com/broadinstitute/cell-health/raw/{commit}"
url = f"{base_url}/1.generate-profiles/data/processed/cell_health_profiles_merged.tsv.gz"

df = pd.read_csv(url, sep="\t")

print(df.shape)
df.head(2)


# In[8]:


# Perform feature selection on all Cell Health profiles
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
    "drop_outliers",
]

df = feature_select(
    profiles=df,
    operation=feature_select_ops,
    na_cutoff=0
)

features = infer_cp_features(df)
meta_features = infer_cp_features(df, metadata=True)


# In[9]:


# Define group and replicate IDs
barcode_col = "Metadata_pert_name"
gene_col = "Metadata_gene_name"

replicate_group_grit = {
    "profile_col": barcode_col,
    "replicate_group_col": gene_col
}

control_group_cut = ["Chr2", "Luc", "LacZ"]

control_barcodes_cut = (
    df.loc[
        df[replicate_group_grit["replicate_group_col"]].isin(control_group_cut),
        replicate_group_grit["profile_col"]
    ]
    .unique()
    .tolist()
)

control_barcodes_cut


# In[10]:


# Calculate grit per cell line
grit_results = []
for cell_line in df.Metadata_cell_line.unique():
    result = evaluate(
        profiles=df.query("Metadata_cell_line == @cell_line"),
        features=features,
        meta_features=[barcode_col, gene_col],
        replicate_groups=replicate_group_grit,
        operation="grit",
        grit_control_perts=control_barcodes_cut
    )
    
    result = result.assign(cell_line=cell_line)

    grit_results.append(result)

# Merge results
grit_results = (
    pd.concat(grit_results)
    .sort_values(by="grit", ascending=False)
    .reset_index(drop=True)
    .dropna()
)

print(grit_results.shape)
grit_results.head()


# In[11]:


# Visualize grit results
gg.options.figure_size = (9, 6)

(
    gg.ggplot(grit_results, gg.aes(x="reorder(group, grit)", y="grit", fill="cell_line")) +
    gg.geom_point() +
    gg.theme_bw() +
    gg.theme(axis_text_x=gg.element_text(angle=90)) +
    gg.xlab("Gene") +
    gg.ylab("Grit")
)

