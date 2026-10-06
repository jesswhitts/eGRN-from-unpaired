#!/usr/bin/env python
# coding: utf-8

import os
import mudata
import re
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple
from typing import Optional, Union
import sklearn
from scipy import sparse
from itertools import combinations
import scenicplus.scenicplus_class
from scenicplus.utils import p_adjust_bh, flatten_list
from scenicplus.scenicplus_class import mudata_to_scenicplus
from scenicplus.eregulon_enrichment import get_eRegulons_as_signatures


# Filepaths and Filtering Threshold
scplus_mudata = './scplusmdata.h5mu'
metadata = '../../scATAC/full_metadata.csv'
mdata_metadata = './scplus_mdata_obs.csv'
ctx_h5 = './ctx_results.hdf5'
dem_h5 = './dem_results.hdf5'
filt_mudata = './filtered_scplusmdata.h5mu'
correlation_threshold = 0.5



# Import MuData
scplus_mdata = mudata.read(scplus_mudata)

# Load metadata
obs = pd.read_csv(metadata, index_col=0)
sample_obs = obs[['sample', 'Lambo_et_al_ID', 'Patient_Sample', 'Age_Months', 
                  'Disease_free_days', 'Clinical_Blast_Percent',
                  'Expected_Driving_Aberration', 'Subgroup']].set_index('sample').drop_duplicates()

# RNA Annotations
scplus_mdata.obs['scRNA_counts:Classified_Celltype'] = scplus_mdata.obs.index
scplus_mdata.obs['scRNA_counts:Classified_Celltype'] = scplus_mdata.obs['scRNA_counts:Classified_Celltype'].astype('str')
scplus_mdata.obs['scRNA_counts:Celltype'] = scplus_mdata.obs['scRNA_counts:Classified_Celltype']
scplus_mdata.obs['scRNA_counts:Malignant'] = scplus_mdata.obs['scRNA_counts:Classified_Celltype'] 
scplus_mdata.obs['scRNA_counts:Sample'] = scplus_mdata.obs['scRNA_counts:Classified_Celltype'] 
scplus_mdata.obs['scRNA_counts:Donor'] = scplus_mdata.obs['scRNA_counts:Classified_Celltype'] 

for row in range(0,len(scplus_mdata.obs['scRNA_counts:Classified_Celltype']),1):
    celltype = scplus_mdata.obs['scRNA_counts:Classified_Celltype'][row]
    splitStr = celltype.split('_')
    # Update Celltypes
    scplus_mdata.obs['scRNA_counts:Malignant'][row] = splitStr[0]
    scplus_mdata.obs['scRNA_counts:Celltype'][row] = splitStr[1]
    scplus_mdata.obs['scRNA_counts:Classified_Celltype'][row] = "_".join(splitStr[:-1])
    # Annotate Sample/Donor
    sample = splitStr[2]
    match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", sample, re.I)
    if match:
        items = list(match.groups())
    donor = str(items[0])+str(items[1])
    scplus_mdata.obs['scRNA_counts:Sample'][row] = sample
    scplus_mdata.obs['scRNA_counts:Donor'][row] = donor

# Add Relevant Metadata Columns to Assay
for col in sample_obs.columns:
    # For numeric cols
    if np.issubdtype(sample_obs[col].dtype, np.number):
        scplus_mdata.obs['scRNA_counts:' + col] = 0.0
    # For strings
    else:
        scplus_mdata.obs['scRNA_counts:' + col] = 'Unknown'

# Update DataFrame
for samp in list(scplus_mdata.obs['scRNA_counts:Sample'].unique()):
    cells = scplus_mdata[scplus_mdata.obs['scRNA_counts:Sample'] == samp].obs.index
    for col in sample_obs.columns:
        scplus_mdata.obs.loc[cells, 'scRNA_counts:' + col] = sample_obs.loc[samp][col]


# ATAC Annotations
scplus_mdata.obs['scATAC_counts:Classified_Celltype'] = scplus_mdata.obs.index
scplus_mdata.obs['scATAC_counts:Classified_Celltype'] = scplus_mdata.obs['scATAC_counts:Classified_Celltype'].astype('str')
scplus_mdata.obs['scATAC_counts:Celltype'] = scplus_mdata.obs['scATAC_counts:Classified_Celltype']
scplus_mdata.obs['scATAC_counts:Malignant'] = scplus_mdata.obs['scATAC_counts:Classified_Celltype'] 
scplus_mdata.obs['scATAC_counts:Sample'] = scplus_mdata.obs['scATAC_counts:Classified_Celltype'] 
scplus_mdata.obs['scATAC_counts:Donor'] = scplus_mdata.obs['scATAC_counts:Classified_Celltype'] 

for row in range(0,len(scplus_mdata.obs['scATAC_counts:Classified_Celltype']),1):
    celltype = scplus_mdata.obs['scATAC_counts:Classified_Celltype'][row]
    splitStr = celltype.split('_')
    # Update Celltypes
    scplus_mdata.obs['scATAC_counts:Malignant'][row] = splitStr[0]
    scplus_mdata.obs['scATAC_counts:Celltype'][row] = splitStr[1]
    scplus_mdata.obs['scATAC_counts:Classified_Celltype'][row] = "_".join(splitStr[:-1])
    # Annotate Sample/Donor
    sample = splitStr[2]
    match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", sample, re.I)
    if match:
        items = list(match.groups())
    donor = str(items[0])+str(items[1])
    scplus_mdata.obs['scATAC_counts:Sample'][row] = sample
    scplus_mdata.obs['scATAC_counts:Donor'][row] = donor

# Add Relevant Metadata Columns to Assay
for col in sample_obs.columns:
    # For numeric cols
    if np.issubdtype(sample_obs[col].dtype, np.number):
        scplus_mdata.obs['scATAC_counts:' + col] = 0.0
    # For strings
    else:
        scplus_mdata.obs['scATAC_counts:' + col] = 'Unknown'
# Update DataFrame
for samp in list(scplus_mdata.obs['scATAC_counts:Sample'].unique()):
    cells = scplus_mdata[scplus_mdata.obs['scATAC_counts:Sample'] == samp].obs.index
    for col in sample_obs.columns:
        scplus_mdata.obs.loc[cells, 'scATAC_counts:' + col] = sample_obs.loc[samp][col]

# Save MuData Obs
obs = scplus_mdata.obs
obs.to_csv(mdata_metadata)


# Make Unfiltered Regulon GMT
subset_direct_metadata = scplus_mdata.uns['direct_e_regulon_metadata'][['Gene', 'eRegulon_name','Gene_signature_name',
                                          'Region_signature_name']].set_index('eRegulon_name').drop_duplicates()
subset_extended_metadata = scplus_mdata.uns['extended_e_regulon_metadata'][['Gene', 'eRegulon_name','Gene_signature_name',
                                          'Region_signature_name']].set_index('eRegulon_name').drop_duplicates()
eregulons = pd.concat([subset_direct_metadata,subset_extended_metadata])
eregulons.to_csv('eRegulons_metadata_unfiltered.csv')
regulons = list(eregulons.index.unique())

stable_regulons = None
for r in regulons:
    subset = eregulons.loc[r]
    gene_list = list(eregulons.loc[r]['Gene'])
    regulon = pd.DataFrame(gene_list, columns=[r])
    stable_regulons = pd.concat([
        stable_regulons, pd.DataFrame({r: gene_list})
    ], axis=1)

# List of target genes per stable regulons
stable_regulons.to_csv('scplus_regulons_unfiltered.csv', index=False)

# Make GMT for AUCell
with open("scplus_regulons_unfiltered.gmt", 'wt') as gmt:
    for c in stable_regulons.columns:
        genes = '\t'.join(stable_regulons.loc[:, c].dropna())
        gmt.write(f'{c}\tRegulon_{c}\t{genes}\n')


# # Regulon Filtering
# Correlation between region based regulons and gene based regulons
direct_df = pd.DataFrame(scplus_mdata["direct_gene_based_AUC"].X, index=scplus_mdata["direct_gene_based_AUC"].obs_names,
                         columns=scplus_mdata["direct_gene_based_AUC"].var_names)
extended_df = pd.DataFrame(scplus_mdata["extended_gene_based_AUC"].X, index=scplus_mdata["extended_gene_based_AUC"].obs_names,
                           columns=scplus_mdata["extended_gene_based_AUC"].var_names)
df1 = direct_df.join(extended_df)

direct_df = pd.DataFrame(scplus_mdata["direct_region_based_AUC"].X, index=scplus_mdata["direct_region_based_AUC"].obs_names,
                         columns=scplus_mdata["direct_region_based_AUC"].var_names)
extended_df = pd.DataFrame(scplus_mdata["extended_region_based_AUC"].X, index=scplus_mdata["extended_region_based_AUC"].obs_names,
                           columns=scplus_mdata["extended_region_based_AUC"].var_names)
df2 = direct_df.join(extended_df)
regulons = list(df1.columns) + list(df2.columns)
print("Total number of regulons detected: " + str(len(regulons)))

df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]
correlations = df1.corrwith(df2, axis = 0)

df = pd.DataFrame(correlations).sort_values(by=0)
df['regulon']=df.index

# Save Correlation Result
df.to_csv('eRegulon_correlations.csv')

correlations = correlations[abs(correlations) > abs(correlation_threshold)]
print("Regulons passing correlation threshold: " + str(len(correlations)))

df['selected'] = 'False'
df.loc[list(correlations.index), 'selected'] = 'True'


sns.set_style("whitegrid")
sns.scatterplot(x="regulon", y=0, data=df, hue="selected", s=10, linewidth=0)
plt.xticks([])
legend = plt.legend(title="",
                    fontsize=10,
                    loc="center left",
                    bbox_to_anchor=(1, 0, 0.5, 1),
                    frameon=False)
plt.xlabel('Regulon', fontsize=12)
plt.ylabel('GEX/ATAC Correlation Scores', fontsize=12)
plt.savefig('gex_atac_correlation_scores.png')


# In[ ]:


# Keep only R2G +
keep = [x for x in correlations.index if '+/+' in x] #+ [x for x in correlations.index if '-/+' in x]
print("Valid R2G+ Regulons: " + str(len(keep)))

# Keep extended if not direct
extended = [x for x in keep if 'extended' in x]
print("Extended Regulons: " + str(len(extended)))
direct = [x for x in keep if not 'extended' in x]
print("Direct Regulons: " + str(len(direct)))
keep_extended = [x for x in extended if not x.replace('extended_', 'direct_') in direct]
keep = direct + keep_extended
print("Regulons to Keep: " + str(len(keep)))

# Filter regulons with fewer than 10 genes
keep_gene = [x for x in scplus_mdata['direct_gene_based_AUC'].var_names if x.split('_(')[0] in keep]
keep_direct = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
keep_gene = [x for x in scplus_mdata['extended_gene_based_AUC'].var_names if x.split('_(')[0] in keep]
keep_extended = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
keep_genes = keep_direct + keep_extended

keep_direct_region = [x for x in scplus_mdata['direct_region_based_AUC'].var_names if x.split('_(')[0] in keep]
keep_extended_region = [x for x in scplus_mdata['extended_region_based_AUC'].var_names if x.split('_(')[0] in keep]
keep_regions = keep_direct_region + keep_extended_region

keep_genes_regions = keep_genes + keep_regions

low_qual = [x for x in regulons if x not in keep_genes_regions]
print("Filtering Out " + str(len(low_qual)) + " Low Quality Regulons")


# In[ ]:


# Save Valid Regulons
direct_metadata = scplus_mdata.uns['direct_e_regulon_metadata']
extended_metadata = scplus_mdata.uns['extended_e_regulon_metadata']

subset = np.array([s in low_qual for s in direct_metadata.Gene_signature_name])
direct_metadata = direct_metadata[~subset]

subset = np.array([s in low_qual for s in extended_metadata.Gene_signature_name])
extended_metadata = extended_metadata[~subset]

scplus_mdata = scplus_mdata[:, ~scplus_mdata.var_names.isin(low_qual)].copy()
scplus_mdata.uns['direct_e_regulon_metadata'] = direct_metadata 
scplus_mdata.uns['extended_e_regulon_metadata'] = extended_metadata

pickle.dump(direct_metadata,
            open(os.path.join('direct_eRegulon_metadata.pkl'), 'wb'))
pickle.dump(extended_metadata,
            open(os.path.join('extended_eRegulon_metadata.pkl'), 'wb'))

scplus_mdata.write(filt_mudata)


# Make Filtered Regulon GMT
subset_direct_metadata = scplus_mdata.uns['direct_e_regulon_metadata'][['Gene', 'eRegulon_name','Gene_signature_name',
                                          'Region_signature_name']].set_index('eRegulon_name').drop_duplicates()
subset_extended_metadata = scplus_mdata.uns['extended_e_regulon_metadata'][['Gene', 'eRegulon_name','Gene_signature_name',
                                          'Region_signature_name']].set_index('eRegulon_name').drop_duplicates()
eregulons = pd.concat([subset_direct_metadata,subset_extended_metadata])
eregulons.to_csv('eRegulons_metadata_filtered.csv')
regulons = list(eregulons.index.unique())

stable_regulons = None
for r in regulons:
    subset = eregulons.loc[r]
    gene_list = list(eregulons.loc[r]['Gene'])
    regulon = pd.DataFrame(gene_list, columns=[r])
    stable_regulons = pd.concat([
        stable_regulons, pd.DataFrame({r: gene_list})
    ], axis=1)

# List of target genes per stable regulons
stable_regulons.to_csv('scplus_regulons_filtered.csv', index=False)

# Make GMT for AUCell
with open("scplus_regulons_filtered.gmt", 'wt') as gmt:
    for c in stable_regulons.columns:
        genes = '\t'.join(stable_regulons.loc[:, c].dropna())
        gmt.write(f'{c}\tRegulon_{c}\t{genes}\n')

# In[ ]:


# # Create SCENIC+ Obj
# Create SCENIC+ for Each Patient (Useful for Pertubation Modelling)
donors = list(scplus_mdata.obs['scRNA_counts:Donor'].value_counts().index)

for d in donors:
    # Subset for Specific Donor
    subset_mdata = scplus_mdata[scplus_mdata.obs['scRNA_counts:Donor']==d].copy()
    subset_mdata.write(os.path.join(d + '_scplus_mdata.h5mu'))
    donor_obs = subset_mdata.obs
    subset_mdata.uns['direct_e_regulon_metadata'] = direct_metadata 
    subset_mdata.uns['extended_e_regulon_metadata'] = extended_metadata
    # Convert to SCPlus
    scplus_obj = mudata_to_scenicplus(subset_mdata)
    scplus_obj.uns["eRegulon_signatures"] = get_eRegulons_as_signatures(scplus_obj.uns["eRegulon_metadata"])
    # Fix Metadata
    scplus_obj.metadata_cell = donor_obs
    pickle.dump(scplus_obj, open(os.path.join(d + '_scplus_obj.pkl'), 'wb'))


# In[ ]:


# Create Full SCENIC+ Obj for Downstream Functions
scplus_obj = mudata_to_scenicplus(
    scplus_mdata,
    path_to_cistarget_h5 = ctx_h5,
    path_to_dem_h5 = dem_h5 
)
scplus_obj.uns["eRegulon_signatures"] = get_eRegulons_as_signatures(scplus_obj.uns["eRegulon_metadata"])

# Fix Metadata
scplus_obj.metadata_cell = obs

# Save
pickle.dump(scplus_obj,
            open(os.path.join('scplus_obj.pkl'), 'wb'))


# Get Metacell Counts

# Get RNA Counts per Metacell
rna = scplus_mdata['scRNA_counts']
counts = pd.DataFrame(rna.X, index=rna.obs_names, columns=rna.var_names)
rna = ad.AnnData(counts)
rna.obs = scplus_mdata.obs
rna.X = sparse.csr_matrix(rna.X)
rna.write_h5ad(os.path.join('scplus_metacells_rna.h5ad'))

# Get ATAC Counts per Metacell
atac = scplus_mdata['scATAC_counts']
counts = pd.DataFrame(atac.X, index=atac.obs_names, columns=atac.var_names)
atac = ad.AnnData(counts)
atac.obs = scplus_mdata.obs
atac.X = sparse.csr_matrix(atac.X)
atac.var['ID'] = atac.var.index
atac.var[['chr', 'region']] = atac.var['ID'].str.split(':', expand=True)
atac.var[['start', 'end']] = atac.var['region'].str.split('-', expand=True)
atac.write_h5ad(os.path.join('scplus_metacells_atac.h5ad'))
