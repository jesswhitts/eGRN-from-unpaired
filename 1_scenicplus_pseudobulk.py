#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
##################################################### AML ATLAS: SCENIC+ ########################################################
#################################### Subset for Relevant Sample and Prepare Non-Multome Data ####################################
#################################################################################################################################

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import re
import pandas as pd
import numpy as np
import scanpy as sc
import pyranges as pr
import requests
import ray
import pickle
import pybiomart as pbm
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
from pycisTopic.iterative_peak_calling import *
from pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.cistopic_class import *
from pycisTopic.clust_vis import run_umap
from pycisTopic.clust_vis import plot_metadata
import multiprocessing
max_cpu = multiprocessing.cpu_count()

cwd = os.getcwd()
work_dir = './'
data_dir = '/data/stemcell/jwhittle/AML_scAtlas/18-lambo/'
tmp_dir = '/scratch/wsspaces/jwhittle-tmp'
macs_path='/data/stemcell/jwhittle/mambaforge/envs/scenicplus/bin/macs2'
adata_path = '/data/stemcell/jwhittle/AML_scAtlas/18-lambo/lambo_filtered_RNA.h5ad.gz'
atac_path = '/data/stemcell/jwhittle/AML_scAtlas/18-lambo/lambo_filtered_ATAC.h5ad.gz'
nr_cells_per_metacell = 10

if not os.path.exists(work_dir):
    os.makedirs(work_dir)
if not os.path.exists(os.path.join(work_dir, 'scATAC')):
    os.makedirs(os.path.join(work_dir, 'scATAC'))

# Path to blacklist regions
path_to_blacklist = '/data/stemcell/jwhittle/ref/hg38_resources/hg38-blacklist.v2.bed'


#################################################################################################################################
####################################################### Prepare GEX Data ########################################################
#################################################################################################################################

print('Preparing GEX Object')
adata = sc.read_h5ad(adata_path)

# Fill Raw Field with Raw Counts
adata.layers["normalised_counts"] = adata.X.copy()
adata.X = adata.layers["counts"].copy()
adata.raw = adata.copy()
adata.X = adata.layers["normalised_counts"]

# Check raw counts are integers
counts = adata.layers["counts"]
non_integer_mask = counts.data % 1 != 0
print(non_integer_mask.any())

# Select Relevant Subtype
adata = adata[adata.obs['Subgroup']=='KMT2A'].copy()
sc.pp.filter_genes(adata,min_cells=10)
print(adata)

# Select Relevant Samples for Analysis (Samples with Both RNA and ATAC Modalities)
atac = sc.read_h5ad(atac_path)
rna_samples = list(adata.obs['sample'].unique())
atac_samples = list(atac.obs['sample'].unique())
sample_list = list(set(rna_samples) & set(atac_samples))

subset = np.array([s in sample_list for s in adata.obs['sample']])
adata = adata[subset].copy()
print(adata.obs["sample"].value_counts())

adata.obs['Classified_Celltype'] = adata.obs['Malignant'].astype(
    'str') + '_' + adata.obs['Classified_Celltype'].astype(
    'str') + '_' + adata.obs['sample'].astype('str')

sample_list = list(adata.obs['Library_ID'].unique())

# Write full RNA modality
adata.write_h5ad(work_dir + 'Lambo_RNA_filtered.h5ad')

#################################################################################################################################
###################################################### Prepare ATAC Data ########################################################
#################################################################################################################################

print('Preparing Fragment Files')
# Pre-Processing of Fragment Files

# Collect Fragment Files for Processing
fragments = []
for s in sample_list:
    for (root, dirs, file) in os.walk(data_dir):
        for f in file:
            filename = os.path.join(root, f)
            if s.casefold() in filename.casefold():
                if '_filtered_fragments.tsv.gz' in filename:
                    fragments.append(filename)
                    
# Check there is one fragment file per sample
if len(sample_list)!=len(fragments):
    raise Exception("Sample and Fragment Lists are Not Equal")

# Prepare Fragment Files
samples = []

for f in fragments:
    # Add Sample Name to List
    s = os.path.basename(f)[11:]
    s = s.replace('_filtered_fragments.tsv.gz', '')
    s = s.replace('_', '')
    samples.append(s)
    # Fragment file processing (if required)
    basefile = f.replace('.tsv.gz', '')
    basefile = os.path.join(data_dir, basefile)
    if os.path.isfile(basefile + '_bedSorted.tsv.gz.tbi') == True:
        print(s + ' Fragment File Processing Not Required')
    else:
        # Add Cell Barcode Prefix
        fragment_file = pd.read_csv(basefile + '.tsv.gz', sep='\t', index_col=0)
        suffix = "-" + s.replace('_', '')
        fragment_file['RG'] = fragment_file['RG'] + suffix
        fragment_file.to_csv(basefile + '_suffix.tsv.gz', sep='\t', header=False, compression='gzip')

        # Command Line Processing (sort, bgzip and tabix indexing)
        os.system('gunzip ' + basefile + '_suffix.tsv.gz')
        os.system('cut -f1,2,3,6 ' + basefile + '_suffix.tsv > ' + basefile + '_filt_suffix.tsv')
        os.system('sortBed -i ' + basefile + '_filt_suffix.tsv > ' + basefile + '_bedSorted.tsv')
        os.system('rm ' + basefile + '_suffix.tsv')
        os.system('rm ' + basefile + '_filt_suffix.tsv')
        os.system('bgzip ' + basefile + '_bedSorted.tsv')
        os.system('tabix -p bed ' + basefile + '_bedSorted.tsv.gz')
        print(s + ' Fragment File Processing Complete')


print('Collecting Analysis Files')
# Collect Analysis Files
fragments = []
metadata = []
for s in sample_list:
    for (root, dirs, file) in os.walk(data_dir):
        for f in file:
            filename = os.path.join(root, f)
            if s.casefold() in filename.casefold():
                if '_filtered_fragments_bedSorted.tsv.gz' in filename:
                    fragments.append(filename)
                if '_filtered_metadata.tsv.gz' in filename:
                    metadata.append(filename)

fragments = list(filter(lambda x: not re.search(r'.tbi', x), fragments))

# Create Fragment File Dictionary
fragments_dict = {}

for s in samples:
    print(s)
    temp = ''
    for chr in s:
        if chr.isdigit():
            temp = temp + chr
    samp = 'AML' + temp
    samp = samp + "_" + s[len(samp): ]
    samp_fragments = list(filter(lambda x: re.search(samp, x), fragments))
    if len(samp_fragments)<2:
        samp_fragments = samp_fragments[0]
    else:
        raise Exception("Sample has multiple fragment files")
    fragments_dict[s] = os.path.join(data_dir, samp_fragments)

print(fragments_dict)

# Combine annotations
full_metadata = pd.DataFrame()

for m in metadata:
    cell_data = pd.read_csv(os.path.join(data_dir, m), sep='\t', index_col=0)
    cell_data['sample'] = cell_data["Library_ID"]
    for row in range(0,len(cell_data["sample"]),1):
        cell_data["sample"][row] = cell_data["sample"][row].replace('_', '')
    suffix = "-" + cell_data["sample"][0]
    cell_data.index = cell_data.index + suffix
    cell_data['Classified_Celltype'] = cell_data['Malignant'].astype(
        'str') + '_' + cell_data['Classified_Celltype'].astype(
        'str') + '_' + cell_data['sample'].astype('str')
    full_metadata = pd.concat([full_metadata,cell_data])


#################################################################################################################################
###################################################### Harmonise Cell Types #####################################################
#################################################################################################################################

# GEX Data
print(adata.obs["Classified_Celltype"].value_counts())
df = pd.DataFrame(adata.obs['Classified_Celltype'].value_counts())
df = df[df['Classified_Celltype'] < nr_cells_per_metacell]
remove = list(df.index.astype(str))
remove = np.array([s in remove for s in adata.obs.Classified_Celltype])
adata = adata[~remove].copy()
GEX_celltypes = list(adata.obs['Classified_Celltype'].unique())


# ATAC Data
df = pd.DataFrame(full_metadata['Classified_Celltype'].value_counts())
df = df[df['Classified_Celltype'] < nr_cells_per_metacell]
remove = list(df.index.astype(str))
remove = np.array([s in remove for s in full_metadata.Classified_Celltype])
full_metadata = full_metadata[~remove].copy()
ATAC_celltypes = list(full_metadata['Classified_Celltype'].unique())

# Common Cell Types
common_celltypes = list(set(GEX_celltypes) & set(ATAC_celltypes))

keep = np.array([s in common_celltypes for s in adata.obs.Classified_Celltype])
adata = adata[keep].copy()

keep = np.array([s in common_celltypes for s in full_metadata.Classified_Celltype])
full_metadata = full_metadata[keep].copy()

# Save Ouput
adata.write_h5ad(work_dir + 'Lambo_multimodal.h5ad')
full_metadata.to_csv(work_dir + 'scATAC/full_metadata.csv')
del(adata)
del(atac)

#################################################################################################################################
##################################################### Prepare Pseudobulk ########################################################
#################################################################################################################################

print('Generating Pseudobulk Profiles per Celltype')
# Create Pseudobulk Profiles
# Select Fragment Files for Samples Left Post Filtering
sample_keys = list(full_metadata['sample'].value_counts().index)
print(sample_keys)
fragments_dict = dict((k, fragments_dict[k]) for k in sample_keys
                      if k in fragments_dict)

# Get chromosome sizes (for hg38 here)
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]

# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)

ray.shutdown()
bw_paths, bed_paths = export_pseudobulk(input_data = full_metadata,
                 variable = 'Classified_Celltype',                                                          # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
                 sample_id_col = 'sample',
                 chromsizes = chromsizes,
                 bed_path = os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/'),  # specify where pseudobulk_bed_files should be stored
                 bigwig_path = os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bw_files/'),# specify where pseudobulk_bw_files should be stored
                 path_to_fragments = fragments_dict,                                                        # location of fragment fiels
                 n_cpu = max_cpu,                                                                                 # specify the number of cores to use, we use ray for multi processing
                 normalize_bigwig = True,
                 #remove_duplicates = True,
                 temp_dir = os.path.join(tmp_dir, 'ray_spill'),
                 split_pattern = '__')


# Save Pseudobulk Paths
pickle.dump(bed_paths,
            open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'wb'))
pickle.dump(bw_paths,
           open(os.path.join(work_dir, 'scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'wb'))


#################################################################################################################################
######################################################## Peak Calling ###########################################################
#################################################################################################################################

print("Calling Peaks per Pseudobulk Profile")
# Run peak calling
sys.stderr = open(os.devnull, "w")  # silence stderr
narrow_peaks_dict = peak_calling(macs_path,
                                 bed_paths,
                                 os.path.join(work_dir, 'scATAC/consensus_peak_calling/MACS/'),
                                 genome_size='hs',
                                 n_cpu=max_cpu,
                                 input_format='BEDPE',
                                 shift=73,
                                 ext_size=146,
                                 keep_dup = 'all',
                                 skip_empty_peaks = True,
                                 q_value = 0.05,
                                 _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
sys.stderr = sys.__stderr__  # unsilence stderr

# Merge peaks into consensus set
peak_half_width = 250
# Get consensus peaks
sys.stderr = open(os.devnull, "w")  # silence stderr
consensus_peaks=get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes, path_to_blacklist=path_to_blacklist)
sys.stderr = sys.__stderr__  # unsilence stderr

consensus_peaks.to_bed(
    path = os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed'),
    keep=True,
    compression='infer',
    chain=False)


# Create Regions Dict
path_to_regions = {}
for s in samples:
    path_to_regions[s] = os.path.join(work_dir, 'scATAC/consensus_peak_calling/consensus_regions.bed')
print(path_to_regions)

# If we want to specify which cells are valid:
#metadata_bc = pickle.load(open(os.path.join(work_dir, 'scATAC/quality_control/metadata_bc.pkl'), 'rb'))
#bc_passing_filters = pickle.load(open(os.path.join(work_dir, 'scATAC/quality_control/bc_passing_filters.pkl'), 'rb'))


#################################################################################################################################
#################################################### Create Cistopic Obj ########################################################
#################################################################################################################################

print('Creating Cistopic Object')
# Create Object
cistopic_obj_list=[create_cistopic_object_from_fragments(path_to_fragments=fragments_dict[key],
                                               path_to_regions=path_to_regions[key],
                                               path_to_blacklist=path_to_blacklist,
                                               #metrics=metadata_bc[key],
                                               #valid_bc=bc_passing_filters[key],
                                               n_cpu=max_cpu,
                                               project=key) for key in fragments_dict.keys()]

cistopic_obj = merge(cistopic_obj_list)
print(cistopic_obj)

# Merge Metadata with CisTopic Obj 
# Note: Creating a CisTopic object adds the 'sample' as a suffix - we need to make the metadata match
cell_data = pd.read_csv(os.path.join(work_dir, 'scATAC/full_metadata.csv'))
cell_data['cisTopic_barcode'] = cell_data['Cell_Barcode'] +'___'+ cell_data['sample']
print(cell_data['cisTopic_barcode'][0:5])
cell_data = cell_data.set_index('cisTopic_barcode')

# Use Only Cells with Metadata
cells_to_keep = np.array(cell_data.index)
cistopic_obj = cistopic_obj.subset(cells=cells_to_keep, copy=True)
cistopic_obj.add_cell_data(cell_data)
print(cistopic_obj)

# Remove doublets (if required)
#scrub = scr.Scrublet(cistopic_obj.fragment_matrix.T, expected_doublet_rate=0.1)
#doublet_scores, predicted_doublets = scrub.scrub_doublets()
#scrub.call_doublets(threshold=0.3)
#scrub.plot_histogram()
#scrublet = pd.DataFrame([scrub.doublet_scores_obs_, scrub.predicted_doublets_], columns=cistopic_obj.cell_names, index=['Doublet_scores_fragments', 'Predicted_doublets_fragments']).T
#singlets = cistopic_obj.cell_data[cistopic_obj.cell_data.Predicted_doublets_fragments == False].index.tolist()
#cistopic_obj = cistopic_obj.subset(singlets, copy=True)


# Save CisTopic Obj
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))


#################################################################################################################################
##################################################### Find Biomart Host #########################################################
#################################################################################################################################

# Find best biomart host
ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}

def test_ensembl_host(adata, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name','transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in adata.var_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(adata.var_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(adata, ensembl_version_dict[version], 'hsapiens')
    except:
        print('Host not reachable')

v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")
