#!/usr/bin/env python
# coding: utf-8


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import re
import pickle
import pandas as pd
import multiprocessing
max_cpu = multiprocessing.cpu_count()

from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import evaluate_models
from pycisTopic.clust_vis import run_umap
from pycisTopic.clust_vis import plot_metadata
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *

cwd = os.getcwd()
work_dir = './'
tmp_dir = '/scratch/wsspaces/jwhittle-tmp'



# Evaluate Models
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
models = pickle.load(open(os.path.join(work_dir, 'scATAC/models.pkl'), 'rb'))

numTopics = 150
model = evaluate_models(models,
                       select_model = numTopics,
                       return_model = True,
                       metrics = ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics = False,
                       save=os.path.join(work_dir, 'scATAC/models.pdf'))

cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))


run_umap(cistopic_obj, target = 'cell', scale = True)


plot_metadata(
    cistopic_obj,
    reduction_name = 'UMAP',
#    color_dictionary = {'sample': color_dict_line},
    variables = ['sample'],
    figsize = (10, 10),
    save=os.path.join(work_dir, 'scATAC/sample.pdf')
)


plot_metadata(
    cistopic_obj,
    reduction_name = 'UMAP',
#    color_dictionary = {'Expected_Driving_Aberration': color_dict_line},
    variables = ['Expected_Driving_Aberration'],
    figsize = (10, 10),
    save=os.path.join(work_dir, 'scATAC/translocation_partner.pdf')
)


plot_metadata(
    cistopic_obj,
    reduction_name = 'UMAP',
    #color_dictionary = {'Expected_Driving_Aberration': color_dict_line},
    variables = ['Classified_Celltype'],
    figsize = (10, 10),
    save=os.path.join(work_dir, 'scATAC/classified_celltype.pdf')
)


# Binarize Topics
# Otsu Method
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')

os.makedirs(os.path.join(work_dir, "scATAC", "region_sets", "topics_otsu"), exist_ok=True)
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome","Start","End"]
    ).to_csv(
        os.path.join(work_dir, "scATAC", "region_sets", "topics_otsu", f"{topic}.bed"),
        sep = "\t", header=False, index=False
    )
    
region_bin_topics_3k = binarize_topics(cistopic_obj, method='ntop', ntop=3_000)

os.makedirs(os.path.join(work_dir, "scATAC", "region_sets", "topics_3k"), exist_ok=True)
for topic in region_bin_topics_3k:
    region_names_to_coordinates(
        region_bin_topics_3k[topic].index
    ).sort_values(
        ["Chromosome","Start","End"]
    ).to_csv(
        os.path.join(work_dir, "scATAC", "region_sets", "topics_3k", f"{topic}.bed"),
        sep = "\t", header=False, index=False
    )


# Get Differentially Accessible Regions
print("Imputing Accessibility")
imputed_acc_obj = impute_accessibility(cistopic_obj, 
                                       selected_cells=None, 
                                       selected_regions=None, 
                                       scale_factor=10**6
                                      )

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, 
                                              scale_factor=10**4
                                             )

variable_regions = find_highly_variable_features(normalized_imputed_acc_obj,
                                                 min_disp = 0.05,
                                                 min_mean = 0.0125,
                                                 max_mean = 3,
                                                 max_disp = np.inf,
                                                 n_bins = 20,
                                                 n_top_features = None,
                                                 plot = True,
                                                 save = 'variable_features.pdf'
                                                )


print('Calculating DARs for each Celltype...')
markers_dict_celltype = find_diff_features(cistopic_obj, 
                                           imputed_acc_obj, 
                                           variable = 'Classified_Celltype', 
                                           var_features = variable_regions,
                                           contrasts = None,
                                           adjpval_thr = 0.05,
                                           log2fc_thr = np.log2(1.5),
                                           n_cpu = max_cpu,
                                           _temp_dir=tmp_dir,
                                           split_pattern = '-'
                                          )

os.makedirs(os.path.join(work_dir, "scATAC", "region_sets", "DARs_celltype"), exist_ok=True)
for celltype in markers_dict_celltype:
    region_names_to_coordinates(
        markers_dict_celltype[celltype].index
    ).sort_values(
        ["Chromosome","Start","End"]
    ).to_csv(
        os.path.join(work_dir, "scATAC", "region_sets", "DARs_celltype", f"{celltype.replace('/','')}.bed"),
        sep = "\t", header=False, index=False
    )
