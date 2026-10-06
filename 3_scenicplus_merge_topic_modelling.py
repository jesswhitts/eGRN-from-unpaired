#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
import pickle
from pycisTopic.lda_models import evaluate_models

work_dir = './'

# Combine Models
models = []
for (root, dirs, file) in os.walk(work_dir):
    for f in file:
        filename = os.path.join(root, f)
        if ('Topic' in filename) & ('.pkl' in filename):
            model = pickle.load(open(filename, "rb"))
            models.append(model)

# Save Models
pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models.pkl'), 'wb'))

# Evaluate Models
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

numTopics = 100
model = evaluate_models(models,
                       select_model = numTopics,
                       return_model = True,
                       metrics = ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics = False,
                       save=os.path.join(work_dir,'scATAC/models.pdf'))
