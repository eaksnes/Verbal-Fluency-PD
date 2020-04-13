# -*- coding: utf-8 -*-
"""
Script for initializing and saving Semantic-object based on wikipedia-download
@author: Eira Aksnes & Andreas Højlund
@email:  & hojlund@cfin.au.dk
@github: https://github.com/ahoejlund/Verbal-Fluency-PD
"""

#==============================================================================
# IMPORTS
#%%============================================================================

from dasem.semantic import Semantic
import numpy as np
from io import open
import collections
from itertools import compress, groupby
import math
import time
from os.path import join
from sys import argv
#import pickle
import dill

#==============================================================================
# IMPORTS
#%%============================================================================
home_path = "/Users/au183362/Documents/postdoc/Parkinson-DBS/Eira_Aksnes/"
save_path = join(home_path, 'misc')

#==============================================================================
# NB! Run the following cell only once! It downloads the Danish wikipedia...
#%%============================================================================
#python -m dasem.wikipedia download --verbose

#wikipedia-dump should be saved to somewhere like
#~/dasem_data/wikipedia/dawiki-latest-pages-articles.xml.bz2


#==============================================================================
# Initialize Semantic-object and pickle-dump it
#%%============================================================================
iterations = 30000 # input to Semantic-module

t = time.time()
semantic = Semantic(iterations)
elapsed = time.time()-t
print(elapsed)

t = time.time()
with open(join(save_path, 'semantic.pkl'), 'wb') as f:  # Python 3: open(..., 'rb')
    #pickle.dump(semantic, f)
    dill.dump(semantic, f)
elapsed = time.time()-t
print(elapsed)
