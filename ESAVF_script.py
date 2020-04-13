# -*- coding: utf-8 -*-
"""
Script for automatic ESA of verbal fluency responses from Danish speakers
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
from os import path
from os.path import join
from sys import argv
import pickle
import dill
import pandas as pd
import re

#Call structure
#python [script.name.py] [particiant id number, e.g. 0001] [home_path]

#==============================================================================
# VARIABLES
#%%============================================================================
#Get variables from command line ipnuts
if __name__ == '__main__':  # only looks for command line inputs if script is called from command line
    ID          = argv[1]
    home_path   = argv[2]

#Define misc variables
#ID = "0001"
textfile = "id{0:s}.txt".format(str(ID))
#home_path = "/Users/au183362/Documents/postdoc/Parkinson-DBS/Eira_Aksnes/"
file_path = join(home_path, 'raw', textfile)
save_path = join(home_path, 'output')
sem_path = join(home_path, 'misc')

#iterations = 30000 # input to Semantic-module

#From Yund et al.:
#ESA switch threshold is "PercentRESA * MeanAllESA" or the (Minimum sequence ESA)+0.00001
#so that there is at least one ESA switch.
PercentRESA = 0.75


#==============================================================================
# NB! Run the following cell only once! It downloads the Danish wikipedia...
#%%============================================================================
#!python -m dasem.wikipedia download --verbose


#==============================================================================
# NB! The following code take 6-12 mins to run, so be patient...
#%%============================================================================
t = time.time()
#On AH's macbook: Semantic(30000) = 500-535 sec ≈ 8-9 min
#semantic = Semantic(iterations)
with open(join(sem_path, 'semantic.pkl'), 'rb') as f:  # Python 3: open(..., 'rb')
    semantic = dill.load(f)
elapsed = time.time()-t
print(elapsed)


#==============================================================================
# DATA IMPORT + ESA
#%%============================================================================
inputs = [line.rstrip() for line in open(file_path)]
for i, inp in enumerate(inputs):
    print(u"{0:d}: {1:s}".format(i, inp))

ESA_table = semantic.relatedness(inputs).round(3)
#print(raw_ESA_table[:,0])
##### NB! unicode not defined?!
# for i in range(int(math.ceil((len(inputs)/12.)))):
#     print("  ".join(map(unicode,inputs[12*i:12*(i+1)])))
#     print(ESA_table[:,12*i:12*(i+1)])
#     print("\n")


#==============================================================================
# EXTRACTING SUMMARY VALUES
#%%============================================================================
#Total number of words (animals), correct, duplicates, and unknowns
totals = {'N_Totl' : len(inputs),
    'N_Uniq' : len(np.unique(inputs)),
    'N_Dups' : len([dupl for dupl in collections.Counter(inputs).values() if dupl>1]),
    'N_UNKN' : sum([np.isnan(unk).all() for unk in ESA_table])}

print("Participant ID: {0:s}".format(str(ID)))
print("\nTotal number of \nwords: {0:d}\ncorrects: {1:d}"
      "\nduplicates: {2:d}\nunknowns: {3:d}\n".format(totals['N_Totl'],
      totals['N_Uniq'], totals['N_Dups'], totals['N_UNKN']))

#Extract ESA values for continuous (i.e., chronological) pairs
ESA = []
DictESA = {}
SumESA = 0.0

for n, i in enumerate(inputs):
    if n > 0:
        print(u"{0:s}-{1:s} = {2:2.3f}".format(inputs[n-1], i, ESA_table[n, n-1]))
        OneESA = ESA_table[n, n-1]
        DictESA[u"{0:s}-{1:s}".format(inputs[n-1], i)] = OneESA
        ESA = ESA + [OneESA]
        #ESA.append(OneESA)
        SumESA = SumESA + OneESA

# creating an "error-file" in case there are nans in the data
if np.isnan(np.sum(ESA_table)):
    #create/open file for writing error-file
    outfile = open(join(save_path, 'errors', str(ID) +'_nan.txt'), 'w')
    for n, i in enumerate(inputs):
        if n > 0:
            outfile.write(u"{0:s}-{1:s} = {2:2.3f}\n".format(inputs[n-1], i, ESA_table[n, n-1]))

MeanESA = SumESA / len(ESA)

print("\nMean (chronological pairwise) ESA for Participant ID {0:s}: {1:2.3f}".format(str(ID), MeanESA.round(3)))

#Extract ESA values for all possible pairs (still mainly clumsy Yund-code, but it works)
DictAllESA = {}
N_AllESA = 0
SumAllESA = 0.0

for i in range(len(inputs)-1):
    for j in range(i+1, len(inputs)):
            OneESA = ESA_table[i,j]
            DictAllESA[u"{0:s}-{1:s}".format(inputs[i], inputs[j])] = OneESA
            N_AllESA = N_AllESA + 1
            SumAllESA = SumAllESA + OneESA

if (N_AllESA == 0):
    MeanAllESA = 0.0
    SOI = 0.0
else:
    MeanAllESA = SumAllESA / N_AllESA
    if MeanAllESA == 0.0:
        SOI = 0.0
    else:
        SOI = MeanESA / MeanAllESA

#print(inputs)
#print(ESA_table)
#print(DictAllESA)
ESAs = {'MeanESA' : MeanESA,
    'MeanAllESA' : MeanAllESA,
    'SOI' : SOI}
print("\nMean All ESA for Participant ID {0:s}: {1:2.3f}".format(str(ID), MeanAllESA.round(3)))
print("\nSOI Participant ID {0:s}: {1:2.3f}".format(str(ID), SOI.round(3)))


#==============================================================================
# CLUSTERS & SWITCHES
#%%============================================================================
# Setting the ESA_Threshold
ESAsort = sorted(ESA) # creating a sorted copy of ESA
ESA_Threshold = PercentRESA * MeanAllESA  # ESA switch if ESA < ESA_Threshold
if ESA_Threshold < ESAsort[0]:
    ESA_Threshold = ESAsort[0] + 0.000001
print(ESA_Threshold)

#ESA switches
ESAsw = [int(esa<ESA_Threshold) for esa in ESA]

ESAswCount = collections.Counter(ESAsw).most_common() # counting the ones and the zeros of the switch-vector
#CSidx = re.finditer('(1, 0)', str(ESAsw))
CSidx = re.findall('(1, 0)', str(ESAsw))
#print(CSidx)
#print('%02d-%02d: %s' % (CSidx.start(), CSidx.end(), CSidx.group(0)))
NumCS = len(CSidx)

sumESAsw = [x[1] for x in ESAswCount if x[0]==1][0]
print("Number of switches: {0:d}\nNumber of clusters: {1:d}\n".format(sumESAsw, NumCS))
ESAs['sumESAsw'] = sumESAsw
ESAs['NumCS'] = NumCS

print("Switch index:")
for sw, i in zip(ESAsw, inputs):
    print(u"{0:d}: {1:s}".format(sw, i))

SizeCS = [sum(1 for _ in g)+1 for _, g in groupby(ESAsw) if _==0] # we add one to each cluster-sum for 1st word in the cluster
SumCS = sum(SizeCS)
print("Cluster sizes: {}".format("".join(str(SizeCS))))
print("Accumulated cluster size: {}".format(SumCS))
ESAs['SumCS'] = SumCS

ESAcs = ESAsw[:]
print("All 'clusters': {}".format("".join(str(ESAcs))))
count = 0
#[ESAcs.insert(i, SizeCS[count]) for i, sw in enumerate(ESAsw)]
for i, sw in enumerate(ESAsw):
    if (sw == 0 and ESAsw[i-1] == 1):
        # print(str(i))
        # print(str(count))
        ESAcs[i] = SizeCS[count]
        # print("ESAcs & SiceCS: {0:d} & {1:d}".format(ESAcs[i], SizeCS[count]))
        # if SizeCS[count] > 2:
        #     ESAcs[i+1:i+SizeCS[count]-1] = np.tile(0, SizeCS[count]-1)
        count = count + 1
print("All 'clusters': {}".format("".join(str(ESAcs))))

#==============================================================================
# Storing all printable variables in a pandas DataFrame
#%%============================================================================
total_df = pd.DataFrame(totals, index=[ID])
ESA_df = pd.DataFrame(ESAs, index=[ID])
dfs = [total_df, ESA_df]
new_df = pd.concat(dfs, axis=1)
print(new_df)

#==============================================================================
# CREATE & SAVE OUTPUT
#%%============================================================================

#Create output file for each subject and save vairables
# Saving the objects:
with open(join(save_path, 'ESAtable_{0:s}.pkl'.format(str(ID))), 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([ESA_table, DictAllESA, totals, ESAs], f)

# # How to load the objects:
# with open(join(save_path, 'esavf{0:s}.pkl'.format(str(ID))), 'rb') as f:  # Python 3: open(..., 'rb')
#     ESA_table, DictAllESA, totals, ESAs = pickle.load(f)

# (Load existing csv) and save output as csv for handling in R
try:
    if path.exists(join(save_path, 'ESA_sumvals.csv')):
        old_df = pd.read_csv(join(save_path, 'ESA_sumvals.csv'))
        if not(ID in old_df.index):
            new_old_df = [old_df, new_df]
            upd_df = pd.concat(new_old_df)
        if ID in old_df.index:
            old_df.loc[ID] = new_df.loc[ID]
            upd_df = old_df
        print(upd_df)
        upd_df.to_csv(join(save_path, 'ESA_sumvals.csv'))
    if not(path.exists(join(save_path, 'ESA_sumvals.csv'))):
        print(new_df)
        new_df.to_csv(join(save_path, 'ESA_sumvals.csv'))
except Exception:
    pass


#Create printed output file for each subject (.txt)
#NB! not all entries in the below code are relevant for our current analysis

#create/open file for writing output
outfile = open(join(save_path, str(ID) +'_out.txt'), 'w')

#summary values
outfile.write(u'Number of animals:   ')
outfile.write(format(totals['N_Totl'], u'>3d'))
outfile.write(u'\nNumber correct:      ')
outfile.write(format(totals['N_Uniq'], u'>3d'))
outfile.write(u'\nRepetitions:        ')
outfile.write(format(totals['N_Dups'], u'>3d'))
outfile.write(u'\nUnknown Animals:    ')
outfile.write(format(totals['N_UNKN'], u'>3d'))

#ESA analysis
outfile.write(u'\n\nSequence Mean ESA:  ')
outfile.write(format(MeanESA, u'10.4f'))
outfile.write(u'\nTotal Mean ESA:     ')
outfile.write(format(MeanAllESA, u'10.4f'))
outfile.write(u'\nSem Org Index (SOI):')
outfile.write(format(SOI, u'10.4f'))
outfile.write(u'\n\nESA switches:       ')
outfile.write(format(sum(ESAsw), u'>4d'))
outfile.write(u'  (')
outfile.write(format(totals['N_Totl']/sumESAsw, u'5.2f'))
outfile.write(u' words/switch)')
outfile.write(u'\n\nESA clustes:        ')
outfile.write(format(NumCS, u'>4d'))

outfile.write(u'  (')
outfile.write(format(SumCS/NumCS, u'5.2f'))
outfile.write(u' words/cluster)')
