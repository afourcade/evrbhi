##########################################################################################
# Script to analyze the questionnaires data (Valence, Presence, STAI-T, UPPS)
# - descriptive statistics
# - test if valence and presence is different between mov_cond
# - reformat STAI-T and UPPS data
# - descriptive statistics of STAI-T and UPPS
#
# 2024 by Antonin Fourcade
##########################################################################################


#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

#%%
# Data paths
data_path_vp = 'E:/NeVRo/new_HEP_data_filtHP_0_3Hz/ratings/Valence_Presence/'
filename_vp = 'Valence_Presence_reformat.csv'
data_path_st = 'E:/NeVRo/new_HEP_data_filtHP_0_3Hz/ratings/State_Trait/'
filename_stait = 'STAI-T.csv'
filename_upps = 'UPPS_Sensation Seeking.csv'

# Load the data
data_vp = pd.read_csv(data_path_vp + filename_vp)
data_stait = pd.read_csv(data_path_st + filename_stait)
data_upps = pd.read_csv(data_path_st + filename_upps) 
# %% Valence and Presence
# presence (7-points Likert scale, from 1 to 7)
# valence (7-points Likert scale, from -3 to 3)
# descriptive statistics
data_vp.describe().transpose().to_csv(data_path_vp + 'Valence_Presence_descriptive.csv')
# %%
# descriptive statistics by mov_cond
data_vp.groupby('mov_cond').describe().transpose().to_csv(data_path_vp + 'Valence_Presence_descriptive_by_mov_cond.csv')	
# %%
# test if valence and presence is different between mov_cond
# separate data
data_vp_mov = data_vp[data_vp['mov_cond'] == 'mov']
data_vp_nomov = data_vp[data_vp['mov_cond'] == 'nomov']
# test
ttest_ind(data_vp_mov['valence'], data_vp_nomov['valence'])
ttest_ind(data_vp_mov['presence'], data_vp_nomov['presence'])

# %% State and Trait
# (1) the ‘Trait’ subscale of the ‘State-Trait Anxiety Inventory’ (STAI-T; Spielberger, 1983; Spielberger, 1989)
# (2) the ‘Sensation Seeking’ subscale of the ‘UPPS Impulsive Behaviour Scale’ (UPPS; Schmidt et al., 2008; Whiteside and Lynam, 2001) - 1 to 4

# get list of included sjs from data_vp
sj_list = data_vp['sj'].unique()

# reformat data_stait
# remove last row
data_stait = data_stait.iloc[:-1]
# reformat Subject Number column
data_stait['Subject Number'] = data_stait['Subject Number'].str.replace('NVR_S', '').astype(int)
# keep only included sjs
data_stait = data_stait[data_stait['Subject Number'].isin(sj_list)]
# create anxiety_trait_score column by summing the 20 items
data_stait['anxiety_trait_score'] = data_stait.iloc[:, 1:21].sum(axis=1)
# descriptive statistics of anxiety_trait_score
data_stait['anxiety_trait_score'].describe().to_csv(data_path_st + 'STAI-T_descriptive.csv')

# reformat data_upps
# reformat Subject Number column
data_upps['Subject Number'] = data_upps['Subject Number'].str.replace('S', '').astype(int)
# keep only included sjs
data_upps = data_upps[data_upps['Subject Number'].isin(sj_list)]
# all items are reversed scored
data_upps.iloc[:, 1:13] = 5 - data_upps.iloc[:, 1:13]
# create sensation_seeking_score column by averaging the 12 items
data_upps['sensation_seeking_score'] = data_upps.iloc[:, 1:13].mean(axis=1)
# descriptive statistics of sensation_seeking_score
data_upps['sensation_seeking_score'].describe().to_csv(data_path_st + 'UPPS_Sensation_Seeking_descriptive.csv')
