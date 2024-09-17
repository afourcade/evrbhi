#############################################################################################################################################################################################
# The script plots the emotional arousal ratings of all NeVRo participants in the same plot and the mean of all participants
# The script also computes the descriptive statistics of the ratings and the standard deviation of the ratings for each participant
#
# 2024 by Antonin Fourcade
#############################################################################################################################################################################################
# %%
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Parameters
mov_cond = ['mov', 'nomov'] # head movement condition
c_style = 'SBA' # sequence o the rollercosaters: Space-Break-Andes
nb_samples = 270 # 1Hz ratings for SBA = 270s

# Data paths
data_path = 'E:/NeVRo/Data/'

#%%

#savepath
savepath = 'E:/NeVRo/Results'

# sj list
sjs = {
    "nomov": ['02', '04', '05', '06', '07', '08', '09', '11', '13', '14', '15', '17', '18', '20', '21', '22', '24', '25', '26', '27', '28', '29', '30', '31', '34', '36', '37', '39', '44'],
    "mov": ['02', '04', '05', '06', '08', '09', '11', '13', '14', '17', '18', '20', '21', '22', '24', '25', '27', '28', '29', '31', '34', '36', '37', '39', '44']   
      }

# rollercoster events
rc_events = {
            'rapid_curve_spiral':47, #space event
            'looping':60, #space event
            'spiral':73, #space event
            'fire':20+178, #andes event
            'steep_fall':24+178, #andes event
            'jump_1':31+178, #andes event
            'landing_1':33+178, #andes event
            'fire_looping':51+178, #andes event
            'jump_2':67+178, #andes event
            'landing_2':72+178 #andes event
              }

# array of samples
samples = np.arange(0,nb_samples)

# %%
# Loop over mov_cond
for mc in mov_cond:
    # non z-scored arousal ratings
    aro_rat_path = data_path + 'ratings/continuous/not_z_scored/' + mc + '/' + c_style + '/'
    files_aro_z = glob.glob(aro_rat_path + '*.txt')
    # keep only the files corresponding to the sjs list
    files_aro_z = [f for f in files_aro_z if any(sj in f for sj in sjs[mc])]

    # Title of the plot
    if mc == 'mov':
        prefix_mov = 'Free'
    else:
        prefix_mov = 'No'
    title = prefix_mov + ' head movement (' + mc + ')'
    # Plot all the participants in the same plot
    fig, ax = plt.subplots()
    aro_mean = np.zeros(nb_samples)
    for idx, f in enumerate(files_aro_z):
        sj = sjs[mc][idx]
        data = pd.read_csv(f, names=["idx", "ar", "latency"])["ar"]
        # rescale the data from [0, 100, :2] to [0, 50, :1]
        data = data/2
        # append the data to the mean
        aro_mean += data
        # plot the data
        ax.plot(samples, data, label=sj, alpha=0.8, linewidth=0.75)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Emotional arousal ratings (a.u.)')
        ax.set_title(title)
        # legend outside the plot
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    # plot the mean of all participants
    aro_mean /= len(files_aro_z)
    ax.plot(samples, aro_mean, label='Mean', color='black', linewidth=1.5)
    # set xlim, break lines and xticks at 270
    ax.set_xlim([0, nb_samples])
    ax.set_xticks(np.arange(0, nb_samples+1, 30))
    ax.axvline(x=148, color='grey', linestyle='dotted', linewidth=1)
    ax.axvline(x=178, color='grey', linestyle='dotted', linewidth=1)

    # %%
    # Save the plot
    #fig.savefig(savepath + 'aro_ratings_' + mc + '.png', dpi=600)
    fig.savefig(savepath + 'aro_ratings_' + mc + '.svg')

# %%
# Descriptive statistics
# create a dataframe with all the data
aro_all = {'mov':pd.DataFrame(), 'nomov':pd.DataFrame()}
# loop over mov_cond
for mc in ['mov', 'nomov']:
    # non z-scored arousal ratings
    aro_rat_path = data_path + 'ratings/continuous/not_z_scored/' + mc + '/' + c_style + '/'
    files_aro_z = glob.glob(aro_rat_path + '*.txt')
    # keep only the files corresponding to the sjs list
    files_aro_z = [f for f in files_aro_z if any(sj in f for sj in sjs[mc])]
    for idx, f in enumerate(files_aro_z):
        sj = sjs[mc][idx]
        data = pd.read_csv(f, names=["idx", "ar", "latency"])["ar"]
        aro_all[mc][sj] = data
        # rescale the data from [0, 100, :2] to [0, 50, :1]
        aro_all[mc][sj] = aro_all[mc][sj]/2
    # tranform to long format
    aro_all[mc] = aro_all[mc].melt(var_name='sj', value_name='ar')
    # descriptive statistics averaged over participants
    aro_all[mc].describe().transpose().to_csv(savepath + 'aro_ratings_' + mov_cond + '_descriptive.csv')

# create a new column 'mov_cond'
aro_all['mov']['mov_cond'] = 'mov'
aro_all['nomov']['mov_cond'] = 'nomov'
# concatenate the two dataframes
aro_all = pd.concat([aro_all['mov'], aro_all['nomov']])
# compute the std of the ratings for each participant
aro_std = aro_all.groupby('sj')['ar'].std()
# transform to dataframe with columns 'sj' and 'std'
aro_std = aro_std.reset_index()
# split the participants in two equal groups: low and high std
aro_std['group'] = pd.qcut(aro_std['std'], 2, labels=['low', 'high'])
# save the data
aro_std.to_csv(savepath + 'aro_ratings_std.csv', index=False)


# %%
