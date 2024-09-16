%% Convert EEGLAB files into Fieldtrip structure
% For Cluster-based permutation testing
% 2024 by Antonin Fourcade

function NVR_12_eeglab2fieldrip(cropstyle, mov_cond, path_data)
%% 1.Set Variables
%clc
%clear all

%1.1 Set different paths:
% input paths:
path_data_eeg = [path_data 'EEG/'];
path_in_eeg_LA = [path_data_eeg '11_sepHALA/' mov_cond '/' cropstyle '/LA/'];
path_in_eeg_HA = [path_data_eeg '11_sepHALA/' mov_cond '/' cropstyle '/HA/']; 

% output paths:
path_out_eeg = [path_data_eeg '12_eeglab2fieldrip/' mov_cond '/' cropstyle '/'];
if ~exist(path_out_eeg, 'dir'); mkdir(path_out_eeg); end

%1.2 Get data files
files_eeg_LA = dir([path_in_eeg_LA '*.set']);
files_eeg_LA = {files_eeg_LA.name};
files_eeg_HA = dir([path_in_eeg_HA '*.set']);
files_eeg_HA = {files_eeg_HA.name};

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%%
all_subjLA = cell(1,length(files_eeg_LA));
all_subjHA = cell(1,length(files_eeg_HA));
if length(files_eeg_LA) ~= length(files_eeg_HA)
    fprintf(['ERROR: Not same number of SJs for HA and LA']);
    return
end

for isub = 1:length(files_eeg_LA)
    
    %1.5 Set filename:
    filename_LA = files_eeg_LA{isub};
    filename_LA = strsplit(filename_LA, '.');
    filename_LA = filename_LA{1};
    
    filename_HA = files_eeg_HA{isub};
    filename_HA = strsplit(filename_HA, '.');
    filename_HA = filename_HA{1};
    
    
    %% 2.Import EEG data
    [EEG_LA, com] = pop_loadset([path_in_eeg_LA, filename_LA, '.set']);
    EEG_LA = eegh(com,EEG_LA);
    EEG_LA.setname=filename_LA;
        
    [EEG_HA, com] = pop_loadset([path_in_eeg_HA, filename_HA, '.set']);
    EEG_HA = eegh(com,EEG_HA);
    EEG_HA.setname=filename_HA;

    all_subjLA{1,isub} = eeglab2fieldtrip(EEG_LA, 'timelockanalysis','chan_loc');
    all_subjLA{1,isub}.fsample=EEG_LA.srate;
    all_subjHA{1,isub} = eeglab2fieldtrip(EEG_HA, 'timelockanalysis','chan_loc');
    all_subjHA{1,isub}.fsample=EEG_HA.srate;
    
end

save([path_out_eeg mov_cond '_HEP_HA_LA.mat'],'all_subjLA','all_subjHA')
