%% NVR Copy ICA weights from HP-1Hz dataset to dataset 
% 2024 by Antonin Fourcade
% Beforehand, need to run ICA decomposition on a copy of the data high-pass filtered at
% 1Hz!

function NVR_06_copy_ICA(cropstyle, mov_cond, path_data)
%% 1.Set Variables
%clc
%clear all

%1.1 Set different paths:
path_data_eeg = [path_data 'EEG/'];
path_data_eegICA = [path_data 'ICA_1Hz/'];
path_in_eegICA = [path_data_eegICA '05_cleanICA/' mov_cond '/' cropstyle '/'];
path_in_eegNEW = [path_data_eeg '05_add_ECG/' mov_cond '/' cropstyle '/'];

% output paths:
path_out_eeg_new = [path_data_eeg '06_ICA/' mov_cond '/' cropstyle '/'];
if ~exist(path_out_eeg_new, 'dir'); mkdir(path_out_eeg_new); end

%1.2 Get data files
files_eegICA = dir([path_in_eegICA '*.set']);
files_eegICA = {files_eegICA.name};
files_eegNEW = dir([path_in_eegNEW '*.set']);
files_eegNEW = {files_eegNEW.name};

for isub = 1:length(files_eegNEW)

%Launch EEGLAB:
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%Get subj name for getting the right event file later:
thissubjectNEW = files_eegNEW{isub};
thissubjectNEW = strsplit(thissubjectNEW, '.set');
thissubjectNEW = thissubjectNEW{1};

thissubjectICA = files_eegICA{isub};
thissubjectICA = strsplit(thissubjectICA, '.set');
thissubjectICA = thissubjectICA{1}; 

%Set filename:
filename_eegNEW = thissubjectNEW; 
filename_eegNEW = char(filename_eegNEW);

filename_eegICA = thissubjectICA; 
filename_eegICA = char(filename_eegICA);
    
%% Import EEG datasets
[EEGNEW, com] = pop_loadset([path_in_eegNEW, filename_eegNEW '.set']);
EEGNEW = eegh(com,EEGNEW);

[EEGICA, com] = pop_loadset([path_in_eegICA, filename_eegICA '.set']);
EEGICA = eegh(com,EEGICA);

% Add ICA weights to NEW dataset:
EEGNEW.icawinv = EEGICA.icawinv;
EEGNEW.icasphere = EEGICA.icasphere;
EEGNEW.icaweights = EEGICA.icaweights;
EEGNEW.icachansind = EEGICA.icachansind;
EEGNEW.icaact = EEGNEW.icaweights*EEGNEW.icasphere*EEGNEW.data(1:32,:); %don't include ECG (chan33)
EEGNEW.etc.icaweights_beforerms = EEGICA.etc.icaweights_beforerms;
EEGNEW.etc.icasphere_beforerms = EEGICA.etc.icasphere_beforerms;

% Save file
EEGNEW = eegh(com,EEGNEW);
EEGNEW.setname=filename_eegNEW;
EEGNEW = pop_saveset(EEGNEW, [filename_eegNEW  '_ICA.set'] , path_out_eeg_new);
end
end