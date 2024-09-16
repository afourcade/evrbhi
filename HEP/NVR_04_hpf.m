%% NVR High-Pass Filtering
% 2017 (edited 2018) by Felix Klotzsche
% 2024 update by Antonin Fourcade

function NVR_04_hpf(cropstyle, m_cond, path_data, locutoff)

%% 1.Set Variables
%clc
%clear all
mov_cond = m_cond;
%1.1 Set different paths:
% input paths:
path_data_eeg = [path_data 'EEG/'];
path_in_eeg = [path_data_eeg '03_eventsAro/' mov_cond '/' cropstyle '/']; 

% output paths:

path_out_eeg = [path_data_eeg '04_HPF/' mov_cond '/' cropstyle '/']; 
if ~exist(path_out_eeg, 'dir'); mkdir(path_out_eeg); end

%1.2 Get data files
files_eeg = dir([path_in_eeg '*.set']);
files_eeg = {files_eeg.name};

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%%

for isub = 1:length(files_eeg)
    
    %1.5 Set filename:
    filename = files_eeg{isub};
    filename = strsplit(filename, '.');
    filename = filename{1};
    
    
    %% 2.Import EEG data
    [EEG, com] = pop_loadset([path_in_eeg, filename, '.set']);
    EEG = eegh(com,EEG);
    EEG.setname=filename;
    
    %% 3.High-pass filtering

    % Hamming windowed sinc FIR filter
    EEG = pop_eegfiltnew(EEG,locutoff);


    %% 4.Save file
    EEG = eegh(com,EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);  
    EEG = pop_saveset( EEG, [filename  '_hpf' num2str(locutoff) '.set'] , path_out_eeg);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end
