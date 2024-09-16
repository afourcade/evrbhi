%% NVR epoching around R-peaks and exclusion noisy epochs
% 2017 (edited 2018) by Felix Klotzsche
% 2024 update by Antonin Fourcade

function NVR_07_epoch(cropstyle, mov_cond, path_data, epoch_begin, epoch_end)
%% 1.Set Variables
%clc
%clear all

% if needed, set voltage threshold for epoch exclusion
% th_exclusion = 100; %in uV

%1.1 Set different paths:
% input paths:
path_data_eeg = [path_data 'EEG/'];
path_in_eeg = [path_data_eeg '06_ICA/' mov_cond '/' cropstyle '/']; 

% output paths:

path_out_eeg = [path_data_eeg '07_epoch/' mov_cond '/' cropstyle '/']; 
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
    
     EEG = pop_epoch( EEG, {'ECG'}, [epoch_begin epoch_end], ...
    'epochinfo', 'yes');
    
    % If needed:
    % run identification of noisy epochs applying a mere threshold based
    % criterion: epochs with measurement values outside of the range[-100;
    % 100][uV] are rejected as "noisy"
    % Exclude ECG HEOG and VEOG from the rejection process
%     ignore_chans = [];
%     ignore = {'ECG', 'HEOG', 'VEOG'};
%     for igc = 1:numel(ignore)
%         idx = find(strcmp({EEG.chanlocs.labels}, ignore{igc}));
%         if ~isempty(idx)
%             ignore_chans(end+1) = idx;
%         end
%     end
% 
%     prep_chans = setdiff([1:33], ignore_chans);
%     EEG = pop_eegthresh(EEG,1,prep_chans ,-th_exclusion,th_exclusion,epoch_begin,epoch_end,0,1);
   
    EEG = eegh(com,EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);  
    EEG = pop_saveset( EEG, [filename  '_epoch.set'] , path_out_eeg);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end
end