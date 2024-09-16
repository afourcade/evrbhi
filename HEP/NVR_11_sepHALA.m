%% NVR Select HA and LA epochs
% Select epochs with two or 1 ECG events and two (identical) or l arousal
% events (LA or HA)
% 2024 by Antonin Fourcade

function NVR_11_sepHALA(cropstyle, mov_cond, path_data)
%% 1.Set Variables
%clc
%clear all

%1.1 Set different paths:
% input paths:
path_data_eeg = [path_data 'EEG/'];
path_in_eeg = [path_data_eeg '10_rmbase/' mov_cond '/' cropstyle '/']; 

% output paths:
path_out_eeg = [path_data_eeg '11_sepHALA/' mov_cond '/' cropstyle '/'];
if ~exist(path_out_eeg, 'dir'); mkdir(path_out_eeg); end
path_out_eeg_LA = [path_out_eeg 'LA/'];
if ~exist(path_out_eeg_LA, 'dir'); mkdir(path_out_eeg_LA); end
path_out_eeg_HA = [path_out_eeg 'HA/'];
if ~exist(path_out_eeg_HA, 'dir'); mkdir(path_out_eeg_HA); end

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
    
    for i = 1:length(EEG.epoch)
        nb_events = length(EEG.epoch(i).event);
        idx_ECG = logical(count(EEG.epoch(i).eventtype,'ECG'));
        idx_aro = logical(ones(1,nb_events) - idx_ECG);
        nb_events_ECG = sum(idx_ECG);
        nb_events_aro = sum(idx_aro);
        
        if nb_events_ECG < 3
            aro_types = EEG.epoch(i).eventtype(logical(idx_aro));
            switch nb_events_aro
                case 2 % 2 arousal events in one epoch                   
                    if aro_types{1} == aro_types{2} % 2 arousal events are the same -> good
                        if aro_types{1} == '1'
                            % Classify as LA
                            eventLA = EEG.epoch(i).event(idx_aro);
                            EEG.event(eventLA(1)).type = 'LA';
                            EEG.event(eventLA(2)).type = 'LA';
                        elseif aro_types{1} == '3'
                            % Classify as HA                      
                            eventHA = EEG.epoch(i).event(idx_aro);
                            EEG.event(eventHA(1)).type = 'HA';
                            EEG.event(eventHA(2)).type = 'HA';
                        end
                    end
                case 1
                    if aro_types{1} == '1'
                            % Classify as LA
                            eventLA = EEG.epoch(i).event(idx_aro);
                            EEG.event(eventLA).type = 'LA';
                        elseif aro_types{1} == '3'
                            % Classify as HA
                            eventHA = EEG.epoch(i).event(idx_aro);
                            EEG.event(eventHA).type = 'HA';
                    end
            end
        end
    end
    
    EEG_LA = pop_selectevent(EEG, 'type', 'LA');
    EEG_HA = pop_selectevent(EEG, 'type', 'HA');

    %Save Files
    EEG = pop_saveset(EEG_LA, [filename  '_LA.set'] , path_out_eeg_LA);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    EEG = pop_saveset(EEG_HA, [filename  '_HA.set'] , path_out_eeg_HA);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end
end