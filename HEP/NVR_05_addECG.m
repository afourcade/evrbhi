%% NVR add ECG to dataset
% 2024 by Antonin Fourcade
% Run HEPLAB beforehand to get R-peak events (ecg .mat files)!

function NVR_05_addECG(c_style, mov_cond, path_data)

% 1.1 Set different paths:
% input paths:
path_data_eeg = [path_data 'EEG/'];
path_in_eeg = [path_data_eeg '04_HPF/' mov_cond '/' c_style '/'];
path_data_ecg = [path_data 'ECG/'];
path_in_ecg = [path_data_ecg mov_cond '/' c_style '/'];

% output paths:
path_out_eeg = [path_data_eeg '05_add_ECG/' mov_cond '/' c_style '/'];
if ~exist(path_out_eeg, 'dir'); mkdir(path_out_eeg); end

%1.2 Get data files
files_ecg = dir([path_in_ecg '*.mat']);
files_ecg = {files_ecg.name};
files_eeg = dir([path_in_eeg '*.set']);
files_eeg = {files_eeg.name};

for isub=1:length(files_eeg)
filename_ecg = files_ecg{isub};

%1.3 Launch EEGLAB:
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% 1.4 Get subj name for getting the right event file later:
thissubject = files_eeg{isub};
thissubject = strsplit(thissubject, '.set');
thissubject = thissubject{1};    

% 1.5 Set filename:
filename_eeg = thissubject; 
filename_eeg = char(filename_eeg);
    
%% 2.Import EEG data
[EEG, com] = pop_loadset([path_in_eeg, filename_eeg '.set']);
EEG = eegh(com,EEG);

%% 3. Import ECG data
load([path_in_ecg, filename_ecg]); % Variable HEP

%% 4. Add ECG to EEG data
EEG.nbchan = EEG.nbchan +1;
EEG.data(EEG.nbchan,:) = single(HEP.ecg); %convert from double to single
EEG.chanlocs(EEG.nbchan).labels = 'ECG';

%% 5. Add ECG events to EEG events
k = length(EEG.event); % number of existing events in the EEG structure
urevents = num2cell(k+1:k+length(HEP.qrs));
evt = num2cell(HEP.qrs);
types = repmat({'ECG'},1, length(evt));

[EEG.event(1,k+1:k+length(HEP.qrs)).latency] = evt{:}; % assign latencies
[EEG.event(1,k+1:k+length(HEP.qrs)).type] = types{:}; % assign types
[EEG.event(1,k+1:k+length(HEP.qrs)).urevent] = urevents{:}; % assign event index

%% 6. Change arousal events (1,2,3) into string (otherwise EEGLAB not happy with 'ECG'...)
for i = 1:length(EEG.event)  
    if isa(EEG.event(i).type,'double')
        EEG.event(i).type = int2str(EEG.event(i).type); 
    end
end

%% 7.Save file
EEG = eegh(com,EEG);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);  
EEG = pop_saveset( EEG, [filename_eeg  '_ECG.set'] , path_out_eeg);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
clear k evt urevents types HEP
end