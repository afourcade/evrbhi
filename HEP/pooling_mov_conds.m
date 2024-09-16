%% Create a new dataset with all subjects averaged across mov-nomov
% If a SJ do not have a mov condition, take only nomov

% BEFOREHAND: Load mov and nomov datasets ('12_eeglab2fieldrip') and rename
% *_mov and *_nomov, respectively

% output paths:
path_data_eeg =  'D:/NeVRo/new_HEP_data_filtHP_0_3Hz/';
path_out_eeg = [path_data_eeg '12_eeglab2fieldrip_new/avg_mov_nomov/SBA/'];
if ~exist(path_out_eeg, 'dir'); mkdir(path_out_eeg); end

% All mov are included in nomov, but nomov have more SJs
SJs_nomov = [2 4 5 6 7 8 9 11 13 14 15 17 18 20 21 22 24 25 26 27 28 29 30 31 34 36 37 39 44];
SJs_mov = [2 4 5 6 8 9 11 13 14 17 18 20 21 22 24 25 27 28 29 31 34 36 37 39 44];

all_subjHA = cell(1,length(SJs_nomov));
all_subjLA = cell(1,length(SJs_nomov));

for i = 1:length(SJs_nomov)
    sj = find(SJs_mov==SJs_nomov(i)); %get index SJ_nomov in mov cond
    if isempty(sj) % If no mov cond, take only nomov cond
        all_subjHA{1,i} = all_subjHA_nomov{1,i};
        all_subjLA{1,i} = all_subjLA_nomov{1,i};
    else % compute average across mov_conds
        all_subjHA{1,i}.elec = all_subjHA_nomov{1,i}.elec;
        all_subjHA{1,i}.time = all_subjHA_nomov{1,i}.time;
        all_subjHA{1,i}.label = all_subjHA_nomov{1,i}.label;
        all_subjHA{1,i}.cfg = all_subjHA_nomov{1,i}.cfg;
        all_subjHA{1,i}.fsample = all_subjHA_nomov{1,i}.fsample;
        all_subjHA{1,i}.avg = (all_subjHA_nomov{1,i}.avg+all_subjHA_mov{1,sj}.avg)/2; %average mov-nomov
        all_subjHA{1,i}.var = (all_subjHA_nomov{1,i}.var+all_subjHA_mov{1,sj}.var)/2; %average mov-nomov
        all_subjLA{1,i}.elec = all_subjLA_nomov{1,i}.elec;
        all_subjLA{1,i}.time = all_subjLA_nomov{1,i}.time;
        all_subjLA{1,i}.label = all_subjLA_nomov{1,i}.label;
        all_subjLA{1,i}.cfg = all_subjLA_nomov{1,i}.cfg;
        all_subjLA{1,i}.fsample = all_subjLA_nomov{1,i}.fsample;
        all_subjLA{1,i}.avg = (all_subjLA_nomov{1,i}.avg+all_subjLA_mov{1,sj}.avg)/2; %average mov-nomov
        all_subjLA{1,i}.var = (all_subjLA_nomov{1,i}.var+all_subjLA_mov{1,sj}.var)/2; %average mov-nomov
    end
end

save([path_out_eeg 'avg_HEP_HA_LA.mat'],'all_subjLA','all_subjHA');