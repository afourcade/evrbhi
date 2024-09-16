%% Plot cortical sources for NeVRo study
% Project the spatial patterns from SSD, CSP, and SPoC into source space 
% by using eLORETA.
%
% Launch this script from top level of the NeVRo repository.
% BEFORE RUNNING: 
% - Download the "New York Head" Leadfield from https://www.parralab.org/nyhead/
% - make sure the file is called "sa_nyhead.mat"
% - Place it in "./Analysis/SourceReconstruction/eLORETA_MJ/"
% - Please note that since that file is very large (>600MB) it is on the
% .gitignore list and will not be committed/downloaded if you use this 
% within our repository or a forked version. 

% CREDITS:
% The centerpiece of this script are the calls to `m_eLORETA_nvr.m` which
% is a slightly adapted version of and a wrapper to code written by 
% MINA JAMSHIDI IDAJI
% GUIDO NOLTE
% STEFAN HAUFE
% 
% For the original work, please refer to:
% http://bbci.de/supplementary/EEGconnectivity/BBCB.html
% https://github.com/minajamshidi 
% and: 
% METH toolbox 
% https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html 
% (as of 06/11/2019)
%
% If you make use of code in `SourceReconstruction`, please cite the 
% following paper:
% Haufe, S., & Ewald, A. (2016):
% A simulation framework for benchmarking EEG-based brain connectivity 
% estimation methodologies. 
% Brain topography, 1-18
% 
% 
% 2020 -- Felix Klotzsche -- 

%%
% set path to source plotting tools provided by MJ Idaji: 

path_MJ_tb = fullfile('D:','NeVRo', 'Analysis', 'SourceReconstruction');
addpath(genpath(path_MJ_tb))
% addpath(genpath(fullfile('.', 'Analysis', 'Statistics', 'Utils')));

SJs_nomov = [2 4 5 6 7 8 9 11 13 14 15 17 18 20 21 22 24 25 26 27 28 29 30 31 34 36 37 39 44];
SJs_mov = [2 4 5 6 8 9 11 13 14 17 18 20 21 22 24 25 27 28 29 31 34 36 37 39 44];

dirPatterns = 'D:\NeVRo\Data\Patterns\';
dirOut = [dirPatterns 'Sources\'];
data_path = 'D:\NeVRo\Data\12_eeglab2fieldrip\';
if ~exist(dirOut, 'dir')
    mkdir(dirOut)
    disp(['Created folder: ', dirOut])
end


% loading the model is slow, so I just do it once.

found_headmodel = exist(fullfile(path_MJ_tb, 'eLORETA_MJ','sa_nyhead.mat'), ...
                        'file') == 2;
assert(found_headmodel, ... 
       ['Looks like you have not loaded the New York Head model with the ', ...
       'name "sa_nyhead.mat" in the correct place. Check the header of ', ...
       '"plot_sources.m" for instructions.']);
mymodelfile= 'sa_nyhead';
load(mymodelfile);

m_conds = {'avg'}; %{'mov', 'nomov'} 

% Specify a single subject or empty vec for all subjects:
sub_id = SJs_nomov; %'NVR_S35'; %'NVR_S08'; %[];

patterns = {'HAvsLA_mean_t_0.328_0.372', 'HA_mean_t_0.328_0.372', 'LA_mean_t_0.328_0.372'};
% patterns = {'HAvsLA_mean_epoch'}; % {'HAvsLA_t_0.348'} %'SSD_1', 'SSD_2', 'SSD_3', 'SSD_4', ...
%patterns_time = 0.344;
       
save_plot = true;
save_format = 'svg'; %'epsc'
show_plot = false;
save_mats = true;
load_mats = true;

%% Source reconstruction

for mc=1:numel(m_conds) 
    file = fullfile(data_path, [m_conds{mc} '_mov_nomov'], 'SBA', [m_conds{mc} '_HEP_HA_LA.mat']);
    load(file);
    Pattern_mats = struct();
    for i=1:length(sub_id)
        subjHA = all_subjHA{SJs_nomov == sub_id(i)};
        subjLA = all_subjLA{SJs_nomov == sub_id(i)};
        subject_name = ['SJ' num2str(sub_id(i),'%02.f')];
        
        fprintf('---- Running subject  %s  -- condition:  %s  ----\n', ...
            subject_name, m_conds{mc});
       
        subjHAvsLA = subjHA.avg - subjLA.avg;
        
        % Note that channels TP9 and TP10 were not recorded in the NY head
        % model.
        replace_tpchans = false;
        tpchans = {'TP9', 'TP10'};
        replacements = {'TP7', 'TP8'};
        % drop channels not in the NY standard head model:
        if replace_tpchans
%             for tpchan_idx = 1:length(tpchans)
%                 tpchan = tpchans(tpchan_idx);
%                 idx = ismember(patt_tab.chanlocs, tpchan);
%                 patt_tab.chanlocs{idx} = replacements{tpchan_idx};
%                 bads = {'HEOG', 'VEOG'};
%                 idx_bads = ismember(subjHA.label, bads);
%             end
        else
            bads = {'HEOG', 'VEOG', 'TP9', 'TP10', 'ECG'};
            idx_bads = ismember(subjHA.label, bads);
        end
        
        subjHAvsLA = subjHAvsLA(~idx_bads, :);
        chan_labels = subjHA.label(~idx_bads);
        subjHA = subjHA.avg(~idx_bads, :);
        subjLA = subjLA.avg(~idx_bads, :);
        
        
        P_mats = struct();
        
        for j = 1:length(patterns)
            fprintf('Plotting: %s\n', patterns{j});
            %Look at pattern at t = pattern_time 
            %0.344 =  subjHA.time(162)
            %0.348 = subjHA.time(163)
            %pattern_time = subjHA.time(163); %double(0.344000000000000); Idk why this does not work....
%             idx_time = find(subjHA.time == pattern_time);
%             subjHAvsLA_pattern = subjHAvsLA(:, idx_time);
            switch patterns{j}
                case 'HAvsLA_mean_t_0.328_0.372'
                    pattern = mean(subjHAvsLA(:, 158:169),2); %mean between t = 0.328 and 0.372
      %             subjHAvsLA_pattern = mean(subjHAvsLA,2); %mean during whole epoch
                case 'HA_mean_t_0.328_0.372'
                    pattern = mean(subjHA(:, 158:169),2);
                case 'LA_mean_t_0.328_0.372'
                    pattern = mean(subjLA(:, 158:169),2);
            end                           

            P_mat = m_eLORETA_nvr(sa, pattern, ...
                chan_labels, patterns{j}, subject_name, m_conds{mc}, ...
                save_plot, save_format, dirOut, show_plot, 'viridis');
            if (~show_plot)
                close all;
            end
            P_mats(j).pattern = patterns{j};
            P_mats(j).weights = P_mat;
        end
        Pattern_mats(i).subject = subject_name;
        Pattern_mats(i).pattern_weights = P_mats;
    end
    
    switch m_conds{mc}
        case 'avg'
            P_avg = Pattern_mats;
            if (save_mats)
                dirResults = fullfile(dirPatterns, 'Sources', 'avg');
                if ~exist(dirResults, 'dir') 
                    mkdir(dirResults);
                    fprintf('Created dir: %s.\n', dirResults);
                end
                filename = fullfile(dirResults, 'source_patterns_avg.mat');
                save(filename, 'P_avg');
            end
        case 'mov'
            P_mov = Pattern_mats;
            if (save_mats)
                dirResults = fullfile(dirPatterns, 'Sources', 'mov');
                if ~exist(dirResults, 'dir') 
                    mkdir(dirResults);
                    fprintf('Created dir: %s.\n', dirResults);
                end
                filename = fullfile(dirResults, 'source_patterns_mov.mat');
                save(filename, 'P_mov');
            end
        case 'nomov'
            P_nomov = Pattern_mats;
            if (save_mats)
                dirResults = fullfile(dirPatterns, 'Sources', 'nomov');
                if ~exist(dirResults, 'dir') 
                    mkdir(dirResults);
                    fprintf('Created dir: %s.\n', dirResults);
                end
                filename = fullfile(dirResults, 'source_patterns_nomov.mat');
                save(filename, 'P_nomov');
            end
    end             
end

%% Plot averages:

% get color maps:
load cm17;
cm_vir = viridis();

patterns = {'HAvsLA_mean_t_0.328_0.372', 'HA_mean_t_0.328_0.372', 'LA_mean_t_0.328_0.372'};

%patterns = {'SSD_1', 'SSD_3'}; %{'SPOC', 'CSP_max', 'CSP_min'};

m_conds = {'avg'}; %{'mov', 'nomov'}

for mc=1:numel(m_conds) 
    switch m_conds{mc}
        case 'avg'
           if load_mats 
                fprintf('Loading data from disk.\n');
                if exist('P_avg', 'var')
                    fprintf('Overwriting data in memory w/ data from disk.\n');
                end
                dirResults = fullfile(dirPatterns, 'Sources', 'avg');
                filename = fullfile(dirResults, 'source_patterns_avg.mat');
                load(filename);
            end
            Pattern_mats = P_avg; 
        case 'mov'
            if load_mats 
                fprintf('Loading data from disk.\n');
                if exist('P_mov', 'var')
                    fprintf('Overwriting data in memory w/ data from disk.\n');
                end
                dirResults = fullfile(dirPatterns, 'Sources', 'mov');
                filename = fullfile(dirResults, 'source_patterns_mov.mat');
                load(filename);
            end
            Pattern_mats = P_mov;
                
        case 'nomov'
            if load_mats 
                fprintf('Loading data from disk.\n');
                if exist('P_nomov', 'var')
                    fprintf('Overwriting data in memory w/ data from disk.\n');
                end
                dirResults = fullfile(dirPatterns, 'Sources', 'nomov');
                filename = fullfile(dirResults, 'source_patterns_nomov.mat');
                load(filename);
            end
            Pattern_mats = P_nomov;
    end
            
    for p=1:length(patterns)
        p_idx = ismember({Pattern_mats(1).pattern_weights.pattern}, ...
            patterns{p});
        all_sub_pats = zeros(2004, size(Pattern_mats, 2));
        for sub=1:size(Pattern_mats, 2)
            m = Pattern_mats(sub).pattern_weights(p_idx).weights;
            all_sub_pats(:, sub) = m / norm(m);
        end

        mean_abs_patt = mean(abs(all_sub_pats), 2); 
        P = mean_abs_patt;
        colormap = cm_vir;
        scale_unit = 'a.u.';
        smooth = 1;
        views = 1:8;
        save_plot = true;
        savepath = fullfile(dirOut, m_conds{mc}, patterns{p}, 'mean');
        if (save_plot)
            if ~exist(savepath, 'dir'); mkdir(savepath); end
        end
        savename = fullfile(savepath, 'source');
        allplots_cortex_mina(sa, P(sa.cortex2K.in_to_cortex75K_geod), ... 
                    [min(P) max(P)], colormap, scale_unit, smooth, 'views', views, ...
                    'save', save_plot, 'savename', savename, ...
                    'saveformat', save_format);
        if (~show_plot)
            close all;
        end
        % Find mni coordinates for max in P - added by Antonin
        P_75k = P(sa.cortex2K.in_to_cortex75K_geod);
        [max_P_75k, idx_P_75k_max] = max(P_75k);
        mni_smooth = sa.cortex75K.vc_smooth(idx_P_75k_max,:);
        % Find mni coordinates for 3 max in P - added by Antonin
        [a, ix]=sort(P_75k);
        idx_max_P_75k_3= ix(end-2:end);
        mni_smooth_3 = sa.cortex75K.vc_smooth(idx_max_P_75k_3,:);
        
    end
end