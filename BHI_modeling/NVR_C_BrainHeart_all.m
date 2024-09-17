%% NeVRo BHI modeling - C_Brain-HRV computation within ROI
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to compute Brain-heart couplings, using the model
% developped by Catrambone et al. (2019).
% All the couplings are averaged with an EEG ROI.
% All timepoints corresponding to low (LA) and high (HA) emotional arousal
% are saved into a .csv file, for linear mixed modeling (LMM) analysis in
% R.
%
% 2024 Antonin Fourcade
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

% add analysis path
NVR_path = 'E:\NeVRo\Analysis';
addpath(NVR_path);

% set parameters
m_conds = {'nomov' 'mov'}; % head movement conditions
a_conds = {'LA' 'HA'};  % emotional arousal conditions
c_style = 'SBA'; % sequence of rollercoasters, Space-Break-Andes
%pow_style = 'full'; %'full'; 'flat'
%mtd_hrv = 'cwt'; %'cwt'; 'spwvd'
%mtd_eegpow = 'cwt'; %'cwt'; 'stft'

% set paths
data_path = 'E:\NeVRo\Data\';
fbp_path = [data_path 'Frequency_Bands_Power\'];
savepath_hf = [fbp_path 'BrainHF\'];
savepath_lf = [fbp_path 'BrainLF\'];
savename_hf = [savepath_hf 'C_BrainHF.mat']; 
savename_lf = [savepath_lf 'C_BrainLF.mat']; 
savepath_hf_roi = [fbp_path 'ROI_BrainHF\'];
savepath_lf_roi = [fbp_path 'ROI_BrainLF\'];
savename_hf_roi = [savepath_hf_roi 'C_BrainHF.mat']; 
savename_lf_roi = [savepath_lf_roi 'C_BrainLF.mat']; 
csv_filename = [fbp_path 'nvr_roi_brainheart.csv'];

% Define ROI
ROI = {'Pz'; 'P3'; 'P4'; 'P7'; 'P8'; 'O1'; 'O2'; 'Oz'};

% Initialize BHI variables
C_BrainHF = cell(length(m_conds),1);
C_BrainLF = cell(length(m_conds),1);
C_ROI_BrainHF = cell(length(m_conds),1);
C_ROI_BrainLF = cell(length(m_conds),1);

% All mov are included in nomov, but nomov have more SJs
SJs_nomov = [2 4 5 6 7 8 9 11 13 14 15 17 18 20 21 22 24 25 26 27 28 29 30 31 34 36 37 39 44];
SJs_mov = [2 4 5 6 8 9 11 13 14 17 18 20 21 22 24 25 27 28 29 31 34 36 37 39 44];

%% Loop over the movement conditions:
for mc=1:numel(m_conds)  
    %% get C_brain-heart parameters
    [C_BrainHF{mc}, C_BrainLF{mc}, C_ROI_BrainHF{mc}, C_ROI_BrainLF{mc}] = NVR_C_BrainHeart(m_conds{mc}, c_style, ROI, data_path);  
end

% save BHI variables
save(savename_hf, 'C_BrainHF');
save(savename_lf, 'C_BrainLF');
save(savename_hf_roi, 'C_ROI_BrainHF');
save(savename_lf_roi, 'C_ROI_BrainLF');

%% Prepare data for LMM analysis in R
% load BHI variables
load(savename_hf_roi, 'C_ROI_BrainHF');
load(savename_lf_roi, 'C_ROI_BrainLF');

design = [];
mov_idx = 2; % Position mov condition in C_BrainHF or m_conds
nomov_idx = 1; % Position nomov condition in C_BrainHF or m_conds

for i = 1:length(SJs_nomov)
    sj = find(SJs_mov==SJs_nomov(i)); % get index SJ_nomov in mov cond
    if isempty(sj) % If no mov cond, take only nomov cond
        mov_conds = {'nomov'};
    else
        mov_conds = m_conds;
    end
    for m = 1:length(mov_conds)
        if strcmp(mov_conds{m},'mov')
            sj_idx = sj;
            m_idx = mov_idx;
        else
            sj_idx = i;
            m_idx = nomov_idx;
        end
        for a = 1:length(a_conds)
            if strcmp(a_conds{a},'LA')
                % BrainHF
                delta2hf = C_ROI_BrainHF{m_idx}{sj_idx}.LA_Delta2HF_meanROI;
                theta2hf = C_ROI_BrainHF{m_idx}{sj_idx}.LA_Theta2HF_meanROI;
                alpha2hf = C_ROI_BrainHF{m_idx}{sj_idx}.LA_Alpha2HF_meanROI;
                beta2hf = C_ROI_BrainHF{m_idx}{sj_idx}.LA_Beta2HF_meanROI;
                gamma2hf = C_ROI_BrainHF{m_idx}{sj_idx}.LA_Gamma2HF_meanROI;
                hf2delta = C_ROI_BrainHF{m_idx}{sj_idx}.LA_HF2Delta_meanROI;
                hf2theta = C_ROI_BrainHF{m_idx}{sj_idx}.LA_HF2Theta_meanROI;
                hf2alpha = C_ROI_BrainHF{m_idx}{sj_idx}.LA_HF2Alpha_meanROI;
                hf2beta = C_ROI_BrainHF{m_idx}{sj_idx}.LA_HF2Beta_meanROI;
                hf2gamma = C_ROI_BrainHF{m_idx}{sj_idx}.LA_HF2Gamma_meanROI;
                % BrainLF
                delta2lf = C_ROI_BrainLF{m_idx}{sj_idx}.LA_Delta2LF_meanROI;
                theta2lf = C_ROI_BrainLF{m_idx}{sj_idx}.LA_Theta2LF_meanROI;
                alpha2lf = C_ROI_BrainLF{m_idx}{sj_idx}.LA_Alpha2LF_meanROI;
                beta2lf = C_ROI_BrainLF{m_idx}{sj_idx}.LA_Beta2LF_meanROI;
                gamma2lf = C_ROI_BrainLF{m_idx}{sj_idx}.LA_Gamma2LF_meanROI;
                lf2delta = C_ROI_BrainLF{m_idx}{sj_idx}.LA_LF2Delta_meanROI;
                lf2theta = C_ROI_BrainLF{m_idx}{sj_idx}.LA_LF2Theta_meanROI;
                lf2alpha = C_ROI_BrainLF{m_idx}{sj_idx}.LA_LF2Alpha_meanROI;
                lf2beta = C_ROI_BrainLF{m_idx}{sj_idx}.LA_LF2Beta_meanROI;
                lf2gamma = C_ROI_BrainLF{m_idx}{sj_idx}.LA_LF2Gamma_meanROI;
            elseif strcmp(a_conds{a},'HA')
                % BrainHF
                delta2hf = C_ROI_BrainHF{m_idx}{sj_idx}.HA_Delta2HF_meanROI;
                theta2hf = C_ROI_BrainHF{m_idx}{sj_idx}.HA_Theta2HF_meanROI;
                alpha2hf = C_ROI_BrainHF{m_idx}{sj_idx}.HA_Alpha2HF_meanROI;
                beta2hf = C_ROI_BrainHF{m_idx}{sj_idx}.HA_Beta2HF_meanROI;
                gamma2hf = C_ROI_BrainHF{m_idx}{sj_idx}.HA_Gamma2HF_meanROI;
                hf2delta = C_ROI_BrainHF{m_idx}{sj_idx}.HA_HF2Delta_meanROI;
                hf2theta = C_ROI_BrainHF{m_idx}{sj_idx}.HA_HF2Theta_meanROI;
                hf2alpha = C_ROI_BrainHF{m_idx}{sj_idx}.HA_HF2Alpha_meanROI;
                hf2beta = C_ROI_BrainHF{m_idx}{sj_idx}.HA_HF2Beta_meanROI;
                hf2gamma = C_ROI_BrainHF{m_idx}{sj_idx}.HA_HF2Gamma_meanROI;
                % BrainLF
                delta2lf = C_ROI_BrainLF{m_idx}{sj_idx}.HA_Delta2LF_meanROI;
                theta2lf = C_ROI_BrainLF{m_idx}{sj_idx}.HA_Theta2LF_meanROI;
                alpha2lf = C_ROI_BrainLF{m_idx}{sj_idx}.HA_Alpha2LF_meanROI;
                beta2lf = C_ROI_BrainLF{m_idx}{sj_idx}.HA_Beta2LF_meanROI;
                gamma2lf = C_ROI_BrainLF{m_idx}{sj_idx}.HA_Gamma2LF_meanROI;
                lf2delta = C_ROI_BrainLF{m_idx}{sj_idx}.HA_LF2Delta_meanROI;
                lf2theta = C_ROI_BrainLF{m_idx}{sj_idx}.HA_LF2Theta_meanROI;
                lf2alpha = C_ROI_BrainLF{m_idx}{sj_idx}.HA_LF2Alpha_meanROI;
                lf2beta = C_ROI_BrainLF{m_idx}{sj_idx}.HA_LF2Beta_meanROI;
                lf2gamma = C_ROI_BrainLF{m_idx}{sj_idx}.HA_LF2Gamma_meanROI;
            end
            for s = 1:length(delta2hf)
                % IMPORTANT remember m=1 -> nomov; m=2-> mov ; a=1 -> LA; a=2 -> HA
                design = [design; ...
                            SJs_nomov(i) m a ... 
                            delta2hf(s) theta2hf(s) alpha2hf(s) beta2hf(s) gamma2hf(s) ...
                            hf2delta(s) hf2theta(s) hf2alpha(s) hf2beta(s) hf2gamma(s) ...
                            delta2lf(s) theta2lf(s) alpha2lf(s) beta2lf(s) gamma2lf(s) ...
                            lf2delta(s) lf2theta(s) lf2alpha(s) lf2beta(s) lf2gamma(s)];
            end
        end
    end
end

tab = array2table(design,'VariableNames',{'SJ','mov_cond','arousal',...
    'Delta2HF','Theta2HF','Alpha2HF','Beta2HF','Gamma2HF',...
    'HF2Delta','HF2Theta','HF2Alpha','HF2Beta','HF2Gamma',...
    'Delta2LF','Theta2LF','Alpha2LF','Beta2LF','Gamma2LF',...
    'LF2Delta','LF2Theta','LF2Alpha','LF2Beta','LF2Gamma'});

% Recode mov_cond and arousal variables
% m=1 -> nomov; m=2-> mov ; a=1 -> LA; a=2 -> HA
for i = 1:length(tab.arousal)
    if tab.arousal(i) == 1
        aro_label{i} = 'LA';
    elseif tab.arousal(i) == 2
        aro_label{i} = 'HA';
    end
    if tab.mov_cond(i) == 1
        mov_label{i} = 'nomov';
    elseif tab.mov_cond(i) == 2
        mov_label{i} = 'mov';
    end
end
tab2 = removevars(tab,{'arousal','mov_cond'});
tab2 = addvars(tab2,mov_label',aro_label','After','SJ','NewVariableNames',{'mov_cond','arousal'});

% save table
writetable(tab2,csv_filename,'Delimiter',',');
