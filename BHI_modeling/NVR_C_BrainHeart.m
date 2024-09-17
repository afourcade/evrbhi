function [C_BrainHF_all, C_BrainLF_all, C_ROI_BrainHF_all, C_ROI_BrainLF_all] = NVR_C_BrainHeart(mov_cond, c_style, ROI, path_data)
% This function computes Brain-Heart couplings, using the model
% developped by Catrambone et al. (2019).
% Also select timepoints corresponding to low (LA) and high (HA) emotional
% arousal ratings.
% 
% Inputs:
% mov_cond = head movement conditions (mov, nomov) 
% c_style = sequence of rollercoaster, Space-Break-Andes
% ROI = region of interest of EEG electrodes
% path_data = path to the data (EEG, ECG)
%
% Outputs:
% C_BrainHF_all = all brain-hf_hrv (both directions) couplings for each
%   participant and each EEG electrode
% C_BrainLF_all = all brain-lf_hrv (both directions) couplings for each
%   participant and each EEG electrode
% C_ROI_BrainHF_all = all brain-hf_hrv (both directions) couplings for each
%   participant averaged within EEG ROI 
% C_ROI_BrainLF_all = all brain-lf_hrv (both directions) couplings for each
%   participant averaged within EEG ROI
%
% 2024 Antonin Fourcade

% set paths to the brain, heart and ratings data
path_in_eeg = [path_data 'Frequency_Bands_Power/EEG_pow/' mov_cond '/' c_style '/'];
path_in_hrv = [path_data 'ECG/HRV/' mov_cond '/' c_style '/'];
path_in_ibi = [path_data 'ECG/IBI/' mov_cond '/' c_style '/'];
path_in_ar = [path_data 'ratings/class_bins/' mov_cond '/' c_style '/'];
files_eeg = dir([path_in_eeg '*.mat']);
files_eeg = {files_eeg.name};
files_hrv = dir([path_in_hrv '*.txt']);
files_hrv = {files_hrv.name};
files_ibi = dir([path_in_ibi '*.txt']);
files_ibi = {files_ibi.name};
files_ar = dir([path_in_ar '*.txt']);
files_ar = {files_ar.name};

% initialize Brain-Heart variables
C_BrainHF_all = cell(length(files_eeg),1);
C_BrainLF_all = cell(length(files_eeg),1);
C_ROI_BrainHF_all = cell(length(files_eeg),1);
C_ROI_BrainLF_all = cell(length(files_eeg),1);

for isub = 1:length(files_eeg)
    
    % get brain, heart and ratings data of participant
    filename_eeg = files_eeg{isub};
    EEG = load([path_in_eeg filename_eeg]);
    filename_hrv = files_hrv{isub};
    [times_hrv, hrv_lf, hrv_hf, ~] = readvars([path_in_hrv filename_hrv]);
    filename_ibi = files_ibi{isub};
    [times_ibi, ibi] = readvars([path_in_ibi filename_ibi]);
    filename_ar = files_ar{isub};
    [~, bins_ar] = readvars([path_in_ar filename_ar]);
    
    % get channels ROI
    chans = cellstr(EEG.chans);
    idx_ROI = ismember(chans, ROI);

    % Preparing Brain inputs for SDG model
    TFR_delta = EEG.delta_power_mirrored';
    TFR_theta = EEG.theta_power_mirrored';
    TFR_alpha = EEG.alpha_power_mirrored';
    TFR_beta = EEG.beta_power_mirrored';
    TFR_gamma = EEG.gamma_power_mirrored';

    % Formating Heart inputs for SDG model
    TFR_HRV_hf = hrv_hf';
    TFR_HRV_lf = hrv_lf';
    RR = ibi';
    t_RR = times_ibi';
    time = times_hrv';
    % Model parameters
    Fs = 1;
    win_RR = 15; %s default 15 
    window = win_RR*Fs;
    nb_ch = size(TFR_delta,1);
    
    % BHI modeling for Delta and HF-HRV
    [LF2Delta, HF2Delta, Delta2LF, Delta2HF] = BHImodel(TFR_delta, TFR_HRV_hf, TFR_HRV_lf, time, t_RR, Fs, RR, window);
    % BHI modeling for Theta and HF-HRV
    [LF2Theta, HF2Theta, Theta2LF, Theta2HF] = BHImodel(TFR_theta, TFR_HRV_hf, TFR_HRV_lf, time, t_RR, Fs, RR, window);
    % BHI modeling for Alpha and HF-HRV
    [LF2Alpha, HF2Alpha, Alpha2LF, Alpha2HF] = BHImodel(TFR_alpha, TFR_HRV_hf, TFR_HRV_lf, time, t_RR, Fs, RR, window);
    % BHI modeling for Beta and HF-HRV
    [LF2Beta, HF2Beta, Beta2LF, Beta2HF] = BHImodel(TFR_beta, TFR_HRV_hf, TFR_HRV_lf, time, t_RR, Fs, RR, window);
    % BHI modeling for Gamma and HF-HRV
    [LF2Gamma, HF2Gamma, Gamma2LF, Gamma2HF] = BHImodel(TFR_gamma, TFR_HRV_hf, TFR_HRV_lf, time, t_RR, Fs, RR, window);

    % Add missing values due to model initialization and save in structure
    C_BrainHF_all{isub}.time = time;
    C_BrainHF_all{isub}.chans = chans;
    C_BrainHF_all{isub}.aro_rat = bins_ar;
    C_BrainHF_all{isub}.Delta2HF = [NaN(nb_ch,2*win_RR) Delta2HF];
    C_BrainHF_all{isub}.HF2Delta = [NaN(nb_ch,win_RR) HF2Delta];
    C_BrainHF_all{isub}.Theta2HF = [NaN(nb_ch,2*win_RR) Theta2HF];
    C_BrainHF_all{isub}.HF2Theta = [NaN(nb_ch,win_RR) HF2Theta];
    C_BrainHF_all{isub}.Alpha2HF = [NaN(nb_ch,2*win_RR) Alpha2HF];
    C_BrainHF_all{isub}.HF2Alpha = [NaN(nb_ch,win_RR) HF2Alpha];
    C_BrainHF_all{isub}.Beta2HF = [NaN(nb_ch,2*win_RR) Beta2HF];
    C_BrainHF_all{isub}.HF2Beta = [NaN(nb_ch,win_RR) HF2Beta];
    C_BrainHF_all{isub}.Gamma2HF = [NaN(nb_ch,2*win_RR) Gamma2HF];
    C_BrainHF_all{isub}.HF2Gamma = [NaN(nb_ch,win_RR) HF2Gamma];
    
    C_BrainLF_all{isub}.time = time;
    C_BrainLF_all{isub}.chans = chans;
    C_BrainLF_all{isub}.aro_rat = bins_ar;
    C_BrainLF_all{isub}.Delta2LF = [NaN(nb_ch,2*win_RR) Delta2LF];
    C_BrainLF_all{isub}.LF2Delta = [NaN(nb_ch,win_RR) LF2Delta];
    C_BrainLF_all{isub}.Theta2LF = [NaN(nb_ch,2*win_RR) Theta2LF];
    C_BrainLF_all{isub}.LF2Theta = [NaN(nb_ch,win_RR) LF2Theta];
    C_BrainLF_all{isub}.Alpha2LF = [NaN(nb_ch,2*win_RR) Alpha2LF];
    C_BrainLF_all{isub}.LF2Alpha = [NaN(nb_ch,win_RR) LF2Alpha];
    C_BrainLF_all{isub}.Beta2LF = [NaN(nb_ch,2*win_RR) Beta2LF];
    C_BrainLF_all{isub}.LF2Beta = [NaN(nb_ch,win_RR) LF2Beta];
    C_BrainLF_all{isub}.Gamma2LF = [NaN(nb_ch,2*win_RR) Gamma2LF];
    C_BrainLF_all{isub}.LF2Gamma = [NaN(nb_ch,win_RR) LF2Gamma];
    
    %% Select data in ROI
    ROI_Delta2HF = C_BrainHF_all{isub}.Delta2HF(idx_ROI,:); 
    ROI_HF2Delta = C_BrainHF_all{isub}.HF2Delta(idx_ROI,:);
    ROI_Delta2LF = C_BrainLF_all{isub}.Delta2LF(idx_ROI,:); 
    ROI_LF2Delta = C_BrainLF_all{isub}.LF2Delta(idx_ROI,:);
    ROI_Theta2HF = C_BrainHF_all{isub}.Theta2HF(idx_ROI,:); 
    ROI_HF2Theta = C_BrainHF_all{isub}.HF2Theta(idx_ROI,:);
    ROI_Theta2LF = C_BrainLF_all{isub}.Theta2LF(idx_ROI,:); 
    ROI_LF2Theta = C_BrainLF_all{isub}.LF2Theta(idx_ROI,:);
    ROI_Alpha2HF = C_BrainHF_all{isub}.Alpha2HF(idx_ROI,:); 
    ROI_HF2Alpha = C_BrainHF_all{isub}.HF2Alpha(idx_ROI,:);
    ROI_Alpha2LF = C_BrainLF_all{isub}.Alpha2LF(idx_ROI,:); 
    ROI_LF2Alpha = C_BrainLF_all{isub}.LF2Alpha(idx_ROI,:);
    ROI_Beta2HF = C_BrainHF_all{isub}.Beta2HF(idx_ROI,:); 
    ROI_HF2Beta = C_BrainHF_all{isub}.HF2Beta(idx_ROI,:);
    ROI_Beta2LF = C_BrainLF_all{isub}.Beta2LF(idx_ROI,:); 
    ROI_LF2Beta = C_BrainLF_all{isub}.LF2Beta(idx_ROI,:);
    ROI_Gamma2HF = C_BrainHF_all{isub}.Gamma2HF(idx_ROI,:); 
    ROI_HF2Gamma = C_BrainHF_all{isub}.HF2Gamma(idx_ROI,:);
    ROI_Gamma2LF = C_BrainLF_all{isub}.Gamma2LF(idx_ROI,:); 
    ROI_LF2Gamma = C_BrainLF_all{isub}.LF2Gamma(idx_ROI,:);

    % Average over the ROI 
    Delta2HF_meanROI = mean(ROI_Delta2HF,1);
    HF2Delta_meanROI = mean(ROI_HF2Delta,1);
    Delta2LF_meanROI = mean(ROI_Delta2LF,1);
    LF2Delta_meanROI = mean(ROI_LF2Delta,1);
    Theta2HF_meanROI = mean(ROI_Theta2HF,1);
    HF2Theta_meanROI = mean(ROI_HF2Theta,1);
    Theta2LF_meanROI = mean(ROI_Theta2LF,1);
    LF2Theta_meanROI = mean(ROI_LF2Theta,1);
    Alpha2HF_meanROI = mean(ROI_Alpha2HF,1);
    HF2Alpha_meanROI = mean(ROI_HF2Alpha,1);
    Alpha2LF_meanROI = mean(ROI_Alpha2LF,1);
    LF2Alpha_meanROI = mean(ROI_LF2Alpha,1);
    Beta2HF_meanROI = mean(ROI_Beta2HF,1);
    HF2Beta_meanROI = mean(ROI_HF2Beta,1);
    Beta2LF_meanROI = mean(ROI_Beta2LF,1);
    LF2Beta_meanROI = mean(ROI_LF2Beta,1);
    Gamma2HF_meanROI = mean(ROI_Gamma2HF,1);
    HF2Gamma_meanROI = mean(ROI_HF2Gamma,1);
    Gamma2LF_meanROI = mean(ROI_Gamma2LF,1);
    LF2Gamma_meanROI = mean(ROI_LF2Gamma,1);

    % Trim the mirrored data and save in structure
    idx_begin = find(times_hrv==0);
    idx_end = find(times_hrv==269);

    C_ROI_BrainHF_all{isub}.time = time(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.roi = ROI;
    C_ROI_BrainHF_all{isub}.aro_rat = bins_ar;
    C_ROI_BrainHF_all{isub}.Delta2HF_meanROI = Delta2HF_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.HF2Delta_meanROI = HF2Delta_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.Theta2HF_meanROI = Theta2HF_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.HF2Theta_meanROI = HF2Theta_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.Alpha2HF_meanROI = Alpha2HF_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.HF2Alpha_meanROI = HF2Alpha_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.Beta2HF_meanROI = Beta2HF_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.HF2Beta_meanROI = HF2Beta_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.Gamma2HF_meanROI = Gamma2HF_meanROI(idx_begin:idx_end);
    C_ROI_BrainHF_all{isub}.HF2Gamma_meanROI = HF2Gamma_meanROI(idx_begin:idx_end);

    C_ROI_BrainLF_all{isub}.time = time(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.roi = ROI;
    C_ROI_BrainLF_all{isub}.aro_rat = bins_ar;
    C_ROI_BrainLF_all{isub}.Delta2LF_meanROI = Delta2LF_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.LF2Delta_meanROI = LF2Delta_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.Theta2LF_meanROI = Theta2LF_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.LF2Theta_meanROI = LF2Theta_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.Alpha2LF_meanROI = Alpha2LF_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.LF2Alpha_meanROI = LF2Alpha_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.Beta2LF_meanROI = Beta2LF_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.LF2Beta_meanROI = LF2Beta_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.Gamma2LF_meanROI = Gamma2LF_meanROI(idx_begin:idx_end);
    C_ROI_BrainLF_all{isub}.LF2Gamma_meanROI = LF2Gamma_meanROI(idx_begin:idx_end);

    % Z-score the time-series
    C_ROI_BrainHF_all{isub}.z_Delta2HF_meanROI = zscore(C_ROI_BrainHF_all{isub}.Delta2HF_meanROI);
    C_ROI_BrainHF_all{isub}.z_HF2Delta_meanROI = zscore(C_ROI_BrainHF_all{isub}.HF2Delta_meanROI);
    C_ROI_BrainHF_all{isub}.z_Theta2HF_meanROI = zscore(C_ROI_BrainHF_all{isub}.Theta2HF_meanROI);
    C_ROI_BrainHF_all{isub}.z_HF2Theta_meanROI = zscore(C_ROI_BrainHF_all{isub}.HF2Theta_meanROI);
    C_ROI_BrainHF_all{isub}.z_Alpha2HF_meanROI = zscore(C_ROI_BrainHF_all{isub}.Alpha2HF_meanROI);
    C_ROI_BrainHF_all{isub}.z_HF2Alpha_meanROI = zscore(C_ROI_BrainHF_all{isub}.HF2Alpha_meanROI);
    C_ROI_BrainHF_all{isub}.z_Beta2HF_meanROI = zscore(C_ROI_BrainHF_all{isub}.Beta2HF_meanROI);
    C_ROI_BrainHF_all{isub}.z_HF2Beta_meanROI = zscore(C_ROI_BrainHF_all{isub}.HF2Beta_meanROI);
    C_ROI_BrainHF_all{isub}.z_Gamma2HF_meanROI = zscore(C_ROI_BrainHF_all{isub}.Gamma2HF_meanROI);
    C_ROI_BrainHF_all{isub}.z_HF2Gamma_meanROI = zscore(C_ROI_BrainHF_all{isub}.HF2Gamma_meanROI);

    C_ROI_BrainLF_all{isub}.z_Delta2LF_meanROI = zscore(C_ROI_BrainLF_all{isub}.Delta2LF_meanROI);
    C_ROI_BrainLF_all{isub}.z_LF2Delta_meanROI = zscore(C_ROI_BrainLF_all{isub}.LF2Delta_meanROI);
    C_ROI_BrainLF_all{isub}.z_Theta2LF_meanROI = zscore(C_ROI_BrainLF_all{isub}.Theta2LF_meanROI);
    C_ROI_BrainLF_all{isub}.z_LF2Theta_meanROI = zscore(C_ROI_BrainLF_all{isub}.LF2Theta_meanROI);
    C_ROI_BrainLF_all{isub}.z_Alpha2LF_meanROI = zscore(C_ROI_BrainLF_all{isub}.Alpha2LF_meanROI);
    C_ROI_BrainLF_all{isub}.z_LF2Alpha_meanROI = zscore(C_ROI_BrainLF_all{isub}.LF2Alpha_meanROI);
    C_ROI_BrainLF_all{isub}.z_Beta2LF_meanROI = zscore(C_ROI_BrainLF_all{isub}.Beta2LF_meanROI);
    C_ROI_BrainLF_all{isub}.z_LF2Beta_meanROI = zscore(C_ROI_BrainLF_all{isub}.LF2Beta_meanROI);
    C_ROI_BrainLF_all{isub}.z_Gamma2LF_meanROI = zscore(C_ROI_BrainLF_all{isub}.Gamma2LF_meanROI);
    C_ROI_BrainLF_all{isub}.z_LF2Gamma_meanROI = zscore(C_ROI_BrainLF_all{isub}.LF2Gamma_meanROI);

% DEBUG
% PLOT
%         toplot = [TFR_HRV_hf; TFR_HRV_lf; mean(log10(TFR_alpha(idx_ROI,:)),1); Alpha2HF_meanROI; HF2Alpha_meanROI; Alpha2LF_meanROI; LF2Alpha_meanROI];
%         titles = {'HF-HRV', 'LF-HRV', 'log alpha power', 'Alpha2HF', 'HF2Alpha', 'Alpha2LF', 'LF2Alpha'};
%         for j=1:size(toplot,1)
%             figure
%             plot(times_hrv,toplot(j,:))
%             xline(0, '--r')
%             xline(148, '--b')
%             xline(178, '--b')
%             xline(270, '--r')
%             xlim([times_hrv(1) times_hrv(length(times_hrv))])
%             title(titles{j})
%         end

%     %Plotting z-scores
%     figure, hold all
%     plot(times_pow, C_BrainHF_all{isub}.z_Delta2HF_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_HF2Delta_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_Theta2HF_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_HF2Theta_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_Alpha2HF_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_HF2Alpha_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_Beta2HF_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_HF2Beta_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_Gamma2HF_meanROI);
%     plot(times_pow, C_BrainHF_all{isub}.z_HF2Gamma_meanROI);
%     plot(times_pow, aro_rat-2)
%     xline(148,'--r');
%     xline(178,'--r');
%     title(['ALPHA - ROI - Z-scores - ' filename_eeg])
%     legend('Delta2HF', 'HF2Delta', 'Theta2HF', 'HF2Theta', 'Alpha2HF','HF2Alpha','Beta2HF','HF2Beta','Gamma2HF','HF2Gamma','BreakStart', 'BreakEnd')

    %% Select HA and LA samples
    % BrainHF
    C_ROI_BrainHF_all{isub}.LA_Delta2HF_meanROI = C_ROI_BrainHF_all{isub}.Delta2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_Delta2HF_meanROI = C_ROI_BrainHF_all{isub}.Delta2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_HF2Delta_meanROI = C_ROI_BrainHF_all{isub}.HF2Delta_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_HF2Delta_meanROI = C_ROI_BrainHF_all{isub}.HF2Delta_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_Theta2HF_meanROI = C_ROI_BrainHF_all{isub}.Theta2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_Theta2HF_meanROI = C_ROI_BrainHF_all{isub}.Theta2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_HF2Theta_meanROI = C_ROI_BrainHF_all{isub}.HF2Theta_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_HF2Theta_meanROI = C_ROI_BrainHF_all{isub}.HF2Theta_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_Alpha2HF_meanROI = C_ROI_BrainHF_all{isub}.Alpha2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_Alpha2HF_meanROI = C_ROI_BrainHF_all{isub}.Alpha2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_HF2Alpha_meanROI = C_ROI_BrainHF_all{isub}.HF2Alpha_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_HF2Alpha_meanROI = C_ROI_BrainHF_all{isub}.HF2Alpha_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_Beta2HF_meanROI = C_ROI_BrainHF_all{isub}.Beta2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_Beta2HF_meanROI = C_ROI_BrainHF_all{isub}.Beta2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_HF2Beta_meanROI = C_ROI_BrainHF_all{isub}.HF2Beta_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_HF2Beta_meanROI = C_ROI_BrainHF_all{isub}.HF2Beta_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_Gamma2HF_meanROI = C_ROI_BrainHF_all{isub}.Gamma2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_Gamma2HF_meanROI = C_ROI_BrainHF_all{isub}.Gamma2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_HF2Gamma_meanROI = C_ROI_BrainHF_all{isub}.HF2Gamma_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_HF2Gamma_meanROI = C_ROI_BrainHF_all{isub}.HF2Gamma_meanROI(bins_ar == 3);
    % z-scores
    C_ROI_BrainHF_all{isub}.LA_z_Delta2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Delta2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_Delta2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Delta2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_HF2Delta_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Delta_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_HF2Delta_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Delta_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_Theta2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Theta2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_Theta2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Theta2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_HF2Theta_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Theta_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_HF2Theta_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Theta_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_Alpha2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Alpha2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_Alpha2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Alpha2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_HF2Alpha_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Alpha_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_HF2Alpha_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Alpha_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_Beta2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Beta2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_Beta2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Beta2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_HF2Beta_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Beta_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_HF2Beta_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Beta_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_Gamma2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Gamma2HF_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_Gamma2HF_meanROI = C_ROI_BrainHF_all{isub}.z_Gamma2HF_meanROI(bins_ar == 3);
    C_ROI_BrainHF_all{isub}.LA_z_HF2Gamma_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Gamma_meanROI(bins_ar == 1);
    C_ROI_BrainHF_all{isub}.HA_z_HF2Gamma_meanROI = C_ROI_BrainHF_all{isub}.z_HF2Gamma_meanROI(bins_ar == 3);
    % BrainLF
    C_ROI_BrainLF_all{isub}.LA_Delta2LF_meanROI = C_ROI_BrainLF_all{isub}.Delta2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_Delta2LF_meanROI = C_ROI_BrainLF_all{isub}.Delta2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_LF2Delta_meanROI = C_ROI_BrainLF_all{isub}.LF2Delta_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_LF2Delta_meanROI = C_ROI_BrainLF_all{isub}.LF2Delta_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_Theta2LF_meanROI = C_ROI_BrainLF_all{isub}.Theta2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_Theta2LF_meanROI = C_ROI_BrainLF_all{isub}.Theta2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_LF2Theta_meanROI = C_ROI_BrainLF_all{isub}.LF2Theta_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_LF2Theta_meanROI = C_ROI_BrainLF_all{isub}.LF2Theta_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_Alpha2LF_meanROI = C_ROI_BrainLF_all{isub}.Alpha2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_Alpha2LF_meanROI = C_ROI_BrainLF_all{isub}.Alpha2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_LF2Alpha_meanROI = C_ROI_BrainLF_all{isub}.LF2Alpha_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_LF2Alpha_meanROI = C_ROI_BrainLF_all{isub}.LF2Alpha_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_Beta2LF_meanROI = C_ROI_BrainLF_all{isub}.Beta2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_Beta2LF_meanROI = C_ROI_BrainLF_all{isub}.Beta2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_LF2Beta_meanROI = C_ROI_BrainLF_all{isub}.LF2Beta_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_LF2Beta_meanROI = C_ROI_BrainLF_all{isub}.LF2Beta_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_Gamma2LF_meanROI = C_ROI_BrainLF_all{isub}.Gamma2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_Gamma2LF_meanROI = C_ROI_BrainLF_all{isub}.Gamma2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_LF2Gamma_meanROI = C_ROI_BrainLF_all{isub}.LF2Gamma_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_LF2Gamma_meanROI = C_ROI_BrainLF_all{isub}.LF2Gamma_meanROI(bins_ar == 3);
    % z-scores
    C_ROI_BrainLF_all{isub}.LA_z_Delta2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Delta2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_Delta2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Delta2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_LF2Delta_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Delta_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_LF2Delta_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Delta_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_Theta2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Theta2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_Theta2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Theta2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_LF2Theta_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Theta_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_LF2Theta_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Theta_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_Alpha2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Alpha2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_Alpha2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Alpha2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_LF2Alpha_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Alpha_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_LF2Alpha_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Alpha_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_Beta2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Beta2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_Beta2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Beta2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_LF2Beta_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Beta_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_LF2Beta_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Beta_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_Gamma2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Gamma2LF_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_Gamma2LF_meanROI = C_ROI_BrainLF_all{isub}.z_Gamma2LF_meanROI(bins_ar == 3);
    C_ROI_BrainLF_all{isub}.LA_z_LF2Gamma_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Gamma_meanROI(bins_ar == 1);
    C_ROI_BrainLF_all{isub}.HA_z_LF2Gamma_meanROI = C_ROI_BrainLF_all{isub}.z_LF2Gamma_meanROI(bins_ar == 3);

end

end