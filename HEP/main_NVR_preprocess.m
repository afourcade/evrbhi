%% NeVRo HEP preprocessing pipeline
%
% Runs the full preprocessing pipeline over all data sets by calling the 
% single helper functions. It does this in seperate iterations for (a) the 
% movement ("mov") vs. the non-movement ("nomov") condition as well as 
% (b) either including ("SBA" := Space Coaster + Break + Andes Coaster) or 
% excluding ("SA") the break.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2024 Antonin Fourcade

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fieldtrip initialization (for PREP pipeline)
addpath(genpath('C:\Users\Antonin\Documents\MATLAB\fieldtrip-master\'));
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath('C:\Users\Antonin\Documents\MATLAB\spm12\external\fieldtrip\')); 
ft_defaults

% Code path
NVR_path = genpath('D:\NeVRo\Analysis\new_HEP_filtHP_0_3Hz\');
addpath(NVR_path);

% Data path
data_path = genpath('D:/NeVRo/Data/');
addpath(data_path);

% Parameters
% head movements conditions
m_conds = {'nomov'; 'mov'};
% sequence of rollercoasters to crop
c_styles = 'SBA'; % Space-Break-Andes
% time to trim from beginning and end of each roller coaster (in s)
trim_s = 2.5;
% high-pass filtering (in Hz)
locutoff = 0.3;
% low-pass filtering (in Hz)
hicutoff = 45;
% epoching (in s)
epoch_begin = -0.3;
epoch_end = 0.6;
% epoch baseline removal (in ms)
baseline = [-125 -25];

%% Loop over the movement conditions:
for mc=1:numel(m_conds)       
    %% Downsample to 250Hz and run PREP pipeline 
    % (Bigdely-Shamlo et al., 2015) for standardized preprocessing:
    % Line-noise (50Hz) removal, robust average rereferencing,
    % interpolation noisy channels
    NVR_01_DS_PREP(m_conds{mc}, data_path)
        
    %% Crop relevant parts and trim
    NVR_02_crop(c_styles,m_conds{mc}, data_path, trim_s);

    %% Add individual arousal events to EEG data
    NVR_03_eventsARO(c_styles,m_conds{mc}, data_path);
    
    %% High-pass filtering 0.3Hz
    NVR_04_hpf(c_styles,m_conds{mc}, data_path, locutoff);
    
    %% Add ECG events to EEG dataset
    % Run HEPLAB beforehand to get R-peak events in a .mat file
    NVR_05_addECG(c_styles,m_conds{mc}, data_path);
    
    %% Copy ICA weights from HP-1Hz dataset to this dataset 
    NVR_06_copy_ICA(c_styles,m_conds{mc}, data_path);

    %% Epoching
    NVR_07_epoch(c_styles,m_conds{mc}, data_path, epoch_begin, epoch_end);
    
    %% Flag artifacts in ICA components
    % Manual visual inspection
    % Eye, muscle and cardiac field artifacts  
      
    %% Rejection ICA artifacts components
    NVR_08_rejcomp(c_styles,m_conds{mc}, data_path)
    
    %% Low-Pass filtering 45Hz
    NVR_09_lpf(c_styles,m_conds{mc}, data_path, hicutoff)
      
    %% Remove epoch baseline
    NVR_10_rmbase(c_styles,m_conds{mc}, data_path, baseline)
    
    %% Separate HA and LA epochs
    NVR_11_sepHALA(c_styles,m_conds{mc}, data_path)
    
    %% Convert EEGLAB files into Fieldtrip structure
    % for statistical analysis (cluster-based permutation testing)
    NVR_12_eeglab2fieldrip(c_styles,m_conds{mc}, data_path)
    
end
