function [LF2B, HF2B, B2LF, B2HF] = BHImodel(TFR_EEG, TFR_HRV_HF, TFR_HRV_LF, time, t_RR, FS, RR, window)
% This function quantifies directional Brain-Heart Interplay (BHI) 
% through the model proposed by Catrambone et al.(2019) [1].

% INPUT variables:
% TFR_EEG = Time course of EEG power spectral density (PSD). This must be a matrix (Dimension: Channels X time)
%           samples. Each series in each row should be filtered in the desired frequency band of
%           interest (psi)
% TFR_HRV_HF = Time course of HF-HRV PSD (Dimension: 1 X time). 
% TFR_HRV_LF = Time course of LF-HRV PSD (Dimension: 1 X time).
% time = time array for EEG and HRV data (should be the same)
% FS      = Sampling Frequency of the two TFRs
% RR = NON interpolated interbeat intervals in seconds
% t_RR = time of each RR in seconds
% win_RR  = windows length (expressed in seconds) in which the heartbeat generation model (IPFM) is
% reconstructed (default = 15s)
% window  = windows length (in seconds) in which the parameters are
% calculated (default: window*FS >= 15 )

% OUTPUT variables:
% - HeartToBrain = Functional coupling index (c_rrTOeeg(T)) from 
% HRV Phi-band to EEG Psi-band
% - BrainToHF, BrainToLF  = Functional coupling indices from  
%  EEG Psi-band to  HRV-LF or  HRV-HF bands
% - HeartToBrain_sigma, HeartToBrain_mc = model parameters to be used for fitting evaluation [1]
% 
% This software assumes that input series 
% are all artifact free, e.g., heartbeat dynamics free of algotirhmic and/or physiological artifacts; e.g.
% EEG series free of artifacts from eye blink, movement, etc.
% ---------------------------------------------------------------------------------------------
%  This code implements the theoretical dissertation published in:
%  [1] Catrambone Vincenzo, Alberto Greco, Nicola Vanello, Enzo Pasquale Scilingo,
%  and Gaetano Valenza. "Time-Resolved Directional Brain–Heart Interplay Measurement 
%  Through Synthetic Data Generation Models." 
%  Annals of biomedical engineering 47, no. 6 (2019): 1479-1489.
% ---------------------------------------------------------------------------------------------
% Copyright (C) 2019 Vincenzo Catrambone, Gaetano Valenza
% 
% This program is a free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
% 
% If you use this program in support of published research, please include a
% citation of the reference above. If you use this code in a software package,
% please explicitly inform the end users of this copyright notice and ask them
% to cite the reference above in their published research.
% ---------------------------------------------------------------------------------------------
% To use the software from Matlab or Octave, simply call the BHImodel function in
% the src/folder. Type 'help BHImodel' from the Matlab/Octave command window for
% help on the command's syntax and input/output arguments.
% 
% The software does not come with a GUI. 
% Assuming the series sampled at 1 Hz, with 'x' as the time-varying PSD of a given EEG channel  
% integrated in the theta band (4-8 Hz) and 'y' as the time-varying PSD of an HRV series integrated in 0.04-0.4 Hz,
% and 'z' as the HRV series, the following example performs the BHImodel 
% analysis (with default parameters) and plots the function outcomes:
% 
% Fs = 1;
% [HeartToTheta, ThetaToHF, ThetaToLF] = BHImodel(x,y,Fs,z);
% figure, hold all
% plot(1:length(HeartToTheta))'/Fs, HeartToTheta);
% plot(1:length(ThetaToHF))'/Fs, ThetaToHF);
% plot(1:length(ThetaToHF))'/Fs, ThetaToLF);
% ---------------------------------------------------------------------------------------------

%% checking input variables

if nargin<5
    disp('Three arguments are needed at least: the PSD time course of EEG and HRV signal; and the HRV time series');
    return
elseif nargin >= 5 && nargin < 6
    switch nargin
        case 5
            win_RR = 15; window = ceil(15/FS);
        case 6
            window = ceil(15/FS);
    end 
elseif nargin > 8
    disp('Too many input arguments');
    return
end

[r,c] = size(TFR_HRV_HF);
if (c==1)&&(r~=1)
    TFR_HRV_HF = TFR_HRV_HF';
elseif ((r~=1)&&(c~=1))||(length(size(TFR_HRV_HF))>2)
    error('Time course of HF-HRV PSD must be a row vector: 1 X time points')
end

[r,c] = size(TFR_HRV_LF);
if (c==1)&&(r~=1)
    TFR_HRV_LF = TFR_HRV_LF';
elseif ((r~=1)&&(c~=1))||(length(size(TFR_HRV_LF))>2)
    error('Time course of LF-HRV PSD must be a row vector: 1 X time points')
end

if length(size(TFR_EEG))>2
    error('Time course of EEG PSD must be a 2D matrix: Channels X time points')
end
[Nch,Nt] = size(TFR_EEG);
if (Nt==1)&&(Nch~=1)
    TFR_EEG = TFR_EEG';
    [Nch,Nt] = size(TFR_EEG);
end

if (c~=Nt)
    error('The two PSDs, i.e. of EEG and HRV, must be homologously sampled and related to the same time vector, so they must have the same length.');
end

if (log10(abs(median(RR)))<-1)||(log10(abs(median(RR)))>0.7)
    error('HRV signal measure unit must be in seconds!')
end

if window < 15
    wind = ceil(15/FS);
    window = wind*FS;
    disp(['The time window used for BHI estimation has been modified to ' num2str(window)...
        'secs, as minimum window allowing robust results with the chosen sampling rate']);
end

%% HRV model based on Poincare plot
w_lf = 2*pi*0.1;                % LF mean frequency
w_hf = 2*pi*0.25;               % HF mean Frequency

ss = window;
sc = 1;
nt = ceil((length(time)-ss)/sc);

Cs = zeros(1,nt);
Cp = zeros(1,nt);
TM = zeros(1,nt);

for i = 1:nt
    ix1 = (i-1)*sc + 1;
    ix2 = ix1 + ss - 1;
%     ixm = floor(mean(ix1:ix2)); % middle of the window   
    ixm = ix1; % beginning of the window
    t1 = time(ix1);
    t2 = time(ix2);
    ix = find(t_RR >= t1 & t_RR<= t2);

    mu_ibi = mean(RR(ix));
    mu_hr = 1/mu_ibi;  
    
    G = sin(w_hf/(2*mu_hr))-sin(w_lf/(2*mu_hr)); 
    
    M_11 = sin(w_hf/(2*mu_hr))*w_lf*mu_hr/(sin(w_lf/(2*mu_hr))*4);
    M_12 = -sqrt(2)*w_lf*mu_hr/(8*sin(w_lf/(2*mu_hr)));
    M_21 = -sin(w_lf/(2*mu_hr))*w_hf*mu_hr/(sin(w_hf/(2*mu_hr))*4);
    M_22 = sqrt(2)*w_hf*mu_hr/(8*sin(w_hf/(2*mu_hr)));
    M = [M_11, M_12; M_21, M_22];
    L = max(RR(ix))-min(RR(ix));         
    W = sqrt(2)*max(abs(RR(ix(2:end))-RR(ix(1:end-1))));
    C = 1/G*M*[L; W];
    Cs(i) = C(1);   Cp(i) = C(2);
    TM(i) = time(ixm);
end

%% normalization 
CSr = Cs/std(Cs);  CPr = Cp/std(Cp);
TFR_EEG = sqrt(TFR_EEG);

%% interpolation (edges are extended to avoid extrapolation)
CSr = interp1([TM time(end)], [CSr CSr(end)], time, 'nearest');
CPr = interp1([TM time(end)], [CPr CPr(end)], time, 'nearest');
% CSr = interp1([time(1) TM time(end)], [CSr(1) CSr CSr(end)], time, 'nearest');
% CPr = interp1([time(1) TM time(end)], [CPr(1) CPr CPr(end)], time, 'nearest');
% CSr = interp1([time(1) TM time(end)], [CSr(1) CSr CSr(end)], time, 'spline');
% CPr = interp1([time(1) TM time(end)], [CPr(1) CPr CPr(end)], time, 'spline');
% CSr = CSr';
% CPr = CPr';

%% model running for each EEG channel

parfor ch = 1:Nch
    [LF2B(ch,:), HF2B(ch,:), B2LF(ch,:), B2HF(ch,:)] = SDG(TFR_EEG(ch,:), TFR_HRV_LF, TFR_HRV_HF, CSr, CPr, window);
end

end

function [LF_to_EEG, HF_to_EEG, EEG_to_LF, EEG_to_HF] = SDG(EEG_ch, HRV_LF, HRV_HF, Cs_i, Cp_i, window)

    Nt = length(EEG_ch);
    %% First time window is calculated separately
    for i = 1 : window
        arx_data = iddata(EEG_ch(i:i+window)', HRV_LF(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]);                                                 
        LF_to_EEG(i) = model_eegP.B(2);

        arx_data = iddata(EEG_ch(i:i+window)', HRV_HF(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]);                                                 
        HF_to_EEG(i) = model_eegP.B(2);
        
        pow_eeg(1,i) = mean(EEG_ch(i:i+window));                                 
    end
    
    
    for i = window+1:min([length(Cp_i),Nt-window, length(HRV_LF)-window])

        %% Heart to brain estimation
        arx_data = iddata(EEG_ch(i:i+window)', HRV_LF(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]); 
        LF_to_EEG(i) = model_eegP.B(2); 

        arx_data = iddata(EEG_ch(i:i+window)', HRV_HF(i:i+window)',1); 
        model_eegP = arx(arx_data,[1 1 1]); 
        HF_to_EEG(i) = model_eegP.B(2); 
        
        pow_eeg(1,i) = mean(EEG_ch(i:i+window));

        %% Brain to heart estimation
        if i-window <= length(Cp_i)-window-1
            EEG_to_HF(i-window) = mean((Cp_i(i-window:i))./pow_eeg(i-window:i));
            EEG_to_LF(i-window) = mean((Cs_i(i-window:i))./pow_eeg(i-window:i));
        else

            EEG_to_HF(i-window) = EEG_to_HF(i-window-1);
            EEG_to_LF(i-window) = EEG_to_LF(i-window-1);
        end
    end

end



