import mne
import scipy
import glob
import numpy as np
import pandas as pd
import fooof
from scipy.stats import zscore
from scipy.io import loadmat, savemat
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def castToList(x):
    # Casts x to a list
    if isinstance(x, list):
        return x
    elif isinstance(x, str):
        return [x]
    try:
        return list(x)
    except TypeError:
        return [x]


def integrate_pow(frs, PSD, band, method='trapezoid'):
    # Integrate power spectral density (PSD) over frequency band, using trapezoid or simpson rule
    #
    # Inputs:
    # frs (array): Array of frequencies.
    # PSD (array): Power spectral density. PSD shape: (nb_ch, nb_fr)
    # band (list): Frequency band of interest.
    # method (str, optional): Integration rule ('trapezoid' or 'simpson'). Default is 'trapezoid'.
    #
    # Outputs:
    # power (array): Integrated power over the frequency band.
    #
    # 2024 by Antonin Fourcade

    power = np.zeros(PSD.shape[0])
    frs_band = frs[np.logical_and(frs >= band[0], frs <= band[1])]
    for i in range(PSD.shape[0]):
        psd = PSD[i]
        psd_band = psd[np.logical_and(frs >= band[0], frs <= band[1])]
        if method == 'trapezoid':
            # The trapezoidal rule approximates the function as a straight line between adjacent points
            power[i] = scipy.integrate.trapezoid(psd_band, frs_band)
        elif method == 'simpson':
            power[i] = scipy.integrate.simpson(psd_band, frs_band)
    return power


def eeg_pow_extract(path_in_eeg, mon_path, path_out_pwr, path_out_pwr_roi, roi, fs, bands, freqs, len_win, overlap, 
                    hpf_freq = 0.3, lpf_freq = 45, mirror_len=80, mirror_break=30, cut=35):
    # processes EEG data files to extract and compute power spectral densities (PSDs) for various frequency bands. 
    # It filters the data,
    # performs time-frequency analysis using continuous wavelet transform (CWT), 
    # mirrors the data to handle edge effects, 
    # and averages the power over specified regions of interest (ROI).
    # The results are saved into .mat files. 
    # 
    # Inputs:
    # path_in_eeg (str): Path to the directory containing EEG .set files.
    # mon_path (str): Path to the montage file for EEG channel locations.
    # path_out_pwr (str): Path to the directory where the power spectral density results will be saved.
    # path_out_pwr_roi (str): Path to the directory where the ROI-averaged power spectral density results will be saved.
    # roi (list): List of channels defining the region of interest.
    # fs (int): Sampling frequency of the EEG data.
    # bands (object): Object containing frequency band definitions (delta, theta, alpha, beta, gamma).
    # freqs (array): Array of frequencies for the wavelet transform.
    # len_win (int): Length of the window for smoothing in seconds.
    # overlap (float): Overlap fraction for the smoothing window.
    # hpf_freq (float, optional): High-pass filter frequency in Hz. Default is 0.3.
    # lpf_freq (float, optional): Low-pass filter frequency in Hz. Default is 45.
    # mirror_len (int, optional): Length of the mirrored data at the beginning and end for symmetric padding. Default is 80.
    # mirror_break (int, optional): Length of the mirrored data for the break section. Default is 30.
    # cut (int, optional): Number of samples to cut from the start and end of the EEG data to remove mirror padding. Default is 35.
    #
    # 2024 by Antonin Fourcade

    # Set tfr parameters
    win = len_win * fs
    noverlap = int(win * overlap)

    # Get eeg files
    files_eeg = glob.glob(path_in_eeg + '*.set')

    # Loop for all subjects
    for isub, filename_eeg in enumerate(files_eeg):
        # Import preprocessed EEG from EEGLAB
        montage = mne.channels.read_custom_montage(mon_path)
        ch_types = ['eeg'] * 33
        ch_types[4] = 'eog'
        ch_types[26] = 'eog'
        ch_types[32] = 'ecg'
        info_eeg = mne.create_info(ch_names=montage.ch_names, sfreq=fs, ch_types=ch_types)
        channel_map = {A: B for A, B in zip(montage.ch_names, ch_types)}
        preprocessed = mne.io.read_raw_eeglab(filename_eeg, eog=['HEOG', 'VEOG'], preload=True, verbose=False, uint16_codec=None)
        preprocessed.info = info_eeg
        preprocessed.set_channel_types(channel_map)
        preprocessed.set_montage(montage)

        # Band-pass filter the data (EEGLAB data not [low-pass] filtered)
        preproc_filtered = preprocessed.filter(l_freq=hpf_freq, h_freq=lpf_freq, verbose=False)

        # pick EEG channels only and get data
        eeg = preproc_filtered.pick_types(eeg=True)
        eeg_data = eeg.get_data()

        # separate S,B,A sections in the data
        eeg_data_S = eeg_data[:, 0:148*fs]
        eeg_data_B = eeg_data[:, 148*fs:178*fs]
        eeg_data_A = eeg_data[:, 178*fs:]
        # Mirror data at beginning and end, symmetric padding
        eeg_data_S_mirrored = np.empty((eeg_data_S.shape[0], eeg_data_S.shape[1] + 2 * mirror_len * fs))
        eeg_data_B_mirrored = np.empty((eeg_data_B.shape[0], eeg_data_B.shape[1] + 2 * mirror_break * fs))
        eeg_data_A_mirrored = np.empty((eeg_data_A.shape[0], eeg_data_A.shape[1] + 2 * mirror_len * fs))
        # Loop for all electrodes
        for e in range(0, eeg_data.shape[0]):
            eeg_data_S_mirrored[e, :] = np.hstack(
                (np.flip(eeg_data_S[e])[-mirror_len * fs:], eeg_data_S[e], np.flip(eeg_data_S[e][-mirror_len * fs:])))
            eeg_data_B_mirrored[e, :] = np.hstack(
                (np.flip(eeg_data_B[e])[-mirror_break * fs:], eeg_data_B[e], np.flip(eeg_data_B[e][-mirror_break * fs:])))
            eeg_data_A_mirrored[e, :] = np.hstack(
                (np.flip(eeg_data_A[e])[-mirror_len * fs:], eeg_data_A[e], np.flip(eeg_data_A[e][-mirror_len * fs:])))

        # TFR computation
        # Continuous wavelet transform (CWT) for each S, B, A sections separately
        Ws = mne.time_frequency.morlet(fs, freqs, n_cycles=7.0, sigma=None, zero_mean=False)
        tfr_S = mne.time_frequency.tfr.cwt(eeg_data_S_mirrored, Ws, use_fft=True, mode='same', decim=1)
        tfr_S_pow = (np.abs(tfr_S)) ** 2
        tfr_B = mne.time_frequency.tfr.cwt(eeg_data_B_mirrored, Ws, use_fft=True, mode='same', decim=1)
        tfr_B_pow = (np.abs(tfr_B)) ** 2
        tfr_A = mne.time_frequency.tfr.cwt(eeg_data_A_mirrored, Ws, use_fft=True, mode='same', decim=1)
        tfr_A_pow = (np.abs(tfr_A)) ** 2

        # Cut the TFR and stack it together again
        tfr_S_pow = tfr_S_pow[:, :, 0:mirror_len * fs + eeg_data_S.shape[1]]
        tfr_B_pow = tfr_B_pow[:, :, mirror_break * fs:mirror_break * fs + eeg_data_B.shape[1]]
        tfr_A_pow = tfr_A_pow[:, :, mirror_len * fs:2 * mirror_len * fs + eeg_data_A.shape[1]]
        tfr_pow_mirrored = np.concatenate((tfr_S_pow, tfr_B_pow, tfr_A_pow), axis=2)

        # Smoothing: Average over 2s with 50% overlap
        PSDs = np.empty(
            (tfr_pow_mirrored.shape[0], tfr_pow_mirrored.shape[1], int(tfr_pow_mirrored.shape[2] / fs) - 1))
        for e in range(0, tfr_pow_mirrored.shape[0]):
            for f in range(0, tfr_pow_mirrored.shape[1]):
                window_avg = [np.mean(tfr_pow_mirrored[e, f, i:i + win]) for i in
                                range(0, len(tfr_pow_mirrored[e, f]), noverlap)
                                if i + win <= len(tfr_pow_mirrored[e, f])]
                PSDs[e, f] = np.asarray(window_avg)

        # Power integration over frequency bands
        # Initialize power variables
        power_delta = np.empty((PSDs.shape[2], PSDs.shape[0],))
        power_theta = np.empty((PSDs.shape[2], PSDs.shape[0],))
        power_alpha = np.empty((PSDs.shape[2], PSDs.shape[0],))
        power_beta = np.empty((PSDs.shape[2], PSDs.shape[0],))
        power_gamma = np.empty((PSDs.shape[2], PSDs.shape[0],))

        # Loop for all time points
        for t in range(0, PSDs.shape[2]):
            PSD = PSDs[:, :, t]
            # Integrate full PSD over the different frequency bands
            power_delta[t] = integrate_pow(freqs, PSD, bands.delta, 'trapezoid')
            power_theta[t] = integrate_pow(freqs, PSD, bands.theta, 'trapezoid')
            power_alpha[t] = integrate_pow(freqs, PSD, bands.alpha, 'trapezoid')
            power_beta[t] = integrate_pow(freqs, PSD, bands.beta, 'trapezoid')
            power_gamma[t] = integrate_pow(freqs, PSD, bands.gamma, 'trapezoid')

        # Cut mirrored data +/- cut s (data w/o artifacts, similar to ECG pipeline)
        power_delta_mirrored = power_delta[mirror_len - 1 - cut:270 + mirror_len - 1 + cut]
        power_theta_mirrored = power_theta[mirror_len - 1 - cut:270 + mirror_len - 1 + cut]
        power_alpha_mirrored = power_alpha[mirror_len - 1 - cut:270 + mirror_len - 1 + cut]
        power_beta_mirrored = power_beta[mirror_len - 1 - cut:270 + mirror_len - 1 + cut]
        power_gamma_mirrored = power_gamma[mirror_len - 1 - cut:270 + mirror_len - 1 + cut]
        times_mirrored = np.arange(0, power_delta_mirrored.shape[0]) - cut

        # Save frequency bands' powers into .mat files
        filename_pwr = filename_eeg.rsplit('\\')[1].rsplit('_')[1] + '_eeg_pow.mat'
        d = {'times': times_mirrored, 'chans': eeg.ch_names, 'delta_power_mirrored': power_delta_mirrored,
             'theta_power_mirrored': power_theta_mirrored,
             'alpha_power_mirrored': power_alpha_mirrored,
             'beta_power_mirrored': power_beta_mirrored, 'gamma_power_mirrored': power_gamma_mirrored}
        savemat(path_out_pwr + filename_pwr, d)
        print(filename_eeg.rsplit('\\')[1].rsplit('_')[1] + ' done ...')

        # Average over ROI
        idx_roi = np.isin(eeg.ch_names, roi)
        roi_delta = power_delta_mirrored[:, idx_roi]
        roi_theta = power_theta_mirrored[:, idx_roi]
        roi_alpha = power_alpha_mirrored[:, idx_roi]
        roi_beta = power_beta_mirrored[:, idx_roi]
        roi_gamma = power_gamma_mirrored[:, idx_roi]
        delta_mean_roi = np.mean(roi_delta, 1)
        theta_mean_roi = np.mean(roi_theta, 1)
        alpha_mean_roi = np.mean(roi_alpha, 1)
        beta_mean_roi = np.mean(roi_beta, 1)
        gamma_mean_roi = np.mean(roi_gamma, 1)

        # Save ROI frequency bands' powers into mat files
        filename_pwr_roi = filename_eeg.rsplit('\\')[1].rsplit('_')[1] + '_eeg_pow_roi.mat'
        d = {'times': times_mirrored, 'roi': roi, 'delta_power_roi_mirrored': delta_mean_roi,
             'theta_power_roi_mirrored': theta_mean_roi,
             'alpha_power_roi_mirrored': alpha_mean_roi,
             'beta_power_roi_mirrored': beta_mean_roi, 'gamma_power_roi_mirrored': gamma_mean_roi}
        savemat(path_out_pwr_roi + filename_pwr_roi, d)
        print(filename_eeg.rsplit('\\')[1].rsplit('_')[1] + ' done ...')


def NVR_eeg_features(path_in_eeg_pow, path_in_ar, roi, cut=35):
    # processes EEG power data and arousal ratings for multiple subjects. 
    # It extracts and z-scores various EEG frequency bands, selects high arousal (HA) and low arousal (LA) samples, 
    # and stores the results in a dictionary for each subject. 
    # 
    # Inputs:
    # path_in_eeg_pow (str): Path to the directory containing EEG power .mat files.
    # path_in_ar (str): Path to the directory containing arousal ratings .txt files.
    # roi (str): Region of interest (electrodes) for EEG data.
    # cut (int, optional): Number of samples to cut from the start and end of the EEG data to remove mirror padding. Default is 35.
    #
    # Outputs:
    # EEGpow_all (list): List of dictionaries containing EEG power data and arousal ratings for each subject.
    #
    # 2024 by Antonin Fourcade

    # Get data files
    files_eeg_pow = glob.glob(path_in_eeg_pow + '*.mat')
    files_ar = glob.glob(path_in_ar + '*.txt')

    # Initialize EEGpow_all
    EEGpow_all = []

    # Loop for all subjects
    for isub, file in enumerate(files_eeg_pow):
        # Set filename:
        filename_eeg_pow = files_eeg_pow[isub]
        filename_ar = files_ar[isub]

        # Import EEG and arousal ratings data
        eeg_pow = loadmat(filename_eeg_pow, simplify_cells=True)
        aro_rat = pd.read_csv(filename_ar, header=0, names=["latency", "class_aro"])

        # Get the different EEG frequency bands power
        # only 269 values for times, because 2s averaging window -> miss 1 sample at the end
        # -> add last value of times pow
        times_pow = np.append(eeg_pow['times'], 270)
        delta_meanROI = eeg_pow['delta_power_roi_mirrored'].T
        theta_meanROI = eeg_pow['theta_power_roi_mirrored'].T
        alpha_meanROI = eeg_pow['alpha_power_roi_mirrored'].T
        beta_meanROI = eeg_pow['beta_power_roi_mirrored'].T
        gamma_meanROI = eeg_pow['gamma_power_roi_mirrored'].T

        # Remove mirror padding
        delta_meanROI = delta_meanROI[cut:270 + cut]
        theta_meanROI = theta_meanROI[cut:270 + cut]
        alpha_meanROI = alpha_meanROI[cut:270 + cut]
        beta_meanROI = beta_meanROI[cut:270 + cut]
        gamma_meanROI = gamma_meanROI[cut:270 + cut]

        # Z-score the time-series
        z_delta_meanROI = zscore(delta_meanROI, nan_policy='omit')
        z_theta_meanROI = zscore(theta_meanROI, nan_policy='omit')
        z_alpha_meanROI = zscore(alpha_meanROI, nan_policy='omit')
        z_beta_meanROI = zscore(beta_meanROI, nan_policy='omit')
        z_gamma_meanROI = zscore(gamma_meanROI, nan_policy='omit')

        # Select the HA and LA samples
        LA_delta_meanROI = delta_meanROI[aro_rat['class_aro'] == 1]
        LA_theta_meanROI = theta_meanROI[aro_rat['class_aro'] == 1]
        LA_alpha_meanROI = alpha_meanROI[aro_rat['class_aro'] == 1]
        LA_beta_meanROI = beta_meanROI[aro_rat['class_aro'] == 1]
        LA_gamma_meanROI = gamma_meanROI[aro_rat['class_aro'] == 1]

        HA_delta_meanROI = delta_meanROI[aro_rat['class_aro'] == 3]
        HA_theta_meanROI = theta_meanROI[aro_rat['class_aro'] == 3]
        HA_alpha_meanROI = alpha_meanROI[aro_rat['class_aro'] == 3]
        HA_beta_meanROI = beta_meanROI[aro_rat['class_aro'] == 3]
        HA_gamma_meanROI = gamma_meanROI[aro_rat['class_aro'] == 3]

        # Select the HA and LA samples from the z-scores
        LA_z_delta_meanROI = z_delta_meanROI[aro_rat['class_aro'] == 1]
        LA_z_theta_meanROI = z_theta_meanROI[aro_rat['class_aro'] == 1]
        LA_z_alpha_meanROI = z_alpha_meanROI[aro_rat['class_aro'] == 1]
        LA_z_beta_meanROI = z_beta_meanROI[aro_rat['class_aro'] == 1]
        LA_z_gamma_meanROI = z_gamma_meanROI[aro_rat['class_aro'] == 1]

        HA_z_delta_meanROI = z_delta_meanROI[aro_rat['class_aro'] == 3]
        HA_z_theta_meanROI = z_theta_meanROI[aro_rat['class_aro'] == 3]
        HA_z_alpha_meanROI = z_alpha_meanROI[aro_rat['class_aro'] == 3]
        HA_z_beta_meanROI = z_beta_meanROI[aro_rat['class_aro'] == 3]
        HA_z_gamma_meanROI = z_gamma_meanROI[aro_rat['class_aro'] == 3]

        # Store relevant variables in eeg_dict
        eeg_dict = {'filename': filename_eeg_pow,
                    'aro_ratings_bins': aro_rat['class_aro'].to_numpy(),
                    'aro_ratings_latency': aro_rat['latency'].to_numpy(),
                    'times': times_pow,
                    'roi': roi,
                    'delta_meanROI': delta_meanROI,
                    'theta_meanROI': theta_meanROI,
                    'alpha_meanROI': alpha_meanROI,
                    'beta_meanROI': beta_meanROI,
                    'gamma_meanROI': gamma_meanROI,
                    'z_delta_meanROI': z_delta_meanROI,
                    'z_theta_meanROI': z_theta_meanROI,
                    'z_alpha_meanROI': z_alpha_meanROI,
                    'z_beta_meanROI': z_beta_meanROI,
                    'z_gamma_meanROI': z_gamma_meanROI,
                    'LA_delta_meanROI': LA_delta_meanROI,
                    'LA_theta_meanROI': LA_theta_meanROI,
                    'LA_alpha_meanROI': LA_alpha_meanROI,
                    'LA_beta_meanROI': LA_beta_meanROI,
                    'LA_gamma_meanROI': LA_gamma_meanROI,
                    'HA_delta_meanROI': HA_delta_meanROI,
                    'HA_theta_meanROI': HA_theta_meanROI,
                    'HA_alpha_meanROI': HA_alpha_meanROI,
                    'HA_beta_meanROI': HA_beta_meanROI,
                    'HA_gamma_meanROI': HA_gamma_meanROI,
                    'LA_z_delta_meanROI': LA_z_delta_meanROI,
                    'LA_z_theta_meanROI': LA_z_theta_meanROI,
                    'LA_z_alpha_meanROI': LA_z_alpha_meanROI,
                    'LA_z_beta_meanROI': LA_z_beta_meanROI,
                    'LA_z_gamma_meanROI': LA_z_gamma_meanROI,
                    'HA_z_delta_meanROI': HA_z_delta_meanROI,
                    'HA_z_theta_meanROI': HA_z_theta_meanROI,
                    'HA_z_alpha_meanROI': HA_z_alpha_meanROI,
                    'HA_z_beta_meanROI': HA_z_beta_meanROI,
                    'HA_z_gamma_meanROI': HA_z_gamma_meanROI,
                    }
        EEGpow_all.append(eeg_dict)
        print(filename_eeg_pow.rsplit('\\')[1].rsplit('_')[0] + ' done ...')
    return EEGpow_all


