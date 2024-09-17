import glob
import neurokit2 as nk
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import zscore

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

def heart_extract(path_in_rpeak, path_out_ibi, path_out_hrv, mov_cond, fs, lf=[0.04, 0.15], hf=[0.15, 0.4], mirror_len=80, cut=35, fig=0):
    # processes R-peak data from ECG recordings to compute Inter-Beat Intervals (IBI) and Heart Rate Variability (HRV) metrics. 
    # It performs resampling, interpolation, and Continuous Wavelet Transform (CWT) to extract low-frequency (LF) and high-frequency (HF) power components. 
    # The function saves the processed IBI and HRV data to the specified output directories.
    #
    # Inputs:
    # path_in_rpeak (str): Path to the directory containing R-peak data files.
    # path_out_ibi (str): Path to the directory to save the IBI data files.
    # path_out_hrv (str): Path to the directory to save the HRV data files.
    # mov_cond (str): Movement condition (e.g., 'mov', 'nomov').
    # fs (int): Sampling frequency for the output ecg features.
    # lf (list, optional): Low-frequency band for HRV analysis. Default is [0.04, 0.15].
    # hf (list, optional): High-frequency band for HRV analysis. Default is [0.15, 0.4].
    # mirror_len (int, optional): Length (s) of mirror padding for CWT. Default is 80.
    # cut (int, optional): Length (s) of mirror padding to keep for the data. Default is 35.
    # fig (int, optional): Flag for generating figures (not used in the function). Default is 0.
    #
    # 2024 by Antonin Fourcade

    
    files_rpeak = glob.glob(path_in_rpeak + '*.txt')

    for isub, file in enumerate(files_rpeak):
        filename_rpeak = files_rpeak[isub].rsplit('\\')[1] # get the SJ nb for saving
        r_peaks = pd.read_csv(files_rpeak[isub], sep='\t', usecols=[0], names=['r_peak'], squeeze=True)
        r_peaks = r_peaks.to_numpy()

        #calculate ibi for each S,B,A sections separately
        r_peaks_S = r_peaks[r_peaks < 148]
        r_peaks_A = r_peaks[r_peaks > 178]
        r_peaks_b = r_peaks[r_peaks > 148]
        r_peaks_B = r_peaks_b[r_peaks_b < 178]
        ibi_S = np.diff(r_peaks_S)
        ibi_A = np.diff(r_peaks_A)
        ibi_B = np.diff(r_peaks_B)
        ibi = np.append(ibi_S, np.append(ibi_B, ibi_A))
        # reformat r_peaks to exclude first R_peak of each S,B,A
        r_peaks_clean = np.append(r_peaks_S[1:], np.append(r_peaks_B[1:], r_peaks_A[1:]))
        # resample ibi to 1Hz
        samples_new_ibi = np.arange(0, 270.25, 1)
        ibi_interp = interp1d(r_peaks_clean, ibi, kind='cubic', bounds_error=False, fill_value=(ibi[0], ibi[-1]))
        rs_ibi = ibi_interp(samples_new_ibi)

        # compute TFR
        # resampling 4Hz
        fs_rs = 4
        samples_new = np.arange(0, 270.25, 1/fs_rs)
        ibi_interp = interp1d(r_peaks_clean, ibi, kind='cubic', bounds_error=False, fill_value=(ibi[0], ibi[-1]))
        rs_ibi2 = ibi_interp(samples_new)

        # CWT
        # mirror data at beginning and end, symmetric padding
        sp_ibi_end = rs_ibi2[-mirror_len*fs_rs:]
        sp_ibi_end = sp_ibi_end[::-1]
        sp_ibi_begin = rs_ibi2[:mirror_len*fs_rs]
        sp_ibi_begin = sp_ibi_begin[::-1]
        rs_ibi2_mirrored = np.concatenate([sp_ibi_begin, rs_ibi2, sp_ibi_end])

        # compute cwt on mirrored data
        freqs, times, tfr = nk.signal_timefrequency(rs_ibi2_mirrored, fs_rs, method="cwt", nfreqbin=50, min_frequency=lf[0], max_frequency=hf[1], show=False)

        # Smoothing: Average over 2s with 50% overlap
        tfr_smooth = np.empty((tfr.shape[0], int(tfr.shape[1] / fs_rs) - 1))
        win = 2
        overlap = 0.5
        win_len = win * fs_rs
        i_win = int(win_len * overlap)
        for f in range(0, tfr.shape[0]):
            window_avg = [np.mean(tfr[f, i:i + win_len]) for i in range(0, len(tfr[f]), i_win)
                            if i + win_len <= len(tfr[f])]
            tfr_smooth[f] = np.asarray(window_avg)
        times = np.arange(0, tfr_smooth.shape[1])
        tfr = tfr_smooth

        # Get LF and HF power -> integrate over frequencies
        power_lf = np.zeros(times.shape[0])
        power_hf = np.zeros(times.shape[0])
        power_hrv = np.zeros(times.shape[0])
        for t, time in enumerate(times): # enumerate because times is float instead of int
            psd = tfr[:, t]
            frs_lf = freqs[np.logical_and(freqs >= lf[0], freqs <= lf[1])]
            frs_hf = freqs[np.logical_and(freqs >= hf[0], freqs <= hf[1])]
            frs_hrv = freqs[np.logical_and(freqs >= lf[0], freqs <= hf[1])]
            psd_lf = psd[np.logical_and(freqs >= lf[0], freqs <= lf[1])]
            psd_hf = psd[np.logical_and(freqs >= hf[0], freqs <= hf[1])]
            psd_hrv = psd[np.logical_and(freqs >= lf[0], freqs <= hf[1])]
            # The trapezoidal rule approximates the function as a straight line between adjacent points
            power_lf[t] = scipy.integrate.trapezoid(psd_lf, frs_lf)
            power_hf[t] = scipy.integrate.trapezoid(psd_hf, frs_hf)
            power_hrv[t] = scipy.integrate.trapezoid(psd_hrv, frs_hrv)

        power_lf = np.append(power_lf, np.nan) # Add one NaN at the end, because only 269 values (0:268)
        power_hf = np.append(power_hf, np.nan) # Add one NaN at the end, because only 269 values (0:268)
        power_hrv = np.append(power_hrv, np.nan) # Add one NaN at the end, because only 269 values (0:268)
        times = np.append(times, times[-1]+1)

        # Plotting
        if fig == 1:
            plt.plot(samples_new_ibi, rs_ibi, c='orange')
            plt.title('rs_IBI 1Hz')
            plt.xlim([0, samples_new_ibi[-1]])
            plt.axvline(x=148, color='grey', linestyle='dotted', linewidth=1)
            plt.axvline(x=178, color='grey', linestyle='dotted', linewidth=1)
            plt.show()
            plt.plot(np.arange(0, rs_ibi2_mirrored.shape[0]), rs_ibi2_mirrored, c='orange')
            plt.title('rs_IBI 4Hz mirrored')
            plt.xlim([0, rs_ibi2_mirrored.shape[0]])
            plt.show()
            for i, hrv in enumerate(['power_lf', 'power_hf']):
                if hrv == 'power_lf':
                    plt.plot(times, power_lf, c='orange')
                    plt.title('Low Frequency HRV')
                elif hrv == 'power_hf':
                    plt.plot(times, power_hf, c='orange')
                    plt.title('High Frequency HRV')
                plt.xlim([0, power_hf.shape[0]])
                plt.ylim([0, 0.14])
                plt.axvline(x=mirror_len-1, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=270+mirror_len-1, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=148+mirror_len-1, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=148+30+mirror_len-1, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=mirror_len-1-cut, color='red', linestyle='dotted', linewidth=1)
                plt.axvline(x=270+mirror_len-1+cut, color='red', linestyle='dotted', linewidth=1)
                plt.ylabel('Time [sec]')
                plt.xlabel('Time [sec]')
                plt.show()

        # Cut mirrored data +/- cut s (data w/o artifacts)
        power_lf_mirrored = power_lf[mirror_len-1-cut:270+mirror_len-1+cut]
        power_hf_mirrored = power_hf[mirror_len-1-cut:270+mirror_len-1+cut]
        power_hrv_mirrored = power_hrv[mirror_len-1-cut:270+mirror_len-1+cut]
        times_mirrored = times[mirror_len-1-cut:270+mirror_len-1+cut]
        times_mirrored = times_mirrored - times_mirrored[0] - cut

        # mirror rs_ibi 1Hz data at beginning and end, symmetric padding +/- cut s
        rs_ibi_sl = rs_ibi[:-1] # remove last dummy value (not important)
        sp_ibi_end = rs_ibi_sl[-cut:]
        sp_ibi_end = sp_ibi_end[::-1]
        sp_ibi_begin = rs_ibi_sl[:cut]
        sp_ibi_begin = sp_ibi_begin[::-1]
        rs_ibi_mirrored = np.concatenate([sp_ibi_begin, rs_ibi_sl, sp_ibi_end])

        #Plotting
        if fig == 1:
            plt.plot(times_mirrored, rs_ibi_mirrored, c='orange')
            plt.title('rs_IBI 1Hz mirrored cut')
            plt.xlim([times_mirrored[0], times_mirrored[-1]])
            plt.axvline(x=0, color='grey', linestyle='dotted', linewidth=1)
            plt.axvline(x=270, color='grey', linestyle='dotted', linewidth=1)
            plt.axvline(x=148, color='grey', linestyle='dotted', linewidth=1)
            plt.axvline(x=178, color='grey', linestyle='dotted', linewidth=1)
            plt.show()
            for i, hrv in enumerate(['power_lf', 'power_hf']):
                if hrv == 'power_lf':
                    plt.plot(times_mirrored, power_lf_mirrored, c='orange')
                    plt.title('clean mirrored Low Frequency HRV')
                elif hrv == 'power_hf':
                    plt.plot(times_mirrored, power_hf_mirrored, c='orange')
                    plt.title('clean mirrored High Frequency HRV')
                plt.xlim([times_mirrored[0], times_mirrored[-1]])
                plt.ylim([0, 0.14])
                plt.axvline(x=0, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=270, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=148, color='grey', linestyle='dotted', linewidth=1)
                plt.axvline(x=178, color='grey', linestyle='dotted', linewidth=1)
                plt.ylabel('Time [sec]')
                plt.xlabel('Time [sec]')
                plt.show()

        #Save IBI and the LF and HF power into files
        filename_hrv = filename_rpeak.rsplit('_')[0] + '_hrv_' + mov_cond + '.txt'
        d = {'times': times_mirrored, 'LF_mirrored': power_lf_mirrored, 'HF_mirrored': power_hf_mirrored, 'HRV_all_mirrored': power_hrv_mirrored}
        df = pd.DataFrame(data=d)
        df.to_csv(path_out_hrv + filename_hrv, sep='\t', na_rep=np.nan, index=False)
        filename_ibi = filename_rpeak.rsplit('_')[0] + '_ibi_' + mov_cond + '.txt'
        d = {'times': times_mirrored, 'rs_ibi_mirrored': rs_ibi_mirrored}
        df = pd.DataFrame(data=d)
        df.to_csv(path_out_ibi + filename_ibi, sep='\t', na_rep=np.nan, index=False)

def NVR_ecg_features(path_in_rpeak, path_in_ibi, path_in_hrv, path_in_ar, cut=35, fig=0):
    # reads ECG data from multiple text files, processes the data to extract ibi and hrv features 
    # (e.g remove mirror padding, z-score, select HA and LA samples) 
    # and returns a list of dictionaries containing these features for each subject.
    #
    # Inputs:
    # path_in_rpeak (str): Path to the directory containing R-peak data files.
    # path_in_ibi (str): Path to the directory containing Inter-Beat Interval (IBI) data files.
    # path_in_hrv (str): Path to the directory containing Heart Rate Variability (HRV) data files.
    # path_in_ar (str): Path to the directory containing arousal ratings data files.
    # cut (int, optional): Number of initial samples to cut from the data. Default is 35.
    # fig (int, optional): Flag for generating figures (not used in the function). Default is 0.
    #
    # Outputs:
    # ECG_all (list): A list of dictionaries, each containing the following keys:
    # filename: The filename of the R-peak data file.
    # r_peaks: Array of R-peak values.
    # aro_ratings_bins: Array of arousal rating classes.
    # aro_ratings_latency: Array of arousal rating latencies.
    # ibi: Array of IBI values.
    # hrv_lf: Array of low-frequency HRV values.
    # hrv_hf: Array of high-frequency HRV values.
    # rs_samples: Array of time samples after cutting.
    # rs_ibi: Array of IBI values after cutting.
    # rs_hrv_lf: Array of low-frequency HRV values after cutting.
    # rs_hrv_hf: Array of high-frequency HRV values after cutting.
    # rs_ratioLFHF: Array of LF/HF ratio values after cutting.
    # z_ibi: Z-transformed IBI values.
    # z_hrv_lf: Z-transformed low-frequency HRV values.
    # z_hrv_hf: Z-transformed high-frequency HRV values.
    # z_ratioLFHF: Z-transformed LF/HF ratio values.
    # LA_ibi: IBI values for low arousal samples.
    # HA_ibi: IBI values for high arousal samples.
    # LA_hrv_lf: Low-frequency HRV values for low arousal samples.
    # HA_hrv_lf: Low-frequency HRV values for high arousal samples.
    # LA_hrv_hf: High-frequency HRV values for low arousal samples.
    # HA_hrv_hf: High-frequency HRV values for high arousal samples.
    # LA_ratioLFHF: LF/HF ratio values for low arousal samples.
    # HA_ratioLFHF: LF/HF ratio values for high arousal samples.
    # LA_z_ibi: Z-transformed IBI values for low arousal samples.
    # HA_z_ibi: Z-transformed IBI values for high arousal samples.
    # LA_z_hrv_lf: Z-transformed low-frequency HRV values for low arousal samples.
    # HA_z_hrv_lf: Z-transformed low-frequency HRV values for high arousal samples.
    # LA_z_hrv_hf: Z-transformed high-frequency HRV values for low arousal samples.
    # HA_z_hrv_hf: Z-transformed high-frequency HRV values for high arousal samples.
    # LA_z_ratioLFHF: Z-transformed LF/HF ratio values for low arousal samples.
    # HA_z_ratioLFHF: Z-transformed LF/HF ratio values for high arousal samples.
    #
    # 2024 by Antonin Fourcade
    
    # Get data files
    files_rpeak = glob.glob(path_in_rpeak + '*.txt')
    files_ibi = glob.glob(path_in_ibi + '*.txt')
    files_ar = glob.glob(path_in_ar + '*.txt')
    files_hrv = glob.glob(path_in_hrv + '*.txt')

    # Initialize ECG_all
    ECG_all = []

    # Loop for all subjects
    for isub, file in enumerate(files_rpeak):
        # Set filename:
        filename_rpeak = files_rpeak[isub]
        filename_ibi = files_ibi[isub]
        filename_ar = files_ar[isub]
        filename_hrv = files_hrv[isub]

        # Import ECG and arousal ratings data
        r_peaks = pd.read_csv(filename_rpeak, sep='\t', usecols=[0], names=['r_peak'], squeeze=True)
        aro_rat = pd.read_csv(filename_ar, header=0, names=["latency", "class_aro"])
        ibi = pd.read_csv(filename_ibi, header=0, sep='\t', names=["rs_ibi_mirrored"])
        hrv = pd.read_csv(filename_hrv, header=0, sep='\t', names=["LF_mirrored", "HF_mirrored", "HRV_mirrored"], na_values='NaN')

        # Get the different ECG features
        r_peaks = r_peaks.to_numpy()
        ibi = ibi['rs_ibi_mirrored'].to_numpy()
        hrv_lf = hrv["LF_mirrored"].to_numpy()
        hrv_hf = hrv["HF_mirrored"].to_numpy()
        times = hrv["times"].to_numpy()
        #remove mirror padding from IBI, LF & HF-HRV
        ibi_cut = ibi[cut:270+cut]
        hrv_lf_cut = hrv_lf[cut:270 + cut]
        hrv_hf_cut = hrv_hf[cut:270 + cut]
        ratio_lfhf_cut = hrv_lf_cut / hrv_hf_cut
        times_cut = times[cut:270 + cut]

        # Select the HA and LA samples
        LA_ibi = ibi_cut[aro_rat['class_aro'] == 1]
        LA_hrv_lf = hrv_lf_cut[aro_rat['class_aro'] == 1]
        LA_hrv_hf = hrv_hf_cut[aro_rat['class_aro'] == 1]
        LA_ratio_lfhf = ratio_lfhf_cut[aro_rat['class_aro'] == 1]

        HA_ibi = ibi_cut[aro_rat['class_aro'] == 3]
        HA_hrv_lf = hrv_lf_cut[aro_rat['class_aro'] == 3]
        HA_hrv_hf = hrv_hf_cut[aro_rat['class_aro'] == 3]
        HA_ratio_lfhf = ratio_lfhf_cut[aro_rat['class_aro'] == 3]

        # Z-transform the time-series
        # nan policy: ‘propagate’ returns nan; ‘omit’ performs the calculations ignoring nan values
        z_ibi = zscore(ibi_cut, nan_policy='omit')
        z_hrv_hf = zscore(hrv_lf_cut, nan_policy='omit')
        z_hrv_lf = zscore(hrv_hf_cut, nan_policy='omit')
        z_ratio_lfhf = zscore(ratio_lfhf_cut, nan_policy='omit')

        # Select the HA and LA samples from the z-scores
        LA_z_ibi = z_ibi[aro_rat['class_aro'] == 1]
        LA_z_hrv_lf = z_hrv_lf[aro_rat['class_aro'] == 1]
        LA_z_hrv_hf = z_hrv_hf[aro_rat['class_aro'] == 1]
        LA_z_ratio_lfhf = z_ratio_lfhf[aro_rat['class_aro'] == 1]

        HA_z_ibi = z_ibi[aro_rat['class_aro'] == 3]
        HA_z_hrv_lf = z_hrv_lf[aro_rat['class_aro'] == 3]
        HA_z_hrv_hf = z_hrv_hf[aro_rat['class_aro'] == 3]
        HA_z_ratio_lfhf = z_ratio_lfhf[aro_rat['class_aro'] == 3]

        # Store relevant variables in ecg_dict
        ecg_dict = {'filename': filename_rpeak,
                    'r_peaks': r_peaks,
                    'aro_ratings_bins': aro_rat['class_aro'].to_numpy(),
                    'aro_ratings_latency': aro_rat['latency'].to_numpy(),
                    'ibi': ibi,
                    'hrv_lf': hrv_lf,
                    'hrv_hf': hrv_hf,
                    'rs_samples': times_cut,
                    'rs_ibi': ibi_cut,
                    'rs_hrv_lf': hrv_lf_cut,
                    'rs_hrv_hf': hrv_hf_cut,
                    'rs_ratioLFHF': ratio_lfhf_cut,
                    'z_ibi': z_ibi,
                    'z_hrv_lf': z_hrv_lf,
                    'z_hrv_hf': z_hrv_hf,
                    'z_ratioLFHF': z_ratio_lfhf,
                    'LA_ibi': LA_ibi,
                    'HA_ibi': HA_ibi,
                    'LA_hrv_lf': LA_hrv_lf,
                    'HA_hrv_lf': HA_hrv_lf,
                    'LA_hrv_hf': LA_hrv_hf,
                    'HA_hrv_hf': HA_hrv_hf,
                    'LA_ratioLFHF': LA_ratio_lfhf,
                    'HA_ratioLFHF': HA_ratio_lfhf,
                    'LA_z_ibi': LA_z_ibi,
                    'HA_z_ibi': HA_z_ibi,
                    'LA_z_hrv_lf': LA_z_hrv_lf,
                    'HA_z_hrv_lf': HA_z_hrv_lf,
                    'LA_z_hrv_hf': LA_z_hrv_hf,
                    'HA_z_hrv_hf': HA_z_hrv_hf,
                    'LA_z_ratioLFHF': LA_z_ratio_lfhf,
                    'HA_z_ratioLFHF': HA_z_ratio_lfhf,
                    }

        ECG_all.append(ecg_dict)
        print(filename_rpeak.rsplit('\\')[1] + ' done ...')
    return ECG_all
