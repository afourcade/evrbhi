# evrbhi (emotion-vr-brain_heart_interaction)

Analysis scripts for **Fourcade et al., (2024) - Linking brain-heart interactions to emotional arousal in immersive virtual reality**. *Psychophysiology*. DOI:10.1111/psyp.14696

Folders:
* BHI_modeling - contains the brain-heart interaction (BHI) modeling and analysis
* ECG - contains the ECG feature extraction (interbeat interval, IBI; heart rate variability, LF/HF-HRV) and analysis
* EEG_power - contains the EEG frequency band (delta, theta, alpha, beta, gamma) power extraction and analysis
* HEP - contains the EEG preprocessing and heart-evoked potential (HEP) analysis
* Plots - contains plotting scripts for emotional arousal ratings and mean time-series of all features across participants
* SourceReconstruction - contains the HEP source reconstruction
* State_traits_reports - contains the preprocessing and analysis of questionnaires (valence, presence, stai-t, upps)

The EEG preprocessing (HEP folder) should be run first, before EEG_power, ECG or BHI_modeling.
