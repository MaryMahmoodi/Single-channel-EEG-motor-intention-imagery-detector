# Single-channel-EEG-motor-intention-imagery-detector


# A new method for accurate detection of movement intention from single channel EEG for online BCI
Computer Methods and Programs in Biomedicine Update
Volume 1, 2021, 100027



https://doi.org/10.1016/j.cmpbup.2021.100027



Low frequency readiness potential (RP) is elicited in electroencephalograms (EEGs) as one intends to perform an imagery (IMI) or real movement (RMI). 

While in most brain-computer-interface (BCI) applications the challenge is to classify RPs of different limbs from the given EEG trials, the objective of this study is fast and automatic detection of RPs from the entire single channel EEG signal. 

The proposed algorithm has two threshold blocks based on the nonlinear Teager-Kaiser energy operator (TEO) in the first block and the morphological properties of the RP waveform as constraints in the second block. The performance is strongly influenced by the abrupt energy changes due to transients and artefacts. As the major contribution, the proposed nonlinear convex optimization algorithm enables separation of transients from low frequency components by providing a fast thresholding mechanism. 

Application of the proposed method to Physionet RMI dataset, BCI competitionIV-1 IMI dataset and our own left hand movement datasets of healthy subjects led to true positive rates (TPRs) of 76.58%, 83.85%, and 81.15%, number of FPs/min of 2.41, 1.4, and 1.6 and accuracy rates of 85.43%, 90%, and 91%. Movement onset detection latency from our automatic RP detector was -384.5 ms.

# As a conclusion, the proposed method outperforms state-of-the-art techniques using as low as single channel EEG making it suitable for real-time neuro-rehabilitation of paralyzed subjects suffering from stroke.



 Run the code named "RunSingleChannel_BPdetection_TVDteagerenergy.m"


