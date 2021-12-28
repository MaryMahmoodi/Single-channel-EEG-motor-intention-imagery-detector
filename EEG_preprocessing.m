function [Cz1,x_d, x_dLP, timescale,spikes,spikes_index ]=EEG_preprocessing (Cz,freqrange,fs, step_EMG,threshold_EMG,step_EMG2,threshold_EMG2, useDWT, BCI_compet,physionet)
%%% Bandpass filter design %%%
if BCI_compet || physionet
% the Parks-McClellan method is used via the ‘remez’ function of MATLAB
rp = 0.01; % Passband ripple
rs = 26; % Stopband ripple
f = freqrange; % Cutoff frequencies
a = [1 0]; % Desired amplitudes
% Compute deviations
dev = [(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)];
[n,fo,ao,w] = remezord(f,a,dev,fs);
B = remez(n,fo,ao,w);
A=1;
% % freqz(B,A);
Cz1=filtfilt(B,A,Cz);

else
d1= fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2,C',50,freqrange(1),freqrange(1)+0.1,freqrange(2),freqrange(2)+0.5,fs);%50,0.01,0.16,45,50.5,fs);%36
Hd1=design(d1,'equiripple');
Cz1=filter(Hd1, Cz);
end

timescale=(1:length(Cz1))./fs;

%%% myogenic rejection
Cz1=myogenic_rejection(Cz1, fs,step_EMG ,threshold_EMG);

%%%  blinking (EOG) artefact rejection%%%
threshold=3*(1/length(Cz1)*sum(abs(Cz1))) ;% 3*

figure(4) ;clf;subplot(3,1,1); plot(timescale, Cz1,'b');
title('x_T_V_D TVD signal and signal of detected EOG triangular artefacts (EOG_x) ')
spikes=zeros(size(Cz1));
xlabel ('time (s)')
[ Cz1, spikes,spikes_index ] = blinkingrejection( Cz1,fs );


hold on; plot(timescale,spikes,'r--' )
title ('x_d denoised signal ')
xlabel('time (s)')
hold on; plot(timescale, threshold*ones(size(spikes)),'k--')
hold on; plot(timescale, -threshold*ones(size(spikes)),'k--')
ylim([-70 70])

legend ('x_T_V_D: TVD and  1-30 Hz bandpass filtered signal', 'EOG_x: signal of blinking artefacts', 'positive threhold on x_T_V_D to reconstruct  EOG_x', 'negative threhold on x_T_V_D to reconstruct  EOG_x  ')


fprintf ('eye blink removed...\n')

subplot(3,1,2);plot(timescale,Cz1,'b');title('after EOG spike rejection');
xlabel('time (s)')

ylabel ('Amplitude (\muV)')

legend ('x_d: x_T_V_D after blinking artefact rejection')

if useDWT
    % wavelet decomposition with the wavelet with most similarity to signal(sym5)
    [C,L] = wavedec(Cz1,3,'sym5');%3 sym5
    Cz1= wrcoef('a',C,L,'sym5',3);%3 sym5
    % %%%  'a' 3 3-->0-6.5
    Cz1=myogenic_rejection(Cz1, fs,step_EMG2 ,threshold_EMG2);
else
    Cz1=myogenic_rejection(Cz1, fs,step_EMG2 ,threshold_EMG2);
end

x_d = Cz1;
%%%low pass filter%%%
[a2,b2]=butter(4,5.5/(fs/2),'low');%5/(fs/2)
x_dLP=filter(a2,b2,Cz1);%%% out or Cz1

subplot(3,1,3); plot(timescale,x_dLP,'b');title(' spikeEOG rejected and lowpass filtered signal')

xlabel('time (s)')
ylabel ('Amplitude (\muV)')
title ('x_d_L_P denoised and 5Hz lowpass filtered signal ')
legend ('x_d_L_P: x_d after 5Hz lowpass filtering ')


end

