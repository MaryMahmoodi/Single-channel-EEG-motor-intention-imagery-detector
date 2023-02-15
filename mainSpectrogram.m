function [T, F, P]=mainSpectrogram(Signal,MaxFreq,SignalCoef,FcFilter,fs,show)
 
if nargin<5; show=0;end 
% Calculating Spectrogram of the input signal
% the function should be run as this:
%
% >> SpectrogramDisplay(SIGNAL,MAXFREQ,SIGNALCOEF);
%
% Where 
%
% SIGNAL is the raw data exported from the software
% MAXFREQ is the maximun frequency that is going to be displayed in the spectrogram
% SIGNALCOEF is the coefficient to convert data values to microvolt
% for ECG Port SIGNALCOEF was -1*1e-3*1/4*20
% for other Ports SIGNALCOEF was -1*1e-3*1/4*10
 
%% STFT parameter
FcHighPass = FcFilter(1);
FcLowPass = FcFilter(2);
%win_size = 100;                         % window size
factor1=5;
win_size = 100;%round(fs/factor1);%100,400
ovlp = 99;%(round(fs/factor1) -1) ;%99, 399
%ovlp = 99 ;                             % number of overlap
%nfft = 256 ;                            % number of fft 
nfft = fs;%round(fs/factor1);%100 ; %fftlength
Fs = fs ;                              % sampling rate
beta = 20 ;%20                             % beta parameter in Kaiser window
T = 1/Fs;                               % Sample time
L = length(Signal);                     % Length of signal
t = (0:L-1)*T;                          % Time vector
MaxFreq = MaxFreq*2;
if (MaxFreq>Fs/2)
    MaxFreq = Fs/2;
end
% getd = @(p)path(p,path); % scilab users must *not* execute this
% getd('toolbox_signal/');
% getd('toolbox_general/');
 
Signal = Signal*SignalCoef;
signal1 = Signal;
 
 
%% Filtering Signal
 
% % apply lowpass filter
% d = fdesign.lowpass('Fp,Fst,Ap,As',68,72,1,20,500);
% Hd = design(d,'IIR');
% signal2 = filter(Hd,signal1);
signal2=signal1;
% apply notchfilter
w0 = 50/(Fs/2);
[num,den] = iirnotch(w0,w0/35);
signal3 = filter(num,den,signal2);
y = signal3;
 
%% Low pass filter
fc = FcLowPass ;
Wn = (2/Fs)*fc;
b = fir1(40,Wn,'low',kaiser(41,3)); 
y1 = filter(b,1,y);
 
%% High pass Filter
if(FcHighPass>0)
Wn = FcHighPass/(Fs/2);
[bBH,aBH] = butter(3,Wn,'high');
y1 = filter(bBH,aBH,y1);
end
 
%% Downsampling for improving contrast
rate = round(Fs/(MaxFreq/1));
y1 = decimate(y1,rate) ;
 

%% Spectrogram Calculation and Display
 %figure 
 %subplot(2,1,1)
 %plot(t,Yfiltered)
 xlim([t(1),t(end)]);
 title('Time Domain Signal')
 xlabel('Time (S)')
%subplot(2,1,2)

% spectrogram(y1,kaiser(win_size,beta),ovlp,nfft,Fs/(rate),'yaxis');
 [S, F, T, P]=spectrogram(y1,kaiser(win_size,beta),ovlp,nfft,Fs/(rate));
 a1=find(F>=0);
 if max(F)>20
 b1=find(F>=max(F)-2);
 else
  b1=length(F);
 end   
 a1=a1(1);b1=b1(1);
 F=F(a1:b1);
 P=P(a1:b1,:);
 for i=1:size(P,1);
 P(i,:)=P(i,:).*sqrt(i);
 end
P=20*log10((P)+1e-6);


if show 
imagesc( T, F,P ); colormap(jet); colorbar off;
% % % % % % % % figure;a=imshow(P,[]); colormap(jet); colorbar off;
%  xlabel('Time (s)');ylabel('Frequency (Hz)');
 
%  title('Spectrogram of The Signal');
 xlim([min(T),max(T)]);
colormap(jet)
colorbar off;
set(gca,'YDir','normal');
ylim([min(F),MaxFreq/2]);
 end
end
