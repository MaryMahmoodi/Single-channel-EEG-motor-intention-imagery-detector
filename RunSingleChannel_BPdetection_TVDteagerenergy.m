%% RP detection accroding to the following paper:

% Mahmoodi, M.; Makkiabadi, B.; Mahmoudi, M.; Sanei, S.
% A New Method for Accurate Detection of Movement Intention from Single Channel EEG for Online BCI.



%% 1- input and initial parameters
clc
clear
subjectnumber=2;
coefstd=0.05; %threshold of coefficient for standard deviation of TEO     % 0.1 0.15 0.25


BCI_compet=0;% if using BCI_competitionIV datasets  BCI_compet=1
physionet=0;%if using phyionet data physionet=1
% else make them zero to consider our database of voluntary left hand movements (voluntary hand movement of pushbutton pressing)

duration1=2*60; % (seconds) duration of signal for RP detection

duration=0.25;% (s) lower bound of RP duration (0.3 s or 0.25) for morphological consideration

th=5;% 5uv lowerband of slope
thUP=45;% 45uv upper bound of slope

threshold_EMG=180; 
step_EMG=2;

threshold_EMG2=80; %after blinking rejection
step_EMG2=0.25;

maxfreq=30;% for spectrogram
freqrange=[1,30];

useDWT=0;
usetemplatematching=0;
% if you want to increase the SNR of signal you can calculate the correlation of signal
% with the EEG template
landau=70;




%% load input data

if BCI_compet==1
    
    fs=100;
%     cd('D:\code_electronicdevice\imRPdetection_TEO\data_nbml_movement_BCICIV_MI\BCICIV_1_mat');
    load BCICIV_calib_ds1b
%     cd('D:\code_electronicdevice\imRPdetection_TEO\Mahmoodi2021CMPB_docs_figuresmfile');
    
elseif physionet
    
%     cd('D:\code_electronicdevice\imRPdetection_TEO\data_nbml_movement_BCICIV_MI\physionet');
    [hdr, record] = edfread( ['S00',num2str(subjectnumber),'R03.edf'], 'assignToVariables',true);
    load allRPsamplemarker_durarion
%     cd('D:\code_electronicdevice\imRPdetection_TEO\Mahmoodi2021CMPB_docs_figuresmfile');
    fs=hdr.frequency (1,11);%Cz
    
    
else
    
    fs=256;
    %     cd('D:\code_electronicdevice\imRPdetection_TEO\data_nbml_movement_BCICIV_MI')
    load('EEG_label.mat');
    load (['EEG', num2str(subjectnumber)]);%EEG% s1_l1 l2 template
    load (['markernumber',num2str(subjectnumber)]);%markernumber%
    %     cd('D:\code_electronicdevice\imRPdetection_TEO\Mahmoodi2021CMPB_docs_figures-mfile');
    
end




%% 2- Preprocessing (baseline correction, total variation denoising (TVD), BP filter, myogenic rejection, EOG rejection)
if ~ exist('tvd.mexw64')
    mex tvd.c
end



if BCI_compet==1
    
[Cz,Cz0,markernumber,EEG_marker,fs]=prepare_TVDdatamarkersBCIC (nfo, cnt, mrk, fs,landau,duration1);
   
    
elseif physionet
    
[Cz,Cz0,markernumber,EEG_marker,samplemarker_durarion,fs]=prepare_TVDdatamarkersPhysionet (Cz,allRPsamplemarker_durarion, subjectnumber, landau,fs);
  
else
    
[Cz,Cz0,markernumber,EEG_marker,fs]=prepare_TVDdatamarkers(EEG, markernumber,fs,landau, duration1);
    
   
end



% BP filter and artefact rejection
[Cz1,x_d, x_dLP, timescale,spikes,spikes_index ]=EEG_preprocessing (Cz,freqrange,fs, step_EMG,threshold_EMG,step_EMG2,threshold_EMG2,useDWT,BCI_compet,physionet);

[TEOsignal,timescale]=make_TEOsignal(x_dLP, fs,markernumber,usetemplatematching);
TEO_x_dLP =TEOsignal;
%Teager energy operated signal


%% 3-RP detection with thresholding on TEO signal and morphological consideration on x_dLP

tic,

% % % % % for ROC curve
%   coefstd=[0 0.02 0.03 0.04  0.05 0.06 0.07 0.08 0.09 0.1 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19   0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31 0.32 0.33 ];
% coefstd=[0:0.001:0.3];
%%%%%%%%%%%%%%%

count2=1;
%   for k2=1:length(coefstd)
%   a=  mean(TEOsignal)+(-coefstd(k2))*std(TEOsignal);
% thresh=  mean(TEOsignal)+(-coefstd)*std(TEOsignal);

[EEG_marker_detection2, EEG_marker_detection,markernumber_detection,markernumber_detection2,detectedsamplemarker_duration,primary_binary,binary]=...
RPdetection (TEOsignal,coefstd,fs,x_dLP,spikes_index, th, thUP,EEG_marker);



%% quantify results%%%
if BCI_compet
    samplemarker_durarion=[];
elseif physionet
    samplemarker_durarion;
else
    samplemarker_durarion=[];
 
end

    [sensitivity,  FPR,  accuracy,   dt,  FPperMin,  Fscore]=quantitative_analysis(Cz1, markernumber,markernumber_detection,markernumber_detection2,detectedsamplemarker_duration, fs, BCI_compet,physionet,samplemarker_durarion) ;

Fscore=2*(sensitivity*accuracy)/(sensitivity+accuracy) /100;

%%%%% By-sample  quantification of binary outputs%%%
%EEG_marker_detection;
%EEG_marker;
% [ Score ] = F1score( EEG_marker_detection2, EEG_marker, EEG_marker );
% recall=Score{2}(5);%sensitivity
% precision=Score{2}(6);% somehow accuracy
% sensitivity1=Score{2}(1)/(Score{2}(1)+Score{2}(4))%= recall
% F1score=Score{2}(7) %2*(recall*precision)/(recall+precision)
% Specificity=Score{2}(8)*100
% Accuracy=Score{2}(10)*100;%accuracy

% for ROC curve
sensitivity2(1,count2)=sensitivity;
FPperMin2(1,count2)=FPperMin;
count2=count2+1;
%   end %for ROC curve


if length(coefstd)>1
    figure(7);
    % plot(FPperMin2,sensitivity2,'r'); hold on;
    plot(FPperMin2,sensitivity2,'r');
    % legend(num2str(coefstd))
    xlabel ('Number of FPs/minute'); ylabel(' True Positive Rate(%)');
    title('ROC curve')
    for i=1:length(coefstd)
        text(FPperMin2(i),sensitivity2(i),num2str(coefstd(i)))
    end
    
end

%% plot  results %%
figure(8);clf;
subplot(2,1,1);  [T1, F1, P1]=mainSpectrogram(Cz1,maxfreq,1,freqrange,fs,1);ylim(freqrange);
title('Spectrogram of denoised signal ')
if physionet ;xlim([62 102]);else xlim([1 130]);end

hold on; subplot(2,1,2); %EEG1
plot(timescale, x_dLP,'b');hold on;plot(timescale,60*EEG_marker,'k');% EEG1
text(1,40,'Denoised and 5Hz lowpass filtered signal and movement onsets (black vertical bars)')
threshTEO=(mean(TEOsignal)-coefstd*std(TEOsignal))*ones(size(TEOsignal));
hold on;plot(timescale,15*TEOsignal-20,'g');
text(1,-30,'Teager energy operator on denoised and low pass filtered signal and threshold (black dashed line) ')
hold on;plot(timescale,15*threshTEO-20,'k--');
%   hold on;plot(timescale,10*EEG_marker_detection1-20,'m')
ylabel('Amplitude (\muV)');xlabel('time (s)');
if physionet ;xlim([62 102]);else xlim([1 130]);end

hold on;
plot(timescale,Cz0-80,'b');text(1,-120,'Original signal and spikes/blinks (red) and final movement detections (pink bars)')
if BCI_compet ||  physionet
    hold on; plot(timescale,40*EEG_marker_detection2-80,'m')
else
    hold on; plot(timescale,40*EEG_marker_detection-80,'m')
end
ylim([-150, 70])

num2=find((markernumber)<=length(x_dLP));
markernumber3=markernumber(num2);


template=zeros(1,round(3*fs));
monset=zeros(1,round(3*fs));monset(1,round(2.04*fs))=1;
for i1=2:length(markernumber3)
    template =template+ x_dLP(1,(markernumber3(1,i1)-round(2*fs)):  (markernumber3(1,i1)+round(1*fs))-1      )  ;
end
%out or Cz1

template=template/(i1-1);
time=[1:length(template)]./fs;
time=time-2*ones(size(time));
%   template=template/i1;
figure(3);clf;  hold on; plot(time,template,'b');%[1:length(template)]./fs
hold on;plot(time,10*monset,'k');%[1:length(template)]./fs
hold on; plot(time, 20*T(template),'g')
legend(' RP template','movement onset','TEO_x_d_l_p');xlabel('time(s)');ylabel('Amplitude (\muV)')
title('RP template')
% hold on ;plot(time, 20*T(template),'g')% for S1_R3
th1=ones(1,length(template))* (mean(T(template)));
% hold on; plot(time,100*th1,'k--')

%% RP detection steps
figure(6);clf
hold on;
plot (timescale, Cz0,'b'); hold on;plot(timescale,60*EEG_marker,'k');
text(5,25, 'Raw EEG (Cz) and movement onset markers'); pause(2);
hold on; plot(timescale,x_d-60,'b');
text(5,-50, 'Denoised signal (x_d)');pause(2);
hold on; plot(timescale,x_dLP-100,'b');
text(5, -90, 'Denoised and lowpass filtered signal (x_d_L_P)'); pause(2);

threshTEO=(mean(TEO_x_dLP)+coefstd*std(TEO_x_dLP))*ones(size(TEO_x_dLP));
hold on;plot(timescale,5*TEO_x_dLP-135,'g');
text(5,-130,'TEO applied to x_d_L_P (TEOx_d_L_P) ');pause(2);
hold on;plot(timescale,threshTEO-128,'k--');
hold on; plot(timescale, 5*primary_binary-160,'k');
text (5, -148,'Primary RP detections (after threshold on TEOx_d_L_P)');pause(2);
hold on; plot(timescale, 5*binary-180,'m');
text (5, -174,'Final RP detections (after duration and slope control)');pause(2);

ylabel('Amplitude (\muV)');xlabel('time (s)');
ylim([-180 40])
xlim([5 20])

%% correlation of template with x_dLP signal
template=zeros(1,round(3*fs));
monset=zeros(1,round(3*fs));monset(1,round(2*fs))=1;
for i1=2:length(markernumber3)
    template =template+ x_dLP(1,(markernumber3(1,i1)-round(2*fs)):  (markernumber3(1,i1)+round(1*fs))-1      )  ;
    %out or Cz1
end
template=template/(i1-1);

corr_EEG1=xcorr(template, x_dLP); % corrcoef(template, EEG1(1:length(template)))
corr_EEG1=corr_EEG1(1:length(x_dLP))/(max(x_dLP).^2);
corr_EEG1=T(corr_EEG1); % non-linear teager energy changes 
% % % % % % % % corr_EEG1=abs(hilbert(corr_EEG1));%envelope detection

% figure(9) ; 
% plot(timescale, x_dLP,'b');hold on;plot(timescale,60*EEG_marker,'k');% EEG1
% text(1,40,'Denoised and 5Hz lowpass filtered signal and movement onsets (black vertical bars)')
% thresh=10*(mean(corr_EEG1)+coefstd*std(corr_EEG1))*ones(size(x_dLP));
% hold on;plot(timescale,50*corr_EEG1-20,'g');
% text(1,-30,'signal of correlation of RP and EEG signal (Cz) and threshold (black dashed line) ')
% hold on;plot(timescale,thresh-20,'k--');
% %   hold on;plot(timescale,10*EEG_marker_detection1-20,'m')
% ylabel('Amplitude (\muV)');xlabel('time (s)');
% ylim([-60 70]);
