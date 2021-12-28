function [EEG_marker_detection2, EEG_marker_detection,markernumber_detection,markernumber_detection2,detectedsamplemarker_duration,primary_binary,binary]=RPdetection (TEOsignal,coefstd,fs,x_dLP,spikes_index,th, thUP,EEG_marker)
thresh=  mean(TEOsignal)+(-coefstd)*std(TEOsignal);

thresholdTEO=thresh*ones(size(TEOsignal));

markernumber_detection=[];
binary = TEOsignal >=thresh ;
primary_binary=binary;
detectedsamplemarker_duration=[];

% find consequtive detected RP samples
E = binary(2:end)-binary(1:end-1);
sise = size(binary);

begins = find(E==1)+1;

if binary(1) == 1
    if sise(1) > 1
        begins = [1; begins];
    elseif sise(2) > 1
        begins = [1 begins];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(begins) == 0 && binary(1) == 0
    begins = NaN;
end

ends = find(E==-1);
if binary(end) == 1
    if sise(1) > 1
        ends = [ends; length(binary)];
    elseif sise(2) > 1
        ends = [ends length(binary)];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(ends) == 0 && binary(end) == 0
    ends = NaN;
end


EEG_marker_detection1=binary;
markernumber_detection1=begins;
fprintf('primary detection by thresholding on TEO\n')
EEG_marker_detection=zeros(size(TEOsignal));
binary2=zeros(size(TEOsignal));
EEG_marker_detection2=zeros(size(TEOsignal));
counter=1;
markernumber_detection=[];

for i=1:length(begins);
    
    if begins(1,i)-round(0.1*fs)<=0
        ss=1;
    else
        ss=begins(1,i);
    end
    
    st=ends(1,i);
    
    ramp1=x_dLP(1,st)-x_dLP(1,ss);
    ramp2=abs(ramp1);
    wave1=x_dLP(1,ss:st);
    len=length(wave1);
    xmax1=find(wave1==max(wave1));
    xmin1=find(wave1==min(wave1));
    ramp3=((x_dLP(1,ss)-min(wave1)));%/(xmin1-1));
    ramp4=((x_dLP(1,st)-min(wave1)));%/(length(wave1)-xmin1));
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % fft_length=256;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % xdft = fft(Cz1(1,ss:st),fft_length);% win_length or 256 point fft
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % xdft1 = xdft(1,1:fft_length/2+1);%N/2+1
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % N=length(xdft1);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % psdx = (1/(N*N)).*abs(xdft1).^2;%Power Spectrum (spectrogram (freq,power)) Estimates Using FFT
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % psdx = 2*real(psdx);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % freq = 0: fs/fft_length: fs/2;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % % a1=find(freq>=1);b1=find(freq<=20);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     % %  [mf, sef, df] = specstats(freq(a1(1):b1(end))', psdx(a1(1):b1(end))')
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % %     th=5;
    % % % %     thUP=50;
    % % % %
 
    % check if the detected area overlaps with artefacts areas (spikes_index)
    overlapp=zeros(size(spikes_index));
    overlap2=zeros(size(spikes_index));
    
    len1=length(find((overlapp==overlap2)));
    
    if (len1==length(spikes_index))  && len>=(duration*fs) && (ramp3>=(th) && ramp4>=(th)) && ramp3<=thUP && ramp4<=thUP
        detectedsamplemarker_duration(counter,:)=[begins(1,i)  (ends(1,i)-begins(1,i))];
        markernumber_detection(1,counter)=ends(1,i); % last point of detection
        
        binary(1,begins(1,i):ends(1,i))=1;
        if begins(1,i)-round(0.6*fs)>0
            binary2(1,begins(1,i)-round(0.6*fs) :ends(1,i))=1; % movement imagery detection %for Berlin BCIC ...
            markernumber_detection2(1,counter)=begins(1,i)-round(0.6*fs);
        else
            binary2(1,begins(1,i) :ends(1,i))=1; % movement imagery detection %for Berlin BCIC ...
            markernumber_detection2(1,counter)=begins(1,i);
        end
        
        counter=counter+1;
        
    else
        binary(1,begins(1,i):ends(1,i))=0;
        
        binary2(1,begins(1,i) :ends(1,i))=0; % movement imagery detection %for Berlin BCIC ...
        
        
    end
end
EEG_marker_detection=binary;
EEG_marker_detection2=binary2;
fprintf('RP calculation done ... \n');
fprintf('removing FPs by some morphological constraints on denoised and lowpass filtered signal\n');


toc,


RMSE=sqrt(mean(((EEG_marker_detection2-EEG_marker).^2)));

end

