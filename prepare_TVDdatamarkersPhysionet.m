function [Cz,Cz0,markernumber,EEG_marker,samplemarker_durarion,fs]=prepare_TVDdatamarkersPhysionet (Cz, allRPsamplemarker_durarion, subjectnumber, landau, fs)
    samplemarker_durarion=allRPsamplemarker_durarion{1,subjectnumber};
    Cz=-Cz;
    
    Fs=fs;
    N = length(Cz);
    time=(1:length(Cz))./fs;
    
    samplemarker_durarion(:,1)=samplemarker_durarion(:,1)./fs;
    % %  convert to second
    fs2=256;
    [p,q]=rat(fs2/fs);
    Cz=resample(Cz,p,q);
    
    fs=fs2;Fs=fs2;
    
    w0 = 50/(fs/2);
    [num,den] = iirnotch(w0,w0/35);
    Cz = filter(num,den,Cz);
    Cz0=Cz;
    
    
    timescale=(1:length(Cz))./fs;
    samplemarker_durarion=round(samplemarker_durarion.*fs); %sample
    markernumber= samplemarker_durarion(:,1)';
    
    EEG_marker=zeros(1,length(Cz));
    for i=1:size(samplemarker_durarion,1)
        ss=samplemarker_durarion(i,1);
        st=ss+samplemarker_durarion(i,2);
        EEG_marker(1,ss:st)=1  ;
    end
    
    
    
    % TVD
       
    Cz=Cz-mean(Cz);
    Cz = tvd(Cz,length(Cz),landau);

end

