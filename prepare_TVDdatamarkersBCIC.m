function [Cz,Cz0,markernumber,EEG_marker,fs]= prepare_TVDdatamarkersBCIC(nfo, cnt, mrk, fs,landau,duration1)

num=find(strcmp(lower(nfo.clab),lower('Cz')));
num2=find(strcmp(lower(nfo.clab),lower('Fz')));

cnt= 0.1*double(cnt);
Cz=cnt(:,num)';
Fz=cnt(:,num2)';
if nargin<6; duration1=length(Cz);end

if nargin<4;duration1=length(Cz);end

Cz=-Cz;
Cz=Cz(1:duration1*fs);
Fz=-Fz;

a=find(mrk.y  ==-1); % -1 left hand & +1 foot
b=find(mrk.y  ==1); % -1 left hand & +1 foot
time=mrk.pos(a)./fs; % convert to seconds
time2=mrk.pos(b)./fs; % convert to seconds

fs2=256;
[p,q]=rat(fs2/fs);
Cz=resample(Cz,p,q);
Fz=resample(Fz,p,q);
Cz=Cz-mean(Cz);
Fz=Fz-mean(Fz);

fs=fs2;

% w0 = 50/(fs/2);
% [num,den] = iirnotch(w0,w0/35);
% Cz = filter(num,den,Cz);
% Fz = filter(num,den,Fz);

Cz0=Cz;

time=round(time.*fs);
time2=round(time2.*fs);

markernumber=[time time2];


num1=find(markernumber<=round(duration1*fs)); num1=num1(end); markernumber=markernumber(1:num1);

EEG_marker=zeros(1,length(Cz));
for i=1:length(markernumber)
    if markernumber(i)+6*fs<length(Cz)
        EEG_marker(1,markernumber(i):markernumber(i)+6*fs)=1;
    else
        EEG_marker(1,markernumber(i):length(Cz))=1;
        
    end
end

% TVD


Cz = tvd(Cz,length(Cz),landau);


Fz = tvd(Fz,length(Fz),landau);


end

