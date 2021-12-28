function [TEOsignal,timescale]=make_TEOsignal(x_dLP, fs,markernumber,usetemplatematching)
if usetemplatematching
    template=zeros(1,round(3*fs));
    monset=zeros(1,round(3*fs));monset(1,round(2*fs))=1;
    for i1=2:length(markernumber)
        template =template+ x_dLP(1,(markernumber(1,i1)-round(2*fs)):  (markernumber(1,i1)+round(1*fs))-1      )  ;
        %out or Cz1
    end
    template=template/(i1-1);
    
    corr_EEG1=xcorr(template, x_dLP);% corrcoef(template, EEG1(1:length(template)))
    corr_EEG1=corr_EEG1(1:length(x_dLP))/(max(x_dLP).^2);
    
    TEOsignal=(corr_EEG1); %envelope detection
    
    
else
    
    %  calculate teager energy for whole x
    
    
    TEOsignal=T(x_dLP);
    
end
timescale=[1: length(TEOsignal)]./fs;

end

