function [sensitivity,  FPR,  accuracy,   dt,  FPperMin,  Fscore]=quantitative_analysis(Cz1,markernumber,markernumber_detection,markernumber_detection2,detectedsamplemarker_duration, fs, BCI_compet,physionet,samplemarker_durarion) 
counterTP=0; counterFP=0; counterextra=0; latency=[];
counter=1;

for i=1:length(markernumber) 
    for j=1:length(markernumber_detection)
        if BCI_compet
            ss=markernumber(1,i)-0*fs;
        elseif physionet
            ss=markernumber(1,i) -.5*fs;
            
        else
            ss=markernumber(1,i)-(round(2*fs));%2
        end
        if BCI_compet
            st=markernumber(1,i)+(round(6*fs));% cue lasts for 5 seconds
        elseif physionet
            st=markernumber(1,i)+samplemarker_durarion(i,2);% cue lasts for some seconds
            
        else
            st=markernumber(1,i)+(round(0.25*fs));% 0.6
        end
        ss1=markernumber_detection(1,j)-(detectedsamplemarker_duration(j,2)/2);%end (or mid-point) of the detected window of teager change
        if BCI_compet || physionet
            ss2= markernumber_detection2(1,j);
        else
            ss2=[];
        end
        if ~isempty (ss2)
            if (ss<=ss1  && ss1<=st)  || (ss<=ss2  && ss2<=st)% if BCI_compet ||  physionet
                
                
                counterTP=counterTP+1;
                counterextra=counterextra+1;
                
                latency(1,counter)=markernumber_detection(1,j)-markernumber(1,i);
                counter=counter+1;
                %                          j=length(markernumber_detection);
                
            end
            
        else
            
            if (ss<ss1  && ss1<st)
                
                
                counterTP=counterTP+1;
                counterextra=counterextra+1;
                
                latency(1,counter)=ss1-markernumber(1,i);
                counter=counter+1;
                %                          j=length(markernumber_detection);
                
            end
            
            
        end
        
    end
end
dt=mean(latency)/fs*1000;%ms
counterFP=length(markernumber_detection)-counterTP;
FPR=(counterFP)/(counterFP+counterTP)*100;
selectivity=(counterTP)/(counterTP+counterFP)*100;
accuracy=selectivity;

if counterTP>length(markernumber);counterTP=length(markernumber);end
sensitivity=(counterTP)/(length(markernumber))*100;
FPperMin=counterFP*60*fs/length(Cz1) ; 
Fscore=2*(sensitivity*accuracy)/(sensitivity+accuracy) /100;

[sensitivity  FPR  accuracy   dt  FPperMin  Fscore];


end

