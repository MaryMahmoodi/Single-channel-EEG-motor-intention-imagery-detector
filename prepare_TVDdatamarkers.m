function [Cz,Cz0,markernumber,EEG_marker,fs]=prepare_TVDdatamarkers(EEG, markernumber,fs,landau,duration1)
    % channels: 9 11 13 (C3, Cz, C4)
    Cz=EEG(11,:);
    if nargin<5;duration1=length(Cz);end
    Cz=Cz(1:duration1*fs);
    
    w0 = 50/(fs/2);
    [num,den] = iirnotch(w0,w0/35);
    Cz = filter(num,den,Cz);
    Cz0=Cz;
    
    figure (2);clf;
    landau1=[ 2 5 8 10 13 15 17 20 25 30 35 40 45 50 55 60 65 70 75 80 90 95 100 105 110];
    N=length(Cz);
    I = speye(N);%sparse identity matrix
    D = I(2:N, :) - I(1:N-1, :);%I(2:N, :) - I(1:N-1, :);% fist order derivative
    cost=[];cost2=[];cost3=[];cost4=[];
    for i=1:length(landau1)
        x = tvd(Cz,length(Cz),landau1(i));
        % x= tvd_mm(Cz,landau1(i),Nit);
        
        if size(x,1)~=1;x=x';end
        
        cost(i) =20*log10(1/N* 0.5*sum(abs(x/norm(x)-Cz/norm(Cz)).^2) + 1/N*landau1(i)*sum(abs(D*x'/norm(x))));
        cost2(i)= 20*log10(1/N* 0.5*sum(abs(x/norm(x)-Cz/norm(Cz)).^2));
        cost3(i)=20*log10(1/N*landau1(i)*sum(abs(D*x'/norm(x))));
        cost4(i)=20*log10(1/N*sum(abs(D*x'/norm(x))));
        
    end
    xlabel(' \lambda'); ylabel('Estimation error')
    hold on; plot(landau1,cost,'b');
    hold on;plot(landau1,cost2,'b--');
    hold on;plot(landau1,cost3,'b^');
    hold on;plot(landau1,cost4,'b*');
    
    legend ('20*log_1_0(1/N F(x_T_V_D))','20*log_1_0(1/2N ||(y-x_T_V_D)||_2^2)',' 20*log_1_0(\lambda/N ||(Dx_T_V_D)||_1)','20*log_1_0(1/N ||(Dx_T_V_D)||_1)');
    
       
    EEG_marker=zeros(1,length(Cz));
    num1=find(markernumber<=round(duration1*fs)); num1=num1(end); markernumber=markernumber(1:num1);
    EEG_marker(1,markernumber)=1;
    
% fast TVD from paper: a drirect algorithm for 1D total variation denoising
    tic,
    Cz=Cz-mean(Cz);
    Cz = tvd(Cz,length(Cz),landau);
    toc,


end

