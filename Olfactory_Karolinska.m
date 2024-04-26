%INPUT
%TS: Create variable TS (timestamps) size 15x30 from the exploratory behavior
%vizual analysis
%data: time-series from recording system

%OUTPUT
%Mean_power_HP_norm
%Mean_power_OB_norm
%Mean_Coherence
%Mean_MI_1 (theta and beta):
    %Mean_MI_HP_1
    %Mean_MI_OB_1
    %Mean_MI_HP_OB_1
    %Mean_MI_OB_HP_1
%Mean_MI_2 (theta and gamma):
    %Mean_MI_HP_2
    %Mean_MI_OB_2
    %Mean_MI_HP_OB_2
    %Mean_MI_OB_HP_2    


%

%



%*****CFC Base parameters*****
% Define phase bins
nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end
%Band analysis definition for THETA-gamma coupling
%Theta
Pf1 = 7;
Pf2 = 12;

%Beta

Af1_1 = 15;
Af2_1 = 35;

%Gamma

Af1_2 = 60;
Af2_2 = 80;

%End of CFC base parameters

%

%


%Naming hippocampus and olfactory bulb variables.
HP=data(:,1)'; %Naming variables
OB=data(:,2)'; %Naming variables

FS=6250; %Frequency sample

%Designing a bandstop filter for 50Hz noise IIR
BSfilter_IIR=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',49,'HalfPowerFrequency2',51,'SampleRate',FS);

HP_F=filtfilt(BSfilter_IIR,HP); %Filtering
OB_F=filtfilt(BSfilter_IIR,OB); %Filtering

%***PREllocating power and coherence
%PREllocating power_HP
power_HP=cell(15,15); %Cell to allocate arrays
power_HP(1:225) = {NaN((FS/2)+1,1)}; %225 =15*15

%PREllocating power_OB
power_OB=cell(15,15); %Cell to allocate arrays
power_OB(1:225) = {NaN((FS/2)+1,1)} ;%225 =15*15

%PREllocating msc_HP_OB (coherence)
msc_HP_OB=cell(15,15); %Cell to allocate arrays
msc_HP_OB(1:225) = {NaN((FS/2)+1,1)}; %225 =15*15

%Calculating POWER, COHERENCE and CFC
k=1;
for j=1:2:length(TS(1,:)) %Colums
    
    for i=1:length(TS(:,1)) %Lines
        
        if TS(i,j)>0 
            Ini=TS(i,j); Ini=(Ini-0.5)*FS; Ini=round(Ini);% Just to facilitate the script
            End=TS(i,j+1); End=(End+0.5)*FS; End=round(End);% Just to facilitate the script
            
            %Selecting time-windows
            TimeStamp_HP{i,k}=HP_F(1,(Ini:End)); % Temporally EEG variable 
            TimeStamp_OB{i,k}=OB_F(1,(Ini:End)); % Temporally EEG variable
            
            %POWER
            [Pxx,W] = pwelch(TimeStamp_HP{i,k},length(TimeStamp_HP{i,k}),(length(TimeStamp_HP{i,k})/10),FS);
            power_HP{i,k}=Pxx;
            
            [Pxx,W] = pwelch(TimeStamp_OB{i,k},length(TimeStamp_OB{i,k}),(length(TimeStamp_OB{i,k})/10),FS);
            power_OB{i,k}=Pxx;
            
            %Coherence
            [Cxy,F] = mscohere(TimeStamp_HP{i,k},TimeStamp_OB{i,k},FS/2,0,FS);
            msc_HP_OB{i,k}=Cxy;
            
            %CFC
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_HP{i,k},TimeStamp_HP{i,k},FS,Pf1,Pf2,Af1_1,Af2_1,position);
            MI_HP_1(i,k)=MI;
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_OB{i,k},TimeStamp_OB{i,k},FS,Pf1,Pf2,Af1_1,Af2_1,position);
            MI_OB_1(i,k)=MI;
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_OB{i,k},TimeStamp_HP{i,k},FS,Pf1,Pf2,Af1_1,Af2_1,position);
            MI_HP_OB_1(i,k)=MI;
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_HP{i,k},TimeStamp_OB{i,k},FS,Pf1,Pf2,Af1_1,Af2_1,position);
            MI_OB_HP_1(i,k)=MI;
            
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_HP{i,k},TimeStamp_HP{i,k},FS,Pf1,Pf2,Af1_2,Af2_2,position);
            MI_HP_2(i,k)=MI;
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_OB{i,k},TimeStamp_OB{i,k},FS,Pf1,Pf2,Af1_2,Af2_2,position);
            MI_OB_2(i,k)=MI;
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_OB{i,k},TimeStamp_HP{i,k},FS,Pf1,Pf2,Af1_2,Af2_2,position);
            MI_HP_OB_2(i,k)=MI;
            
            [MI,MeanAmp] = ModIndex_v1_DAN(TimeStamp_HP{i,k},TimeStamp_OB{i,k},FS,Pf1,Pf2,Af1_2,Af2_2,position);
            MI_OB_HP_2(i,k)=MI;            
            
        else 
            power_HP{i,k}=NaN; %Complete empty with NaN
            power_OB{i,k}=NaN;
            msc_HP_OB{i,k}=NaN;
            
            MI_HP_1(i,k)=NaN; %Complete empty with NaN
            MI_OB_1(i,k)=NaN;
            MI_HP_OB_1(i,k)=NaN;
            MI_OB_HP_1(i,k)=NaN;
            
            MI_HP_2(i,k)=NaN; %Complete empty with NaN
            MI_OB_2(i,k)=NaN;
            MI_HP_OB_2(i,k)=NaN;
            MI_OB_HP_2(i,k)=NaN;            
                    
       end
        
    end
    k=k+1;
    
end

%CFC 
%Changing zero to NaN of CFC values
MI_HP_1(MI_HP_1 == 0) = NaN;
MI_OB_1(MI_OB_1 == 0) = NaN;
MI_OB_HP_1(MI_OB_HP_1 == 0) = NaN;
MI_HP_OB_1(MI_HP_OB_1 == 0) = NaN;

MI_HP_2(MI_HP_2 == 0) = NaN;
MI_OB_2(MI_OB_2 == 0) = NaN;
MI_OB_HP_2(MI_OB_HP_2 == 0) = NaN;
MI_HP_OB_2(MI_HP_OB_2 == 0) = NaN;

%OUTPUT CFC
%Mean CFC
Mean_MI_HP_1 = nanmean (MI_HP_1);
Mean_MI_OB_1 = nanmean (MI_OB_1);
Mean_MI_OB_HP_1 = nanmean (MI_OB_HP_1);
Mean_MI_HP_OB_1 = nanmean (MI_HP_OB_1);

Mean_MI_HP_2 = nanmean (MI_HP_2);
Mean_MI_OB_2 = nanmean (MI_OB_2);
Mean_MI_OB_HP_2 = nanmean (MI_OB_HP_2);
Mean_MI_HP_OB_2 = nanmean (MI_HP_OB_2);

Mean_MI_1(1,:) = Mean_MI_HP_1;
Mean_MI_1(2,:) = Mean_MI_OB_1;
Mean_MI_1(3,:) = Mean_MI_HP_OB_1;
Mean_MI_1(4,:) = Mean_MI_OB_HP_1;

Mean_MI_2(1,:) = Mean_MI_HP_2;
Mean_MI_2(2,:) = Mean_MI_OB_2;
Mean_MI_2(3,:) = Mean_MI_HP_OB_2;
Mean_MI_2(4,:) = Mean_MI_OB_HP_2;




%Mean POWER and COHERENCE

for j=1:length(power_HP(1,:)) %Colums
    
    temp_HP=NaN((FS/2)+1,length(power_HP(1,:))); %Prellocating temporaly variable 
    temp_OB=NaN((FS/2)+1,length(power_HP(1,:))); %Prellocating temporaly variable
    temp_Coherence=NaN((FS/2)+1,length(power_HP(1,:))); %Prellocating temporaly variable

    for i=1:length(power_HP(:,1)) %Lines
        if  ~isempty (power_HP{i,j})
        temp_HP(:,i)=power_HP{i,j};
        temp_OB(:,i)=power_OB{i,j};
        temp_Coherence(:,i)=msc_HP_OB{i,j};
        end
    end
    
    Mean_power_HP(:,j)=nanmean(temp_HP');
    Mean_power_OB(:,j)=nanmean(temp_OB');
    Mean_Coherence(:,j)=nanmean(temp_Coherence');
    
    clear temp_HP temp_OB temp_Coherence
    
    
end

%Selecting frequencies below 120Hz
Mean_power_HP=Mean_power_HP(1:220,:);
Mean_power_OB=Mean_power_OB(1:220,:);
Mean_Coherence=Mean_Coherence(1:220,:);



%Normalizing POWER
for j=1:(length(Mean_power_HP(1,:))) % Coluns
    for i=1:(length(Mean_power_HP(:,1))) %j = Lines
        Mean_power_HP_norm(i,j)=Mean_power_HP(i,j)./sum(Mean_power_HP(:,j));
        Mean_power_OB_norm(i,j)=Mean_power_OB(i,j)./sum(Mean_power_OB(:,j));
    end
end



%Power and Coherence calculi for the whole period - HP and OB
%HP
[Pxx,W] = pwelch(HP_F,FS,0,FS);
All_Power_HP=Pxx;
%Normalizing
All_Power_HP_norm=All_Power_HP./sum(All_Power_HP); 
All_Power_HP_norm=All_Power_HP_norm(1:220,:); %Selecting frequencies
%OB
[Pxx,W] = pwelch(OB_F,FS,0,FS);
All_Power_OB=Pxx;
%Normalizing
All_Power_OB_norm=All_Power_OB./sum(All_Power_OB);
All_Power_OB_norm=All_Power_OB_norm(1:220,:);%Selecting frequencies
%Coherence
[Cxy,F] = mscohere(HP_F,OB_F,FS,0,FS);
All_Coherence_HP_OB=Cxy;





clear Af1 Af2 Cxy End F i Ini j k MI nbin Pf1 Pf2 position Pxx W winsize MeanAmp







%% 21_05_10 Movement analysis 
close all

%EEG analysis at long distance movements periods and speed, avoiding
%cosssing with any olfactory exploration period.


%Import CSV file from DLC analysis and name it X.

%Importar os timestamps da exploracao olfatoria e nomear como TS

%Importar as variaveis HIP e OB



%Biopac recording start correlated with the video in seconds

% inicio=3; %Naive_296691

% inicio=2; %Naive 296692
% inicio=2; %Naive 296692 PPX

% inicio=2; %Sham_305976
% inicio=2; %Sham_305976 PPX

% inicio=2; %Sham_305977
% inicio=2; %Sham_305977 PPX

% inicio=2; %Sham_305978
% inicio=3; %Sham_305978 PPX

% inicio=2; %Sham_323599
% inicio=2; %Sham_323599 PPX

% inicio=2; %Sham_323600
% inicio=3; %Sham_323600PPX

% inicio=2; %Sham_323601
% inicio=2; %Sham_323601 PPX

% inicio=3; %Sham_323602
% inicio=2; %Sham_323602 PPX

% inicio=2; %Sham_323603
% inicio=2; %Sham_323603 PPX

% inicio=3; %Lesion_322155
% inicio=2; %Lesion_322155 PPX

% inicio=2; %Lesion_322156
% inicio=3; %Lesion_322156 PPX

% inicio=3; %Lesion_322157
% inicio=2; %Lesion_322157 PPX

% inicio=2; %Lesion_322158
% inicio=3; %Lesion_322158 PPX

% inicio=3; %Lesion_322159
% inicio=4; %Lesion_322159 PPX

% inicio=3; %Lesion_324066-Dead
% inicio=3; %Lesion_324066-Dead

% inicio=2; %Lesion_324067
% inicio=2; %Lesion_324067 PPX

% inicio=2; %Lesion_324068
% inicio=2; %Lesion_324068 PPX

%  inicio=2; %Lesion_324069
%  inicio=2; %Lesion_324069 PPX


%Selecionando as janelas (segundos) nas quais o animal estah explorando o estimulo olfatorio 

for j=1:2:length(TS(1,:)) %Colums
    
    for i=1:length(TS(:,1)) %Lines
        
        if TS(i,j)>0 
            Ini=TS(i,j); Ini=round(Ini)-1;% Just to facilitate the script.
            End=TS(i,j+1); End=round(End)+1;% Just to facilitate the script.
            
            Olfactory_win{i,j}=Ini:1:End; %Janelas de exploracao.
   
        end
        
    end
       
end

Olfactory_seconds = [Olfactory_win{:}]; %Segundos nos quais o animal estah explorando.


%

%


%Calculando a movimentacao do animal baseado no CSV do DLC
%Open the csv file and import as matrix and name the imported matrix as X.

Y=X;

%Selecting X and Y based on three body parts
X=[X(:,2) X(:,5) X(:,8)]; 
Y=[Y(:,3) Y(:,6) Y(:,9)];


X=mean(X'); %Mean of all X points to decrease the variability
Y=mean(Y'); %Same with Y.


for i=1:(length(X)-1)
XD=((X(i+1))-(X(i)))^2; %Distance at X axes
YD=((Y(i+1))-(Y(i)))^2; %Distance at Y axes
D(i,1)=sqrt(XD+YD); %Distance at 2D map
end

j=1;
m=9;
for i=1:10:(length(D)-9)
DM1s(j)=mean(D(i:(i+m))); %Distancia Media, em 1s.
j=j+1;
end
DM1s=DM1s';

%Velocidade media (derivada da distancia/tempo)
vel=sqrt(diff(DM1s).^2);


%

%


%*****CFC Base parameters*****
% Define phase bins
nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end
%Band analysis definition for THETA-gamma coupling
%Theta
Pf1 = 7;
Pf2 = 12;

%Beta

Af1_1 = 15;
Af2_1 = 35;

%Gamma

Af1_2 = 60;
Af2_2 = 80;

%End of CFC base parameters

%

%


%Naming hippocampus and olfactory bulb variables.
HP=data(:,1)'; %Naming variables
OB=data(:,2)'; %Naming variables

FS=6250; %Frequency sample

%Designing a bandstop filter for 50Hz noise IIR
BSfilter_IIR=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',49,'HalfPowerFrequency2',51,'SampleRate',FS);

HP_F=filtfilt(BSfilter_IIR,HP); %Filtering
OB_F=filtfilt(BSfilter_IIR,OB); %Filtering



%DEFINICAO DO LIMIAR

threshold=5*std(DM1s); %limar baseado na distancia maxima percorrida 

% threshold=5*std(vel); %limiar baseado na velocidade maxima




clear i j
k=1;
m=0;
for i=inicio:1:length(DM1s)-100 
    
    if DM1s(i,1)>threshold
        
        if (sum(Olfactory_seconds == i ))==0 %i=segundo; logo, se o segundo referente a distancia maior que o limiar (threshold) NAO estiver presente na exploracao, deve-se continuar. 
            
        ini=(i-1)*FS;
        fim=(i+0.5)*FS;
        
        [Pxx,W]=pwelch(HP_F(1,ini:fim),length(HP_F(1,ini:fim)),length(HP_F(1,ini:fim))/10,FS);
        power_HP_MOV(k,:)=Pxx;
    
        [Pxx,W]=pwelch(OB_F(1,ini:fim),length(OB_F(1,ini:fim)),length(OB_F(1,ini:fim))/10,FS);
        power_OB_MOV(k,:)=Pxx; 
    
        [Cxy,F] = mscohere(HP_F(1,ini:fim),OB_F(1,ini:fim),FS/2,0,FS);
        msc_HP_OB_MOV(k,:)=Cxy;
        
        
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_HP_MOV_1(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_OB_MOV_1(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_HP_OB_MOV_1(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_OB_HP_MOV_1(k,:)=MI;
        
        
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_HP_MOV_2(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_OB_MOV_2(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_HP_OB_MOV_2(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_OB_HP_MOV_2(k,:)=MI;         
    
        k=k+1;
        
        else
            m=m+1; %numero de vezes que a maior movimentacao coincidiu com alguma exploracao olfatoria.
            
        end
                
    end

end

clear i j

Mean_power_HP_MOV=mean(power_HP_MOV(:,1:220));
Mean_power_OB_MOV=mean(power_OB_MOV(:,1:220));
Mean_power_HP_MOV=Mean_power_HP_MOV';
Mean_power_OB_MOV=Mean_power_OB_MOV';

Mean_power_HP_MOV_norm=Mean_power_HP_MOV./sum(Mean_power_HP_MOV);
Mean_power_HP_MOV_norm=Mean_power_HP_MOV_norm';

Mean_power_OB_MOV_norm=Mean_power_OB_MOV./sum(Mean_power_OB_MOV);
Mean_power_OB_MOV_norm=Mean_power_OB_MOV_norm';

Mean_msc_HP_OB_MOV=mean(msc_HP_OB_MOV(:,1:220));
Mean_msc_HP_OB_MOV=Mean_msc_HP_OB_MOV';

Mean_MI_MOV_1(1,1)=mean(MI_HP_MOV_1);
Mean_MI_MOV_1(2,1)=mean(MI_OB_MOV_1);
Mean_MI_MOV_1(3,1)=mean(MI_HP_OB_MOV_1);
Mean_MI_MOV_1(4,1)=mean(MI_OB_HP_MOV_1);

Mean_MI_MOV_2(1,1)=mean(MI_HP_MOV_2);
Mean_MI_MOV_2(2,1)=mean(MI_OB_MOV_2);
Mean_MI_MOV_2(3,1)=mean(MI_HP_OB_MOV_2);
Mean_MI_MOV_2(4,1)=mean(MI_OB_HP_MOV_2);






% figure









%% Antigo movement analysis, no qual foram feitas as primeiras analises de oscilacao nos periodos de maior movimento.






close all

%EEG analysis at long distance movements periods and speed.

%Import CSV file from DLC analysis and name it X.

%INPUT
%DM1s = mean distance in one second, calculated by the DLC-csv file

%
%Biopac recording start correlated with the video in seconds

% inicio=3; %Naive_296691

inicio=2; %Naive 296692
% inicio=2; %Naive 296692 PPX

% inicio=2; %Sham_305976
% inicio=2; %Sham_305976 PPX

% inicio=2; %Sham_305977
% inicio=2; %Sham_305977 PPX

% inicio=2; %Sham_305978
% inicio=3; %Sham_305978 PPX

% inicio=2; %Sham_323599
% inicio=2; %Sham_323599 PPX

% inicio=2; %Sham_323600
% inicio=3; %Sham_323600PPX

% inicio=2; %Sham_323601
% inicio=2; %Sham_323601 PPX

% inicio=3; %Sham_323602
% inicio=2; %Sham_323602 PPX

% inicio=2; %Sham_323603
% inicio=2; %Sham_323603 PPX

% inicio=3; %Lesion_322155
% inicio=2; %Lesion_322155 PPX

% inicio=2; %Lesion_322156
% inicio=3; %Lesion_322156 PPX

% inicio=3; %Lesion_322157
% inicio=2; %Lesion_322157 PPX

% inicio=2; %Lesion_322158
% inicio=3; %Lesion_322158 PPX

% inicio=3; %Lesion_322159
% inicio=4; %Lesion_322159 PPX

% inicio=3; %Lesion_324066-Dead
% inicio=3; %Lesion_324066-Dead

% inicio=2; %Lesion_324067
% inicio=2; %Lesion_324067 PPX

% inicio=2; %Lesion_324068
% inicio=2; %Lesion_324068 PPX

%  inicio=2; %Lesion_324069
%  inicio=2; %Lesion_324069 PPX



%
%
%Open the csv file and import as matrix and name the imported matrix as X.

Y=X;

%Selecting X and Y based on three body parts
X=[X(:,2) X(:,5) X(:,8)]; 
Y=[Y(:,3) Y(:,6) Y(:,9)];


X=mean(X'); %Mean of all X points to decrease the variability
Y=mean(Y'); %Same with Y.


for i=1:(length(X)-1)
XD=((X(i+1))-(X(i)))^2; %Distance at X axes
YD=((Y(i+1))-(Y(i)))^2; %Distance at Y axes
D(i,1)=sqrt(XD+YD); %Distance at 2D map
end

j=1;
m=9;
for i=1:10:(length(D)-9)
DM1s(j)=mean(D(i:(i+m))); %Distancia Media, em 1s.
j=j+1;
end
DM1s=DM1s';

%Velocidade media (derivada da distancia/tempo)
vel=sqrt(diff(DM1s).^2);


%*****CFC Base parameters*****
% Define phase bins
nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end
%Band analysis definition for THETA-gamma coupling
%Theta
Pf1 = 7;
Pf2 = 12;

%Beta

Af1_1 = 15;
Af2_1 = 35;

%Gamma

Af1_2 = 60;
Af2_2 = 80;

%End of CFC base parameters

%

%


%Naming hippocampus and olfactory bulb variables.
HP=data(:,1)'; %Naming variables
OB=data(:,2)'; %Naming variables

FS=6250; %Frequency sample

%Designing a bandstop filter for 50Hz noise IIR
BSfilter_IIR=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',49,'HalfPowerFrequency2',51,'SampleRate',FS);

HP_F=filtfilt(BSfilter_IIR,HP); %Filtering
OB_F=filtfilt(BSfilter_IIR,OB); %Filtering



%DEFINICAO DO LIMIAR

threshold=5*std(DM1s); %limar baseado na distancia maxima percorrida 

% threshold=5*std(vel); %limiar baseado na velocidade maxima




clear i j
k=1;
for i=inicio:1:length(DM1s)-10 
    
    if DM1s(i,1)>threshold
        ini=(i-1)*FS;
        fim=(i+0.5)*FS;
        
        [Pxx,W]=pwelch(HP_F(1,ini:fim),length(HP_F(1,ini:fim)),length(HP_F(1,ini:fim))/10,FS);
        power_HP_MOV(k,:)=Pxx;
    
        [Pxx,W]=pwelch(OB_F(1,ini:fim),length(OB_F(1,ini:fim)),length(OB_F(1,ini:fim))/10,FS);
        power_OB_MOV(k,:)=Pxx; 
    
        [Cxy,F] = mscohere(HP_F(1,ini:fim),OB_F(1,ini:fim),FS/2,0,FS);
        msc_HP_OB_MOV(k,:)=Cxy;
        
        
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_HP_MOV_1(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_OB_MOV_1(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_HP_OB_MOV_1(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_1,Af2_1,position);
        MI_OB_HP_MOV_1(k,:)=MI;
        
        
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_HP_MOV_2(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_OB_MOV_2(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(HP_F(1,ini:fim),OB_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_HP_OB_MOV_2(k,:)=MI;
            
        [MI,MeanAmp] = ModIndex_v1_DAN(OB_F(1,ini:fim),HP_F(1,ini:fim),FS,Pf1,Pf2,Af1_2,Af2_2,position);
        MI_OB_HP_MOV_2(k,:)=MI;         
    
        k=k+1;
        
    end

end
clear i j

Mean_power_HP_MOV=mean(power_HP_MOV(:,1:220));
Mean_power_OB_MOV=mean(power_OB_MOV(:,1:220));

Mean_power_HP_MOV_norm=Mean_power_HP_MOV./sum(Mean_power_HP_MOV);
Mean_power_HP_MOV_norm=Mean_power_HP_MOV_norm';

Mean_power_OB_MOV_norm=Mean_power_OB_MOV./sum(Mean_power_OB_MOV);
Mean_power_OB_MOV_norm=Mean_power_OB_MOV_norm';

Mean_msc_HP_OB_MOV=mean(msc_HP_OB_MOV(:,1:220));
Mean_msc_HP_OB_MOV=Mean_msc_HP_OB_MOV';

Mean_MI_MOV_1(1,1)=mean(MI_HP_MOV_1);
Mean_MI_MOV_1(2,1)=mean(MI_OB_MOV_1);
Mean_MI_MOV_1(3,1)=mean(MI_HP_OB_MOV_1);
Mean_MI_MOV_1(4,1)=mean(MI_OB_HP_MOV_1);

Mean_MI_MOV_2(1,1)=mean(MI_HP_MOV_2);
Mean_MI_MOV_2(2,1)=mean(MI_OB_MOV_2);
Mean_MI_MOV_2(3,1)=mean(MI_HP_OB_MOV_2);
Mean_MI_MOV_2(4,1)=mean(MI_OB_HP_MOV_2);






figure