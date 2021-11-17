function [V3,Phase]=Boosters_over_Time_MTP_UK(UpTake, Mnths, V2)

%%

load Main_Parameters.mat

if nargin<2
    Separation=170;   % six month post dose 2.
else
    Separation=round(365*Mnths/12)-14;
end

Delay=7;
Inf_to_Symp=5;

% Assume booster doses are constant at 1.6million per week
v(datenum(2021,9,20)+1-datenum(2020,1,1):1000)=1.6e6/7;

propBUT=UpTake;

a80=[17:21]; p80=Region_PP(1:11,a80)/sum(Region_PP(2:11,a80),'all');  %Note average for UK data !!
aW=[4:13]; pW=Region_PP(1:11,aW)/sum(Region_PP(2:11,aW),'all');  %Note average for UK workers data !!

V3=zeros(11,1000,21);

%FORECASTING FORWARDS
T=(datenum(2021,9,20)+1-datenum(2020,1,1)):(datenum(2023,7,1) +1-datenum(2020,1,1));
v(end:max(T))=v(end);


RPPall=sum(Region_PP(1:11,:),2)/sum(Region_PP(2:11,:),'all');   % RPPall = proportion of population in UK
for R=2:11
    RPP(R,1:21)=Region_PP(R,:)./sum(Region_PP(2:11,:),1);
end

for R=2:11
    
    UpTake(1:21)=min(propBUT(1:21).*squeeze(sum(V2(R,:,1:21),2))'./Region_PP(R,1:21),0.99*ones(1,21));
    
    for t=T
        
        V3(R,t,:)=0;  V2(R,t,:)=0; Phase(t)=0;
        Eligible=sum(squeeze(V2(R,1:(t-Separation),:)),1)-sum(squeeze(V3(R,1:(t-1),:)),1);

        Available=RPPall(R)*v(t);
        
        
        TbV=aW;
        Already=sum(V3(R,1:t,TbV),'all');  Required=mean(0.9*(0.7e6+1.6e6)*sum(Region_PP(R,aW),2)/sum(Region_PP(2:11,aW),'all')); % 90% in Health Care Workers
        if Already<Required & Available>0
            Given=min(Required-Already,Available/2);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=1; 
        end

        TbV=a80;
        Already=sum(V3(R,1:t,TbV),'all');  Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % over 80s & CareHomes
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=2; 
        end
        
        TbV=(75/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');  Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 75-79
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=3; 
        end
        
        TbV=(70/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');  Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 70-74
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=4; 
        end
        
        TbV=[5:14];
        Already=sum(V3(R,1:t,TbV),'all');  Required=mean(UpTake(TbV))*(2.3e6+2.2e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all'); % CEV
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=4.5; 
        end
 
        TbV=(65/5);
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 65-69
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=5; 
        end
        
             
       TbV=[5:12];
       Already=sum(V3(R,1:t,TbV),'all');   Required=mean(UpTake(TbV))*(2.3e6+2.2e6+8.5e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all'); % Health Conditions <65
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=6; 
        end            
        
        TbV=(60/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 60-64
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=7; 
        end                     
        
        TbV=(55/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 55-59
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=8; 
        end 
        
        TbV=(50/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 50-54
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=9; 
        end 
        
        TbV=[9:10];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 40-49
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=10; 
        end 
        
        TbV=[7:8];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 30-39
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=11; 
        end 
        
        TbV=[5:6];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 20-29
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=12;
        end
                                                 
        TbV=[4];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 15-19
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=13; 
        end 
                                                        
        TbV=3;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 10-14
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=14; 
        end 
                                                            
        TbV=2;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 5-9
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=15; 
        end 
                                                                
        TbV=1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 0-4
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=16; 
        end        
    end
end

% Remove any spurious values
for R=2:11
    for A=1:21
        if squeeze(sum(V3(R,:,A),2))>Region_PP(R,A)
            V3(R,:,A)=V3(R,:,A)*0.99*Region_PP(R,A)/squeeze(sum(V3(R,:,A),2));
        end   
    end
end

% Put in the delay to effect
V3(:,[1:end]+Delay,:)=V3;  V3(:,1:Delay,:)=0;

