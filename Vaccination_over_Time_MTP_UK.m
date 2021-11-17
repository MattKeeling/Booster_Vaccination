function [V1,V2,VET1,VET2,RC,RCH,RCD,Phase,RatioPfMd]=Vaccination_over_Time_MTP_UK(UpTake,nVEI,nVES, nVEH, nVED)

%%
%
%%

load Main_Parameters.mat

Separation=11*7;
Delay=10;
Inf_to_Symp=5;

% VE input = Pf1, Pf2, AZ1, AZ2
RatioAZ=1-RatioPfMd;
VET1=ones(11,1).*(nVEI(1)*RatioPfMd + nVEI(3)*RatioAZ);
VE1=ones(11,1).*(nVES(1)*RatioPfMd + nVES(3)*RatioAZ);
VEH1=ones(11,1).*(nVEH(1)*RatioPfMd + nVEH(3)*RatioAZ);
VED1=ones(11,1).*(nVED(1)*RatioPfMd + nVED(3)*RatioAZ);
VET2=ones(11,1).*(nVEI(2)*RatioPfMd + nVEI(4)*RatioAZ);
VE2=ones(11,1).*(nVES(2)*RatioPfMd + nVES(4)*RatioAZ);
VEH2=ones(11,1).*(nVEH(2)*RatioPfMd + nVEH(4)*RatioAZ);
VED2=ones(11,1).*(nVED(2)*RatioPfMd + nVED(4)*RatioAZ);

a80=[17:21]; p80=Region_PP(1:11,a80)/sum(Region_PP(2:11,a80),'all');  %Note average for UK data !!
aW=[4:13]; pW=Region_PP(1:11,aW)/sum(Region_PP(2:11,aW),'all');  %Note average for UK data !!

V1=zeros(11,1000,21);  V2=zeros(11,1000,21);

%FORECASTING FORWARDS
T=300:(datenum(2023,7,1) +1-datenum(2020,1,1));
v(end:max(T))=v(end);

RPPall=sum(Region_PP(1:11,:),2)/sum(Region_PP(2:11,:),'all');   % RPPall = proportion of population in UK
for R=2:11
    RPP(R,1:21)=Region_PP(R,:)./sum(Region_PP(2:11,:),1);
end

Keep_UpTake=UpTake;

for R=2:11
    
    for t=T
        
        if t>500
            Separation=Separation-1;   %reduce until at 8 weeks.
            Separation(Separation<8*7)=8*7;
        end
        
        V1(R,t,:)=0;
        V2(R,t,:)=sum(V1(R,1:(t-Separation),:),2)-sum(V2(R,1:(t-1),:),2);
        % Only 1 dose for 12-17 year olds, that all of age-group 3 and 3/5
        % age group 4.
        V2(R,t,1:3)=0; V2(R,t,4)=sum(0.6*V1(R,1:(t-Separation),4),2)-sum(V2(R,1:(t-1),4),2);
        
        V2(R,t,(V2(R,t,:)<0))=0;
        if sum(V2(R,t,:),3)>RPPall(R)*v(t);
            V2(R,t,:)=V2(R,t,:)*RPPall(R)*v(t)/sum(V2(R,t,:),3);
        end
        
        v1=v(t)*RPPall(R)-sum(V2(R,t,:),3);
        v1(v1<0)=0;
        
        Phase(t)=2; TbV=a80;
        
        if sum(V1(R,1:t,aW),'all')>UpTake*(0.7e6+1.6e6)*sum(Region_PP(R,aW),2)/sum(Region_PP(2:11,aW),'all') % Done HCW
            VHCW=0;
        else
            VHCW=v1/2;
            V1(R,t,aW)=pW(R,:)*(VHCW)/sum(pW(R,:));
        end
        
        Vold=v1-VHCW;
        
        while(Vold>0)
            
        sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
        if sum(V1(R,1:t,TbV),'all')>sUT  % over 80s & CareHomes
            Phase(t)=3; TbV=(75/5)+1; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
            
            if sum(V1(R,1:t,TbV),'all')>sUT % 75-79
                Phase(t)=4; TbV=(70/5)+1; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                
                if sum(V1(R,1:t,TbV),'all')>sUT % 70-74
                    Phase(t)=4.5; TbV=[5:14]; sUT=mean(UpTake(TbV))*(2.3e6+2.2e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all');
                    
                    if sum(V1(R,1:t,TbV),'all')>sUT % Extremely vulnerable
                        Phase(t)=5; TbV=(65/5)+1; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                        
                        if sum(V1(R,1:t,TbV),'all')>sUT % 65-69
                            Phase(t)=6; TbV=[5:12]; sUT=mean(UpTake(TbV))*(2.3e6+2.2e6+8.5e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all');
                            
                            if sum(V1(R,1:t,TbV),'all')>sUT % Health Conditions <65
                                Phase(t)=7; TbV=(60/5)+1;  sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                
                                if sum(V1(R,1:t,TbV),'all')>sUT % 60-64
                                    Phase(t)=8; TbV=(55/5)+1;  sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                    
                                    if sum(V1(R,1:t,TbV),'all')>sUT % 55-59
                                        Phase(t)=9; TbV=(50/5)+1;  sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                        
                                        if sum(V1(R,1:t,TbV),'all')>sUT % 50-54
                                            Phase(t)=10; TbV=[9:10]; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                            
                                            if sum(V1(R,1:t,TbV),'all')>sUT % 40-49
                                                Phase(t)=10; TbV=[7:8]; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                                
                                                if sum(V1(R,1:t,TbV),'all')>sUT % 30-39
                                                    Phase(t)=10; TbV=[5:6]; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                                    
                                                    if sum(V1(R,1:t,TbV),'all')>sUT % 20-29
                                                        Phase(t)=10; TbV=[4]; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                                        
                                                        if sum(V1(R,1:t,TbV),'all')>sUT % 15-19
                                                            Phase(t)=10; TbV=3; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                                            
                                                            if sum(V1(R,1:t,TbV),'all')>sUT % 10-14
                                                                Phase(t)=10; TbV=2; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                                                
                                                                if sum(V1(R,1:t,TbV),'all')>sUT % 5-9
                                                                    Phase(t)=10; TbV=1; sUT=sum(UpTake(TbV).*Region_PP(R,TbV),'all');
                                                                    
                                                                    if sum(V1(R,1:t,TbV),'all')>sUT % 0-5
                                                                        Phase(t)=0; TbV=aW; Vold=0; sUT=0;
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        Vtmp=min(Vold, sUT);       
        V1(R,t,TbV)=squeeze(V1(R,t,TbV)) + squeeze(Vtmp*Region_PP(R,TbV)/sum(Region_PP(R,TbV),'all'))';
        Vold=Vold-Vtmp;
        if sUT==0 Vold=0; end % once we've completed vaccination.
        end
    end
end

% Remove any spurious values
for R=2:11
    for A=1:21
        if squeeze(sum(V1(R,:,A),2))>Region_PP(R,A)
            V1(R,:,A)=V1(R,:,A)*0.99*Region_PP(R,A)/squeeze(sum(V1(R,:,A),2));
        end
        if squeeze(sum(V2(R,:,A),2))>Region_PP(R,A)
            V2(R,:,A)=V2(R,:,A)*0.99*Region_PP(R,A)/squeeze(sum(V2(R,:,A),2));
        end
    end
end

Region_PP(Region_PP==0)=1;
v1=cumsum(V1,2);
v2=cumsum(V2,2);
v1o=v1-v2;
clear W1 W2 RC

for t=(2+Delay+Separation):size(V1,2)
    RC(:,t+Inf_to_Symp,:)=(squeeze(v1o(:,t-Delay,:)).*(VE1-VET1) + squeeze(v2(:,t-Delay,:)).*(VE2-VET2)) ./ (Region_PP-squeeze(v1o(:,t-Delay,:)).*VET1-squeeze(v2(:,t-Delay,:)).*VET2);
    RCD(:,t+Inf_to_Symp,:)=(squeeze(v1o(:,t-Delay,:)).*(VED1-VET1) + squeeze(v2(:,t-Delay,:)).*(VED2-VET2)) ./ (Region_PP-squeeze(v1o(:,t-Delay,:)).*VET1-squeeze(v2(:,t-Delay,:)).*VET2);
    RCH(:,t+Inf_to_Symp,:)=(squeeze(v1o(:,t-Delay,:)).*(VEH1-VET1) + squeeze(v2(:,t-Delay,:)).*(VEH2-VET2)) ./ (Region_PP-squeeze(v1o(:,t-Delay,:)).*VET1-squeeze(v2(:,t-Delay,:)).*VET2);
    
end
RC=1-RC;
RCH=1-RCH;
RCD=1-RCD;

% Put in the delay to effect
V1(:,[1:end]+Delay,:)=V1;  V1(:,1:Delay,:)=0;
V2(:,[1:end]+Delay,:)=V2;  V2(:,1:Delay,:)=0;



