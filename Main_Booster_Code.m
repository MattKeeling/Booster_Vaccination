%%
clear all

Remove_New_Variant=0;

XL_output=0;
Matlab_Save=1;

Rmax=11;
Rtot=[2:11];

RepeatBoost=1;
SixMonthBoost=1;


RUN_UNTIL=datenum(2023,10,20);

load Main_Parameters.mat


for TYPE=0:1
    
    if TYPE==0
        BT='LongerLasting';
    else
        BT='Waning';
    end
    
    load MCMC_Medians.mat
    
    [V1,V2,VET1,VET2,RC,RCH, RCD,Phase, RatioPf]=Vaccination_over_Time_MTP_UK(UT,VIPfAZ,VSPfAZ,VHPfAZ,VDPfAZ);
    
    [V3,Phase]=Boosters_over_Time_MTP_UK(BUT,6,V2);
    
    VET1(isnan(VET1))=mean(VET1,'all','omitnan');
    VET2(isnan(VET2))=mean(VET2,'all','omitnan');
    
    
    %%
    % generate run_stops by region
    Comps = zeros(4,Num_Comp+1,11);
    
    tmpRS=RUN_STOP;  Z=(RUN_STOP(end-1)+7):7:(ceil(RUN_UNTIL+1-datenum(2020,1,1)));
    tmpRS(Num_Comp+[1:length(Z)]-1)=Z;
    
    RUN_STOPs = zeros(11,Num_Comp+1);
    for Region = 2:11
        RUN_STOPs(Region,1:length(tmpRS)) = tmpRS;
    end
    
    
    Keep_RUN_STOPs=RUN_STOPs(2,:);
    mx=max(RUN_STOPs,[],'all')+100;
    
    if size(V1,2)<max(RUN_STOPs,[],'all')
        ed=min([size(V3,2) size(V1,2) size(RCH,2)]);
        for r=1:11  for a=1:21
                V1(r,ed:mx,a)=V1(r,ed,a);
                V2(r,ed:mx,a)=V2(r,ed,a);
                V3(r,ed:mx,a)=V3(r,ed,a);
                RC(r,ed:mx,a)=RC(r,ed,a);
                RCH(r,ed:mx,a)=RCH(r,ed,a);
                RCD(r,ed:mx,a)=RCD(r,ed,a);
            end
        end
    end
    
    
    if SixMonthBoost==2 % Repeated Booster every 6 months
        V3(:,600+182+[1:182],:)=V3(:,600+[1:182],:);
    end
    
    % Boosters in Later Years.
    if RepeatBoost
        Yrs=(size(V3,2)-600)/365;
        for y=1:(Yrs+1)
            V3(:,600+365*y+[1:365],:)=0.7*V3(:,600+[1:365],:)/0.9;   % 70\% in following years.
        end
    end
    
    if SixMonthBoost==1 % Extra Booster at 6 months.
        V3(:,600+182+[1:182],:)=V3(:,600+[1:182],:);
    end
    
    %%
    maxtime = max(RUN_STOPs(end,:))+30; % needs to be big enough to take account of any lags
    nALL = zeros(11,maxtime,63,length(nALPHA));
    nDEATHS = zeros(11,maxtime,63,length(nALPHA));
    nHOSP_AD = zeros(11,maxtime,63,length(nALPHA));
    nHOSP_OCC = zeros(11,maxtime,63,length(nALPHA));
    nICU_AD = zeros(11,maxtime,63,length(nALPHA));
    nICU_OCC = zeros(11,maxtime,63,length(nALPHA));
    
    for QQ=1:3  % Run the three different assumptions about waning efficacy.
        
        H_STRETCH=nH_STRETCH(:,QQ);
        I_STRETCH=nI_STRETCH(:,QQ);
        
        FACTOR=nFACTOR(:,QQ);
        D_FACTOR=nD_FACTOR(:,QQ);
        H_FACTOR=nH_FACTOR(:,QQ);
        I_FACTOR=nI_FACTOR(:,QQ);
        TAU=nTAU(QQ);
        SCALING=nSCALING(:,QQ);
        COMPLIANCE=nCOMPLIANCE(:,:,QQ);
        COMPLIANCEO=nCOMPLIANCEO(:,:,QQ);
        ALL2PP=nALL2PP(:,QQ);
        
        NV_SPEED=nNV_SPEED(:,QQ);
        NV_BETA=nNV_BETA(:,QQ);
        IMPORT_LEVEL=nIMPORT_LEVEL(:,QQ);
        IMPORT_DAY=nIMPORT_DAY(:,QQ);
        
        LAG=nLAG(:,QQ);
        INC_P=nINC_P(QQ);
        ALPHA=nALPHA(QQ);
        GAMMA=nGAMMA(QQ);
        
        DH_SCALING=nDH_SCALING(:,QQ);
        
        Detection=nDETECTION(:,QQ)';
        Susceptibility=nSUSCEPTIBILITY(:,QQ)';
        gamma=nGAMMA(QQ);
        
        if QQ==1
            WaningEfficacy=[0.5 0.5 0.3 0.3 0.3 0.3 0.8 0.8];
            WaningSpeed=[1/360 1/60 1/1500 1/120];
        end
        if QQ==2
            WaningEfficacy=[0.3 0.3 0.3 0.3 0.3 0.3 0.8 0.8];
            WaningSpeed=[1/360 1/60 1/1500 1/250];
        end
        if QQ==3
            WaningEfficacy=[0.0 0.0 0.3 0.3 0.3 0.3 0.8 0.8];
            WaningSpeed=[1/360 1/60 1/1500 1/400];
        end
         
        % generate Comp by region
        clear Comps CompsO tmpO
        
        RUN_STARTs=[82*ones(11,1) RUN_STOPs(:,1:(end-1))];
        
        Christmas=0;
        
        MAX_REGION=11;
        
        NC=Num_Comp;
        for Region = Rtot
            Comps(Region,1:NC)=COMPLIANCE(1:NC,Region);
            Comps(Region,NC:size(RUN_STOPs,2))=Comps(Region,NC);
            CompsO(Region,1:NC)=COMPLIANCEO(1:NC,Region);
            CompsO(Region,NC:size(RUN_STOPs,2))=CompsO(Region,NC);
            
            DT=datenum(2021,10,31)-datenum(2020,1,1);
            DT2=KDT(QQ);
            DTO2=KDTO(QQ);
            
            for Q=find(RUN_STOPs(2,:)>DT,1,'first'):size(RUN_STOPs,2)
                Comps(Region,Q)=Comps(Region,Num_Comp)*(1-(RUN_STOPs(2,Q)-DT)/(DT2-DT));
                tmpO(Region,Q)=(CompsO(Region,Num_Comp)+Comps(Region,Num_Comp))*(1-(RUN_STOPs(2,Q)-DT)/(DTO2-DT));
                if Comps(Region,Q)<0.01 Comps(Region,Q)=0.001; end
                if tmpO(Region,Q)<0.0 tmpO(Region,Q)=0.001; end
                CompsO(Region,Q)=tmpO(Region,Q)-Comps(Region,Q);
            end
            
        end
        
        MAX_REGION=Rmax;
        
        
        TXT_STR=['Booster_Output'];
        
        parfor Region=Rtot
            %for Region=2:MAX_REGION
            if TYPE==0
                [nHospital, inHospital, nICU, inICU, nDeaths, ~, All]=CALL_ODEs_REPEATEDLY(Region, TAU, ALPHA, INC_P, SCALING(Region), FACTOR(Region), H_FACTOR(Region+[0 11 22]), I_FACTOR(Region+[0 11 22]), D_FACTOR(Region+[0 11 22]),  ...
                    H_STRETCH(Region+[0 11 22]), I_STRETCH(Region+[0 11 22]), LAG(Region), START_DATE(Region)+1, 1, Comps(Region,:), CompsO(Region,:), RUN_STOPs(Region,:), NV_BETA(Region+[0 11])', NV_SPEED(Region+[0 11])', IMPORT_DAY(Region+[0 11]),IMPORT_LEVEL(Region+[0 11]),...
                    squeeze(V1(Region,:,:)), squeeze(V2(Region,:,:)), squeeze(V3(Region,:,:)), VET1(Region,:), VET2(Region,:), Transmission_Reduction, squeeze(RC(Region,:,:)), squeeze(RCH(Region,:,:)), squeeze(RCD(Region,:,:)),VoC_Sus, General_Reduction, RatioPf,...
                    Detection, Susceptibility, gamma, WaningSpeed, WaningEfficacy, [0.1 0],@Longer_Booster_ODEs);
            else
                [nHospital, inHospital, nICU, inICU, nDeaths, ~, All]=CALL_ODEs_REPEATEDLY(Region, TAU, ALPHA, INC_P, SCALING(Region), FACTOR(Region), H_FACTOR(Region+[0 11 22]), I_FACTOR(Region+[0 11 22]), D_FACTOR(Region+[0 11 22]),  ...
                    H_STRETCH(Region+[0 11 22]), I_STRETCH(Region+[0 11 22]), LAG(Region), START_DATE(Region)+1, 1, Comps(Region,:), CompsO(Region,:), RUN_STOPs(Region,:), NV_BETA(Region+[0 11])', NV_SPEED(Region+[0 11])', IMPORT_DAY(Region+[0 11]),IMPORT_LEVEL(Region+[0 11]),...
                    squeeze(V1(Region,:,:)), squeeze(V2(Region,:,:)), squeeze(V3(Region,:,:)), VET1(Region,:), VET2(Region,:), Transmission_Reduction, squeeze(RC(Region,:,:)), squeeze(RCH(Region,:,:)), squeeze(RCD(Region,:,:)),VoC_Sus, General_Reduction, RatioPf,...
                    Detection, Susceptibility, gamma, WaningSpeed, WaningEfficacy, [0.1 0],@Waning_Booster_ODEs);
            end
            
            
            nDeaths=nDeaths.* ((1+DH_SCALING(Region)*sum(inHospital,2)./sum(Region_PP(Region,:),2))*ones(1,size(nDeaths,2)));
            
            T=1:size(nHospital,1);
            
            padding = maxtime-length(T);
            nALL_INF(Region,:,:,QQ)=[All; zeros(padding,63)];
            nDEATHS(Region,:,:,QQ)=[nDeaths; zeros(padding,63)];
            nHOSP_OCC(Region,:,:,QQ)=[inHospital; zeros(padding,63)];
            nHOSP_AD(Region,:,:,QQ) = [nHospital; zeros(padding,63)];
            nICU_OCC(Region,:,:,QQ) = [inICU; zeros(padding,63)];
            nICU_AD(Region,:,:,QQ) = [nICU; zeros(padding,63)];
            
        end
    end
    
    %% NEED TO COMBINE ALL THE VARIANTS
    a=1:21; b=21+a; c=42+a;
    nALL_INF=nALL_INF(:,:,a,:)+nALL_INF(:,:,b,:)+nALL_INF(:,:,c,:);
    nDEATHS=nDEATHS(:,:,a,:)+nDEATHS(:,:,b,:)+nDEATHS(:,:,c,:);
    nHOSP_OCC=nHOSP_OCC(:,:,a,:)+nHOSP_OCC(:,:,b,:)+nHOSP_OCC(:,:,c,:);
    nHOSP_AD=nHOSP_AD(:,:,a,:)+nHOSP_AD(:,:,b,:)+nHOSP_AD(:,:,c,:);
    nICU_OCC=nICU_OCC(:,:,a,:)+nICU_OCC(:,:,b,:)+nICU_OCC(:,:,c,:);
    nICU_AD=nICU_AD(:,:,a,:)+nICU_AD(:,:,b,:)+nICU_AD(:,:,c,:);
    
    if Matlab_Save
        save([TXT_STR '_' BT '.mat'],'RepeatBoost','SixMonthBoost','nDEATHS','nHOSP*','nALL_INF','nINC*','nDEATHS');
    end
end
