%lifespan calculation
load tst180v3_FIGoutputsv3.mat
load('STOREXXTide2_v3.mat')

SF_SWB=-STOREXXTide2_v3(4,:);
SF_MF=-STOREXXTide2_v3(1,:);%at mudflat bound (ie 2km from SWB)
SF_7km=-STOREXXTide2_v3(2,:);%7km from seaward bound
SF_12km=-STOREXXTide2_v3(3,:);

% max_val=5*max(abs(SF_MF));%max val for sedflux at a time cannot exceed max val at mudflat boundary by too much
% SF_7km(abs(SF_7km)>max_val)=mean(SF_7km); %set extreme outliers to mean val
% SF_12km(abs(SF_12km)>max_val)=mean(SF_12km);
% 
% t02=[0.2:0.2:3000];
% figure
% plot(t02,cumsum(SF_SWB)*2650*0.3/1000/(10^6)) %cumsum converted from m3 to kg mud to tons mud
% hold on
% plot(t02,cumsum(SF_MF)*2650*0.3/1000/(10^6))
% plot(t02,cumsum(SF_7km)*2650*0.3/1000/(10^6))
% plot(t02,cumsum(SF_12km)*2650*0.3/1000/(10^6))
% xlabel('time (yrs)')
% ylabel('cumulative net sediment flux landward (Mt)')
% legend('seaward boundary','2 km from seaward boundary','7km from seaward boundary','12km from seaward boundary')

%%
SF_MF=-STOREXXTide2_v3(1,:);
yy  = reshape(SF_MF(1:end), 100, []);     % Reshape x to a [n, m/n] matrix
SF_MF_20y = transpose(sum(yy, 1) / 100); %20yr avg %m^3/20y
load tst180v3_newgtr_v3.mat
load STORE180v3_newgtr_v3.mat
load TidetimeseriesWIND
% 
 STORE=STORE180v3_newgtr_v3;
TidetimeseriesWIND=[TidetimeseriesWIND TidetimeseriesWIND TidetimeseriesWIND]; %just for long timeseries

%load in parameters. Make sure these are consistent w/ simulation:
Trange=mean(TidetimeseriesWIND(1,:));TrangePLOT=Trange;
Ttide=mean(TidetimeseriesWIND(2,:));
Trange90_o=0.67;
Ttide90=0.52;
RSLR=1.5;
% dissipation_factor=0.2;%Trange remaining by end %0.25 for tst58, can be as low as 0.1
% %[Trange]=gettidalrange_static(Trange,x,z,dissipation_factor);
% %[Trange90]=gettidalrange_static(Trange90_o,x,z,dissipation_factor);
% [Trange]=gettidalrange_exp(Trange,x,z,dissipation_factor);
% [Trange90]=gettidalrange_exp(Trange90_o,x,z,dissipation_factor);

%fully dynamic

TrangePLOT=zeros(1100,200,150); 
Trange90PLOT=zeros(1100,200,150);
UPLOT=zeros(1100,200,150);
VEGPLOT=zeros(1100,200,150);
UPLANDPLOT=zeros(1100,200,150);
CHANPLOT=zeros(1100,200,150);
OFFSHOREPLOT=zeros(1100,200,150);
UVVRPLOT=[];
UVVR_x_kmPLOT=zeros(22,150);
Chandepth=zeros(1100,150);


for i=1:length(STORE(1,1,:))
    zbed=STORE(:,:,i);
    msl=i*20*RSLR/1000;
    z=zbed+msl;
    TrangeCONST=mean(TidetimeseriesWIND(1,:));%TrangePLOT=TrangeCONST;
    Ttide=mean(TidetimeseriesWIND(2,:));
    Trange90_o=0.67;
    Ttide90=0.52;
    Cb=0.015;%P
    Cv=0.1;%P
%     dx=20;N=2200/2*1;M=400/2*1;%/4; %y    500/2*1
%     x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
%     A=ones(N,M);
%     %sea b.c.
%     A(end,:)=2;Active=A;
    [Trange]=gettidalrange(Ttide,TrangeCONST,Cb,Cv,dx,x,zbed,A);
    [Trange90]=gettidalrange(Ttide90,Trange90_o,Cb,Cv,dx,x,zbed,A);
    TrangePLOT(:,:,i)=Trange;
    Trange90PLOT(:,:,i)=Trange;
    zsill=NaN;
    distdr=NaN;
    zntwrk=0;
    minponddepth=0.05;
    periodic=0;
    dBlo=0;
    MSL90=0;%igniring MSL anomalies for plant growth now
    %dBup=MSL90+Trange90/2+0.1; %%%new formulation
    dBup=MSL90+Trange90_o+0.1; %real high bound so we dont classify levees as unvegetated
    kro=0.01;%P

    Cb=0.015;%P
    Cv=0.1;%P

    alpha_surge=0.25;
%%

%%
%plot elevation map
% zbed=z-msl;
% 
% %Trange=0.7;
% %%plotting z
% figure
% IM=imagesc(x,y,zbed');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
% cmp=demcmap([-3 TrangePLOT/2],256); %-3 1 P.Trange/2
% colormap(cmp)
% caxis([-3 TrangePLOT/2]);
% title('Z')
%%
%find ponds
% Trange=0.7;
% zsill=Trange/2;  
% zntwrk=(Trange/2)*0.1;
% minponddepth=0.05;
%[S,AC,DIF]=findisolatedponds(z-msl,A,N,M,dx,zntwrk,zsill,NaN,minponddepth);
    [S,AC,DIF]=findisolatedponds(zbed,A,N,M,dx,zntwrk,zsill,distdr,minponddepth,Active);%AC is only used for plotting, not in the computation
%     figure
%     IM=imagesc(x,y,S');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
%     title('S')
%%
%find marsh (VEG)
%dBlo=0;
    Zlev=zbed;
    VEG=Zlev>=dBlo & Zlev<=dBup;
    VEG=VEG.*(S==0);
    VEGPLOT(:,:,i)=VEG;
    UPLAND=Zlev>dBup;
    UPLANDPLOT(:,:,i)=UPLAND;
% figure
% IM=imagesc(x,y,VEG');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
% title('VEG')
%%
%find channel
% Trangelow=0.4;%P
% kro=0.1;%P
%Water depth
%[h,ho,fTide,dtide,dsurge,dHW,wl,wlo,rangespatial]=getwaterdepth_vtr(Trange,Trangelow,msl,z,kro,0);
    Hsurge=0;hpRIV=0;tempdeltaMSL=0; %no surge, river correction, or sea level excursion
    [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV,tempdeltaMSL);
% Cb=0.02;%P
% Cv=0.1;%P
    MANN=0*A+Cb;
    MANN(VEG==1)=Cv;


%DIF is the impounded water depth, need to subtract to the imput discharge from the ponded area (THE WATER REMAINS THERE!!!)
    pondleakage=0.2;
    DIF=max(DIF,0); %you cannot impound a negative water depth!!! This happens because of the trick used to swap the cell during the pond floodin
%dtideI=Trange/2-(z-msl)-DIF;%THSI IS JUST TO PLOT
    dtide(S==1)=max(0, mean(Trange/2,'all')-(z(S==1)-msl)-max(0,DIF(S==1)-pondleakage)); %1/19/22: used mean of spatially variable Trange
%reduce the hydperperio din the connecte dpond, becuase they do not exahcneg water as much as their wwater depth. Some of that depth is as if it was made of soil. only the top layer count as moving water!
    [~,~,fTidePOND]=getwaterdepth(Trange,Hsurge,msl,z+DIF,kro,hpRIV,tempdeltaMSL);fTide(S==1)=fTidePOND(S==1);
%CALCULATE TIDAL PRSIM
    DHeff=dtide+alpha_surge*min(Hsurge,dsurge);%the water level excusrion for the total water prismmin

    [U,Uy,Ux]=tidalFLOW(A,MANN,h,ho,dHW,dx,DHeff,Ttide,periodic,kro);
    U(A==10)=0;Ux(A==10)=0;Uy(A==10)=0;
    UPLOT(:,:,i)=U;
    
%[Uwave_sea,Tp_sea,Hsea,Fetch,kwave,PWsea]=SeaWaves(h,0,hwSea_lim,Trange,wind,MASK,64,dx); %72,dx    
%now turn 
    Uthresh=0.1*max(U,[],'all'); %threshold is arbitrary: 0.2 works ok, could be less
    
    CHAN=U;
    CHAN(CHAN<Uthresh)=0;
    CHAN(CHAN>=Uthresh)=1;
%now add a 1 pixel buffer ***may need to adjust depending on grid size
    CHAN=bwmorph(CHAN,'thicken');
    CHAN(zbed>0)=0;
    CHANPLOT(:,:,i)=CHAN;
% figure
% IM=imagesc(x,y,CHAN');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
% title('CHAN')
zrsl=z-msl;
CHANz=zrsl;CHANz(CHAN==0)=NaN;
Chandepth(:,i)=min(CHANz,[],2);
%%
%now we have S, VEG, CHAN--> need to find subtidal and mudflat
    rangespatial=Trange(:,100);
%subt=zbed<(-rangespatial/2);
    MF=zbed>=(-rangespatial/2);
    MF(VEG==1)=0;
    MF(S==1)=0;
    MF(CHAN==1)=0;
    MF(UPLAND==1)=0;
%     figure
%     IM=imagesc(x,y,MF');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
%     title('MF')
    %SUBT
    SUBT=zbed<(-rangespatial/2);
    SUBT(VEG==1)=0;
    SUBT(S==1)=0;
    SUBT(CHAN==1)=0;
% figure
% IM=imagesc(x,y,SUBT');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
% title('SUBT')
    connSUBT=bwlabel(SUBT);
    regions=unique(connSUBT(A==2));regions=regions(2:end);
    boundCONN=ismember(connSUBT,regions); %boundCONN=double(reshape(boundCONN,[size(A)]));
    OFFSHORE=boundCONN; %OFFSHORE(OF=SUBT(boundCONN==1);
    SUBT(OFFSHORE==1)=0;
    OFFSHOREPLOT(:,:,i)=OFFSHORE;
%%
%check that locations aren't double classified
    check=S+VEG+CHAN+MF+SUBT+UPLAND+OFFSHORE;
    max(check,[],'all')
    min(check,[],'all')
% figure
% IM=imagesc(x,y,check');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
% title('check')
%%
%calculate UVVR
    AOI=A*NaN;AOI(end-850:end-100,:)=1; %***remove 2km mudflat and upland 5km

    UVVR=((nansum((1-VEG).*AOI,'all'))-nansum(CHAN.*AOI,'all')-nansum(OFFSHORE.*AOI,'all'))/(nansum(VEG.*AOI,'all')); %exclude channel & offshore
    UVVRPLOT(i)=UVVR;
    UVVR_x=((nansum((1-VEG).*AOI,2))-nansum(CHAN.*AOI,2)-nansum(OFFSHORE.*AOI,2))./(nansum(VEG.*AOI,2)); %UVVR (excluding channel) in x-direction

    n=50; %50 blocks per km
    UVVR_x_km = arrayfun(@(i) mean(UVVR_x(i:i+n-1)),1:n:length(UVVR_x)-n+1)'; 
    UVVR_x_kmPLOT(:,i)=UVVR_x_km;
% figure
% plot(1:22,UVVR_x_km,'-o')
%ylim([0,0.65])

    zrsl=STORE(:,:,i)*AOI;
    area=length(~isnan(zrsl));% area (m2)  **NEED TO FIX UPPER BOUND FOR UPLAND
    area=area-length()
    marsh_m=TOTMARSH(t)*1000*1000;
    RSLR=.003*(marsh_m/area)+.0015*(1-(marsh_m/area))*20; %average RSLR for veg and unveg for 20y
    Qs=SF_MF_20y(t)/area;
    Qb=Qs-RSLR;
    marshZ=

end

for t=1:150
    zrsl=STORE(251:1000,:,t);
    area=length((zrsl(zrsl<0.7)))*20*20;%non upland area (m2)  **NEED TO FIX UPPER BOUND FOR UPLAND
    marsh_m=TOTMARSH(t)*1000*1000;
    RSLR=.003*(marsh_m/area)+.0015*(1-(marsh_m/area))*20; %average RSLR for veg and unveg for 20y
    Qs=SF_MF_20y(t)/area;
    Qb=Qs-RSLR;
    marshZ=