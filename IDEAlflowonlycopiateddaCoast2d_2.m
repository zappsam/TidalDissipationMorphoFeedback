clear; close all;clc
%NOTES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qs is the total lateral transport, E is the vertical erosion flux

%lateral boundary options
% 0 is closed (no flux if nothign specified. no-gradient if AW equal to 1 or -1)
% 1 is periodic

%options for lateral wave condtions
%AW

%A
%0: not in the domain
%1: a normal cell
%2: the open sea boundary
%10: the river boundary
%%%FALSE%3 and -3: no-gradient boundaries (for the lateral). 3 is left (1)  -3 is right (end)

%AP= 0 is not a pond (mostly eveything); 1 if an isolate dpond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Various initiliaztion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first line is bottom %second line is top
% %this is the good one
 cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmap.cpt', cdata, ' ');

 cdata = [-1    50  255 255     0     255    255    255
          0    255 255 255     1     0    0    0];
 cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmapCLR.cpt', cdata, ' ');

cdata = [-1    205  165 0     0     0    0    0];
dlmwrite('mycmap1.cpt', cdata, ' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%to get always the same random numbers
rng(2)

%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=struct;
P.g=9.81; %gravity [m/s2]
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quatz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.kro=0.2;%0.1; % minimum water depth [m]. 0.2 is good NEEDS TO BE SMALLER THAN hwSea_lim
P.DiffS=0.1; %coefficient for tidal dispersion [-]. DO NOT CHANGE

%Sea level rise
P.RSLR=0/1000/365;  %from mm/yr to m/day (the time unit is the day!)

%Tide
P.Ttide=12.5/24; %tidal period [day]
P.Trange=1.5;%1.5; %mean tidal Trange [m]
P.TrangeVEG=P.Trange;%tidal Trange for vegetation [m]. Generally same of tidal range

%Storm surge
P.alpha_surge=0.25;%how much the surge contributes to tidal prsim   0.1;%

%Swell waves
P.gridDIR=1; %1: waves propoagation is top to bottom;   -1: waves propoagation is bottom to top
P.Pswelldir=0.4;  %if 0.5, then is symmetric  (1 is left or right) (0 is rigth to left)
P.Pswellhighangle=0.2; %if zero, only low angle waves
P.Ho=3; %boundary swell height Hs [m]
P.Tp_swell=6;% %boundary swell period Tp [m]
P.nrefrac=4;%either 0,1,2,3,4
P.multifrequency=0;%on/off
P.wavediffraction=1;%on/off
P.Cbr=0.73;%breaking coefficient
P.Cbed=0.038;%wave bed friction
P.wavefrictionCollins=1;

%Wind for sea waves
P.wind=8;%reference wind speed [m/s]

%Edge erosion
P.aw=0.3/365; %wave edge erodability m/yr/W/m2
P.fox=0;%fraction of edge eroded material that is oxidized

%Wind waves and swell numerics
P.hwSwell_lim=0; %limiter water depth for swell
P.hwSea_lim=0.5;%0.5; %limiter water deth of sea waves %THIS IS USED TO FIND THE "EDGE"
%NEEDS TO BE LARGER THAN KO!!!!

%SSC at the sea boundary
P.co1=0/1000; % Sea boundary SSC for sand [g/l]
P.co2=10/1000; %Sea boundary SSC for mud [g/l]
P.co3=0/1000; %Sea boundary SSC for mud [g/l]

%Impose the sediment discharge input at the river mouth
P.Qmouth=5; %river discharge per unit of cell [m2/s]
P.hmouth=5; %water depth [m]
%P.Umouth=P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
P.co2mouth=0;%500/1000; %SSC of mud at the river [g/l]

%Manning coeffinent unvegeated (same for sand and mud)
P.Cb=0.02;
P.alphaSAND=8; %coefficient for bedload downslope of sand

%Sand
P.d50_1=0.25/1000;% %sand grain size [m]
P.ws1=0.02;%sand with D50=500um  0.05 %m/s
P.por1=0.4;P.rbulk1=P.rhos*(1-P.por1);

%Mud
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws2=0.2/1000;%
P.por2=0.7;P.rbulk2=P.rhos*(1-P.por2);

%Mud parameters
P.me=0.1*10^-4*24*3600;  %per day!!!
P.tcrgradeint=0;%only for the small marsh simulations
P.taucr=0.2;
P.crMUD=3.65/365;%creep coeffcinet
P.crMARSH=0.2/365;%creep coeffcinet vegetated
P.DoMUD=100;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)

%Correction for proceeses duration
P.fMFswell=0.1;%0.1;%/16;%0.1; %to scale waves and tidal transport. Waves do not occur all the time
P.fMFsea=0.3;
P.fMFriver=10/365;

%Vegetation parameters
P.dBlo=-(0.237*P.TrangeVEG-0.092);
P.dBup=P.TrangeVEG/2;%-0.2;
P.Cv=0.1;%Manning for vegetated ara
P.wsB=1/1000;%Mud Settling velocity for vegetated ara
P.taucrVEG=0.5;%Critical sheak stress for vegetated areas

%Organic accretion by vegetation
% Bpeak=2.5;
% nuGp=0.0138;%1/day%rate at which organic matter is stored. to convert to AMC, total mass of organic accumalted
% chiref=0.158;%refractory fraction
% rorg=1200;%density of organic matter * 1 - water content
% Worg=0.1;
% org=(Bpeak*365/2*nuGp)*chiref/(rorg*Worg); %in m/hour
P.Korg=8/1000/365;%5/1000/365;%P.Korg=org/365;%5/1000/365;

%ON/OFF processes
P.conces=10;
P.VEGETATION=0;%vegeation resistance and settling (DOES NOT controll organic accretion)
P.computemud=0;
P.computesand=1;
P.computeSwellwave=1;
P.computeSeaWaves=1;
P.computeEdgeErosionSea=1;
P.computeEdgeErosionSwell=0;
P.computetide=1;
P.computeriver=0;
P.riverwaterlevel=0;

P.evolvestratigraphy=0;
P.VEGstratigraphy=0;%if 0 then you put the organic into the mud. If 1 then you calculate the organic as a sediment per se (advection, divergence,etc)

P.periodic=1;

P.imposeseaboundarydepthmorphoALL=0;

%Ebb-flood momentum correction
P.ebbfloodcorrection=1;
P.facMOM=0.4;%0.4;%0.5;% ALMSOT THE SAME

%Curvature flow modifications
P.curvaturecorrection=0;
P.alphacurvaturesmooth=100;
P.curvfactor=40;

P.calculateponddynamics=0;

%Stratigraphy
P.conces=10;%how much to extra erode, a parameter
P.nlyr=20; %max number of layers
P.dlyr=1; %thickenss of layers
P.tlyrU=3; %max depth to add layer %must be larger than dlyr
P.tlyrD=0.5; %min depth merge layers %mus be larger than dlyr
P.tcko=10;%tickness of bed layer
P.levo=15;%intial level occupied
P.YUi=1000;%initial thickess of active layer
P.initialfU=1;%initial composition of the active layer
P.initialf=1;%initial composition of all the layers

%Global numerical limiters
limitdeltaz=5;
limitmaxup=2;

%Time parameters
tmax=10000;%1000;%in days
tINT=1;%plot interval time step
dtO=2*365;%the reference time step
time=[0:tmax-1]*dtO/365;

%Time series
%scale in m. Shape the higher tha fatter, the smaller the more peaked
surge=1*wblrnd(0.3,0.6,tmax,1); surge(surge<1)=0; surge(surge>10)=10; %surge(1:10)=0;surge(surge<1)=0;  SCALE/SHAPE orginal%surge=wblrnd(0.1,0.5,tmax,1);
%surge=1*wblrnd(0.5,0.8,tmax,1);  surge(surge>5)=5; %surge(1:10)=0;surge(surge<1)=0;  SCALE/SHAPE orginal%surge=wblrnd(0.1,0.5,tmax,1);



makevideo=0;





%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,S]=initializegeometry_3sedimentsbarrier(P);
%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,S]=initializegeometry_Delft3DShamim(P);
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,S,'basin');
%load basin
%[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);

%Store value for mass balance check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumY1IN=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1IN=sum(sumY1IN(A==1));
sumY2IN=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2IN=sum(sumY2IN(A==1));
sumY3IN=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3IN=sum(sumY3IN(A==1));
FLX1=zeros(4,1);FLX2=zeros(4,1);FLX3=zeros(4,1);KBTOT=0;Y2OX=0;
FQsW_L=0;FQsW_R=0;
pondloss=0;

IO.S=S;
IO.Y1=Y1;IO.Y2=Y2;IO.Y3=Y3;
IO.flyr1=flyr1;IO.flyr2=flyr2;IO.flyr3=flyr3;
IO.flyrb1=flyrb1;IO.flyrb2=flyrb2;IO.flyrb3=flyrb3;
IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
fIO.pondloss=pondloss;
fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
fIO.FQsW_R=FQsW_R;fIO.FQsW_L=FQsW_L;

z=zb-(Yb+plyr*P.dlyr)-(Y1+Y2+Y3);
Y=Y1+Y2+Y3;
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
flyrU1=max(0,Y1)./Ytot;flyrU1(Ytot==0)=1;
flyrU2=max(0,Y2)./Ytot;flyrU2(Ytot==0)=1;
flyrU3=max(0,Y3)./Ytot;flyrU3(Ytot==0)=1;
zbedo=z+msl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]) %set(gca,'Color','k')
if makevideo==1;v=VideoWriter('BasinR3_P.RSLR=3','Motion JPEG AVI');open(v);end %_fillSeamudco20RSL0 _fillSeamudco20RSL2
s=0;step=0;tic;
for t=1:tmax;  
       
%swell wave direction
randdir=rand(1);if randdir>P.Pswelldir;dirsign=1;else;dirsign=-1;end %dirsign=sign((mod(t,2)-0.5));
rndhl=rand(1);if rndhl>P.Pswellhighangle;angleswell=dirsign*(rand(1)*45);else;angleswell=dirsign*min(80,(rand(1)*45+45));end
%angleswell=0;%

%SeaWave direction
angleWIND=rand(1)*360; %every time step a random direction

if t==1;dto=0.00001;else;dto=dtO;end   
dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=limitdeltaz+1;maxup=limitmaxup+1;
        while maxdeltaz>limitdeltaz | maxup>limitmaxup
        if firstattemp==1;else;dt=dt/2*min(limitdeltaz/maxdeltaz,limitmaxup/maxup);end;firstattemp=0;
        %if t<=2;dt=min(0.1*365,dt);end
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,P,dx,dt,zb,IO,fIO,surge(t),t,angleswell,angleWIND);
        step=step+1;
        end

    %the partial updating step was succefull
    IO=IOtemp;
    fIO=fIOtemp;
    dti=dti+dt;%how much you moved forward
    dt=min(dt*2,max(0,dto-dti));%the remaining time in the time step
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(t,tINT)==0;s=s+1;  
    
%read the variables
names = fieldnames(IO);
for i=1:length(names);eval([names{i} '=IO.' names{i} ';' ]);end

%read the fluxes
names = fieldnames(fIO);
for i=1:length(names);eval([names{i} '=fIO.' names{i} ';' ]);end
    
%read the plot
names = fieldnames(PLT);
for i=1:length(names);eval([names{i} '=PLT.' names{i} ';' ]);end


z=zb-(Yb+plyr*P.dlyr)-(Y1+Y2+Y3);
Y=Y1+Y2+Y3;
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
flyrU1=max(0,Y1)./Ytot;flyrU1(Ytot==0)=1;
flyrU2=max(0,Y2)./Ytot;flyrU2(Ytot==0)=0;
flyrU3=max(0,Y3)./Ytot;flyrU3(Ytot==0)=0;
zbed=z+msl;





ax1 = subplot(2,3,1);%s+1 %set(IM,'alphadata',~(A==0));
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(y,x,-zbed);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-10 3],256); %-3 1 P.Trange/2
colormap(ax1,cmp)
%colormap('jet')
caxis([-10 3]);
%colorbar('hori') 


%ax1 = subplot(2,2,2); 
% plot(x,-zbed(:,1))
% ylim([-20 5])


% figure
% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Hs);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 5]);
%axis([3 7 6 15])
%colorbar
% %  
% figure
% ax1 = subplot(2,2,2); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Ux*pi/2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
%axis([3 7 6 15])
% 
% ax1 = subplot(2,2,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,-VfloodX*pi/2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
% %axis([3 7 6 15])
% % 
% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Vebb*pi/2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
%axis([3 7 6 15])


% ax1 = subplot(2,2,2); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,VebbX);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
% %axis([3 7 6 15])*pi/2
% 
% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,-VfloodX);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
% %axis([3 7 6 15])
% 
% 
% ax1 = subplot(2,2,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Ux);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
% %axis([3 7 6 15])

% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,U*pi/2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
% %axis([3 7 6 15])
%  


% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,-VfloodX*pi/2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 3/1]);
% %axis([3 7 6 15])
% %  

% ax1 = subplot(1,4,1); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Ux);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([-2 2]);



% ax1 = subplot(2,2,2); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,U);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 1]);


%  ax1 = subplot(1,3,2); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Hs);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% %cmp=demcmap([-3 P.Trange/2],256); %-3 1
% %colormap(ax1,cmp)
% %caxis([-3 P.Trange/2]);
% %colorbar('hori') 
%  %colorbar
%  
 
%    ax2 = subplot(2,2,3);
%   IM=imagesc(y,x,UR);set(IM,'alphadata',~(A==0));axis equal;set(gca,'YDir','normal');%colormap('jet');%colorbar('hori') 
%   caxis([0 0.3]);colormap('jet')
%  colormap(ax2,flipud(hot))
 

 %max(1000*SSC(:))
%   ax2 = subplot(2,2,3);
%   IM=imagesc(y,x,1000*SSC);set(IM,'alphadata',~(A==0));axis equal;set(gca,'YDir','normal');%colormap('jet');%colorbar('hori') 
%   caxis([0 10]);colormap('jet')
%  colormap(ax2,flipud(hot))
%  %colorbar
%  
%  
% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Vebb);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 1]);
 

%   ax2 = subplot(2,3,2);
%  IM=imagesc(y,x,h);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% %cptcmap('mycmap1', 'mapping', 'direct'); 
% colormap('jet')
%  caxis([0 10])

% ax3 = subplot(3,1,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 =4]);%colorbar('hori') 
% IM=imagesc(y,x,h);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-20 3]);
% colormap(ax3,cmp)
% caxis([-20 3]);
% %colorbar('hori') 

%  ax3 = subplot(2,3,3);
%  IM=imagesc(y,x,hriver);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% %cptcmap('mycmap1', 'mapping', 'direct'); 
% colormap('jet')
%  caxis([0 3])
% 
%  ax3 = subplot(2,2,4);
%  IM=imagesc(y,x,1-flyrU1);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% cptcmap('mycmap1', 'mapping', 'direct'); 
%  caxis([0 1])
% % % 

%S(AC==1)=2;
% ax3 = subplot(2,3,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,U);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% caxis([0 0.5]);%colormap('jet')
%  %colorbar
 
%  ax4 = subplot(2,3,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,AC);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% caxis([0 1]);%colormap('jet')
%  %colorbar

% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,PLT.Ux);axis equal;set(IM,'alphadata',~(A==0));%set(gca,'YDir','normal');%colormap('jet');
% 
% ax1 = subplot(2,2,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,PLT.Uy);axis equal;set(IM,'alphadata',~(A==0));%set(gca,'YDir','normal');%colormap('jet');



%  ax3 = subplot(2,2,4);
%  IM=imagesc(y,x,Y2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
%  colormap(ax3,flipud(gray))
%  caxis([0 1])
 
%  ax1 = subplot(2,2,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Hs);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% caxis([0 4])
% colorbar
%  




% subplot(3,1,3);%MUD vs total
% %transx_y=1;Mi=20; %cross section along direction x(1) or y(2)
% transx_y=2;Mi=N; %cross section along direction x(1) or y(2)
% [si,zi,Ai]=getstrat2plot(-zb,flyrU2,flyrb2,flyr2,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
% hold off;pcolorCENTER(si-dx/2/1000,zi-msl,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
% cptcmap('mycmap', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.001]);shading flat
% %colorbar  
% hold on;plot(si,si*0+msl*0+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+0*msl+P.Trange/2,'--k',si,si*0+0*msl-P.Trange/2,'--k')
% ylim([-5 2]);%caxis([0 1]) si,-zbed(Mi,:)*NaN,'.-b',
% 

% subplot(4,1,3);%MUD vs total
% transx_y=1;Mi=20; %cross section along direction x(1) or y(2)
% [si,zi,Ai]=getstrat2plot(-zb,flyrU1+flyrU3,flyrb1+flyrb3,flyr1+flyr3,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
% hold off;pcolorCENTER(si-dx/2/1000,zi,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
% cptcmap('mycmap', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.01]);shading flat
% %colorbar
% hold on;plot(si,si*0+msl+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+msl+P.Trange/2,'--k',si,si*0+msl-P.Trange/2,'--k')
% ylim([-15 3]);%caxis([0 1])

% subplot(4,1,4);%ORGANIC vs total
% transx_y=1;Mi=20; %cross section along direction x(1) or y(2)
% [si,zi,Ai]=getstrat2plot(-zb,flyrU1+flyrU2,flyrb1+flyrb2,flyr1+flyr2,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
% hold off;pcolorCENTER(si-dx/2/1000,zi,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
% cptcmap('mycmapCLR', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.01]);shading flat
% hold on;plot(si,si*0+msl+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+msl+P.Trange/2,'--k',si,si*0+msl-P.Trange/2,'--k')
% ylim([-15 3]);%caxis([0 1])


title(strcat(num2str(time(t)),' years ',num2str(step)))
%title(strcat(num2str(time(t)),' years '))



%TERM1; Qouthriver: if postive it enters
%TERM2; Qseatide: if postive it exits
%TERM3; Qseariver: if postive it exits
%TERM4 Qmouth tide. THIS IS IMPOSED ZERO BY setting D=0 at the mouth in sedtran

%SAND
QmouthRiver=FLX1(1);QseaTide=FLX1(2);QseaRiver=FLX1(3);QmouthTide=FLX1(4);
sumFLUX1=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%MUD
QmouthRiver=FLX2(1);QseaTide=FLX2(2);QseaRiver=FLX2(3);QmouthTide=FLX2(4);
sumFLUX2=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%ORG
QmouthRiver=FLX3(1);QseaTide=FLX3(2);QseaRiver=FLX3(3);QmouthTide=FLX3(4);
sumFLUX3=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;

sumY1=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1=sum(sumY1(A==1));
sumY2=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2=sum(sumY2(A==1));
sumY3=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3=sum(sumY3(A==1));

%NOTE: Thsi is the equivalent volumetric flux, not the mass flux
checksum=[[(sumY1IN-sumY1)+sumFLUX1/dx^2]+fIO.FQsW_L+fIO.FQsW_R  [(sumY2IN-sumY2)+sumFLUX2/dx^2]+pondloss+KBTOT  [(sumY3IN-sumY3)-Y2OX+sumFLUX3/dx^2]]
%[fIO.FQsW_R fIO.FQsW_R]
%if abs(checksum(1))>1;pause;end
%RfIO.FQsW_L
%fIO.FQsW_R
%+fIO.FQsW_L+fIO.FQsW_R

if makevideo==1;V=getframe(figure(1));writeVideo(v,V);end
pause(0.1)
end

end

if makevideo==1;close(v);end %UNCOMMENT THIS TO CREATE A VIDEO



