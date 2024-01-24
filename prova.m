clear; close all;clc
%things chnaged March 27 2019
%fetch smoothing SeaWave
%fetch trhsold distance SeaWave
%fetch throdlld method  SeaWave
%ko for waves    in total sedimenersoionMUD
%fw 0.015 instead of variable    in total sedimenersoionMUD
%EXTRA DOWNLSOPE FOR MUD PROPORTiONAL TO Qs2
%wind intensity and freqnecuy fWsea


%NOTES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Various initiliaztion for plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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





%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=struct;
%STORE_t200_r07_ALL=zeros(500,300,16);
% 
%gR=[1 1 1 1 4 4 4 4 7 7 7 7 10 10 10 10];
% gC=[50 25 10 5 50 25 10 5 50 25 10 5 50 25 10 5];

gR=[2.5 2.5 2.5 2.5 5 5 5 5 7.5 7.5 7.5 7.5 10 10 10 10];
gC=[50 25 10 5 50 25 10 5 50 25 10 5 50 25 10 5]*1.5;
%gC=[60 30 15 5 60 30 15 5 60 30 15 5 60 30 15 5];
figure('units','normalized','outerposition',[0 0 1 1]) %set(gca,'Color','k')
for contatore=1
    contatore

%to get always the same random numbers
rng(2)

P.g=9.81; %gravity [m/s2]
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quatz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.kro=0.1;%0.05;%0.1;%0.2;%1;%0.2;%0.02;%02;%0.5;%0.2;%0.2;%0.1; % minimum water depth [m]. 0.2 is good NEEDS TO BE SMALLER THAN hwSea_lim
P.DiffSsand=1; %coefficient for tidal dispersion [-]. 1 DO NOT CHANGE
P.DiffSmud=1;%0.5;%0.2; %coefficient for tidal dispersion [-]. 0.1 DO NOT CHANGE
P.DoMUD=1;%10;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)
%P.DiffSmud=0.2; %coefficient for tidal dispersion [-]. 0.1 DO NOT CHANGE
%P.DoMUD=1;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)

%Sea level rise
P.RSLR=gR(contatore)/1000/365;  %from mm/yr to m/day (the time unit is the day!)
%P.RSLR=1/1000/365;  %from mm/yr to m/day (the time unit is the day!)

%Tide
P.Ttide=24/24;%12.5/24; %tidal period [day]
P.Trange=0.7;%1.4;%0.5;%0.8;%1.4; %mean tidal Trange [m]
P.TrangeVEG=P.Trange;%tidal Trange for vegetation [m]. Generally same of tidal range

%Storm surge
P.alpha_surge=0.25;%1;%2*0.25;%0.25;%0.25;%how much the surge contributes to tidal prsim   0.1;%

%Swell waves
P.gridDIR=1; %1: waves propoagation is top to bottom;   -1: waves propoagation is bottom to top
P.Pswelldir=0.5;%0.5;  %if 0.5, then is symmetric  (1 is left or right) (0 is rigth to left)
P.Pswellhighangle=0.0; %if zero, only low angle waves
P.Ho=3;%1.7;%2; %boundary swell height Hs [m]
P.Tp_swell=8;%8;%6;% %boundary swell period Tp [m]
P.nrefrac=4;%either 0,1,2,3,4  Wave refraction. If zero there is no wave refraction
P.multifrequency=0;%on/off
P.wavediffraction=1;%on/off
P.Cbr=0.73;%73;
P.Cbed=0.038;%NaN;%wave bed friction if you use Jonswap (0.067 or 0.038)
P.wavefrictionCollins=0;

%Wind for sea waves
P.wind=7;%6.5;%6.9;%reference wind speed [m/s]

%Edge erosion
P.aw=0.3/365; %wave edge erodability m/yr/W/m2
P.maxedgeheight=1;%2;
P.fox=1;%0.5;%0.2;%fraction of edge eroded material that is oxidized. 20% is oxydices (=lost)

%Wind waves and swell numerics
P.hwSwelltransport_lim=1;
P.hwSwell_lim=0.2; %limiter water depth for swell: swells are imposed zero below this limit
P.hwSea_lim=0.2;%0.5;%0.5; %limiter water deth of sea waves %THIS IS USED TO FIND THE "EDGE"%NEEDS TO BE LARGER THAN KO!!!!

%SSC at the sea boundary
P.co1=0/1000; % Sea boundary SSC for sand [g/l]
P.co2=gC(contatore)/1000;%40/1000; %Sea boundary SSC for mud [g/l]
%P.co2=50/1000;%40/1000; %Sea boundary SSC for mud [g/l]
P.co3=0/1000; %Sea boundary SSC for mud [g/l]

%Impose the sediment discharge input at the river mouth
P.Qmouth=5; %river discharge per unit of cell [m2/s]
P.hmouth=5; %water depth [m]
%P.Umouth=P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
P.co2mouth=0;%500/1000; %SSC of mud at the river [g/l]

%Manning coeffinent unvegeated (same for sand and mud)
P.Cb=0.02;

%Sand
P.d50_1=0.25/1000/1;% %sand grain size [m]
P.ws1=0.02/1;%sand with D50=500um  0.05 %m/s
P.por1=0.4;P.rbulk1=P.rhos*(1-P.por1);

%Sand Parameters: Downslope paramters for sand and mud (proportional to the sediment transport Qs!!!!)
P.alphaSAND=15; %coefficient for bedload downslope of sand. Calibrated with JMSE 2018, do not change
P.hlimC=1;%limit to apply to downslope of sand (and maybe also mud) Calibrated with JMSE 2018, do not change

%Mud
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws2=0.2/1000;%0.2/1000;% m/s
P.por2=0.7;P.rbulk2=P.rhos*(1-P.por2);%P.por2=0.7
%P.por2=0.8;P.rbulk2=P.rhos*(1-P.por2);

%Mud parameters
P.me=0.1*10^-4*24*3600;  %per day!!!
P.taucr=0.2;
P.tcrgradeint=0;% Pa/m
P.leveltauincrease=P.TrangeVEG/2;%1;
P.crMARSH=0.2/365;%creep coefficient vegetated
P.crMUD=3.65/365;%3.65/365;%3.65/365;%creep coeffcinet
P.alphaMUD=0.5;%0.2;%0.5;%*0.1;%0.1;%1; %coefficient for bedload downslope of mud. added April 2019. Similar to P.alphaSAND

%Correction for proceeses duration (to scale waves and tidal transport, because waves do not occur all the time!)
P.fMFswell=1; %use 1 if you use the equvilanet wave height
P.fMFsea=1; %use 1 if you use the equvilanet wind speed
P.fMFriver=10/365;

%Vegetation parameters
P.dBlo=0;%-0.2;%-0.2;%0;%-(0.237*P.TrangeVEG-0.092);
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
P.AccreteOrganic=1;
P.Korg=6/1000/365;%8/1000/365;%P.Korg=org/365;%5/1000/365;

%ON/OFF processes
P.VEGETATION=1;%vegeation resistance and settling (DOES NOT controll organic accretion)
P.depthlimiterflow_withVEG=0;%P.depthlimiterflow_withVEGVALUE=?????
P.computemud=1;
P.computesand=0;
P.computeSwellwave=0;
P.computeSeaWaves=1;
P.computeEdgeErosionSwell=0;
P.computeEdgeErosionSea=1;
P.compute_currentbankerosion=0;
P.computetide=1;
P.computeriver=0;

%Correction for second-order river dynamics
P.riverwaterlevel=0;
P.rivermomemntumcorrection=0;

%Various boundary condtions
P.periodic=0;
P.imposeseaboundarydepthmorphoALL=1; %to use when a channel mouth is at a boundary

%Ebb-flood momentum correction
P.ebbfloodcorrection=1;
P.residualcurrents=0;

%Curvature flow modifications
P.curvaturecorrection=0;
%P.alphacurvaturesmooth=100;
%P.curvfactor=40;

%Pond dynamics
P.calculateponddynamics=0;
    P.Epondform=4*10^-4*0.01*10;%probabiliy of new pond formation per year (area/area)
    %P.Epondform=4*10^-4*0.1/1000;%probabiliy of new pond formation per year (area/area)
    %P.zpondcr=0.01;%P.Trange/4;%base of new pond formation with respect to MSL
    %P.maxdpond=999;%0.5;%maximum depth scour of a new pond
    P.zpondcr=-0.2;%P.Trange/4;%base of new pond formation with respect to MSL
    P.minponddepth=0.1;%1; %minimum depth to define a pond after you identified the "lakes"
    P.maxdpond=max(0.2,max(P.minponddepth*1.1,0.15*P.Trange));%0.5;%0.5;%maximum depth scour of a new pond
    %%
    P.zntwrk=(P.Trange/2)*0.2;%P.Trange/2*0.9;%P.Trange/2-0.3;%0.3; %depth above msl that defines the channel network.  the smaller the harder to drain!
    P.distdr=NaN;%Clealry if nan is nto used4; %m, distance of extra drainage
    %%
    P.aPEXP=0.015;%isolated pond expansion rate m/yr
    P.ponddeeprate=0.002;%m/yr

%Stratigraphy
P.evolvestratigraphy=0;
P.VEGstratigraphy=0;%if 0 then you put the organic into the mud. If 1 then you calculate the organic as a sediment per se (advection, divergence,etc)
P.VEGonsand=0;

%Stratigraphy parameters
P.conces=10;%how much to extra erode, a parameter
P.nlyr=20; %max number of layers
P.dlyr=1; %thickenss of layers
P.tlyrU=3; %max depth to add layer %must be larger than dlyr
P.tlyrD=0.5; %min depth merge layers %mus be larger than dlyr
P.tcko=10;%tickness of bed layer
P.levo=15;%intial level occupied
P.YUi=1000;%initial thickess of active layer
P.initialfU=0;%initial composition of the active layer
P.initialf=0;%initial composition of all the layers

P.reducefractionsediment=1;%this should be 1 unless you to strange stuff. ADDED JUNE 2019@@@@@@@@@@@

%Global numerical limiters
% limitdeltaz=10/4;
% limitmaxup=5;
limitdeltaz=5/2;%10;%/4;
limitmaxup=2/2;%0.5*4;%5/4;%0.5;%

%Time parameters
tmax=2400;%10000;%2000;%400;%2*199;%1000;%149;%1000;%149;%2000;%3000;%50;%42/2-1;%2000;%1250;% number of time steps
tINT=1;%how many time steps you want to do the plot (if 1 you plot every time step). Does not affect the computation
%dtO=2*365;%the reference time step [days]. The reference unit of time is days!!!
%time=[0:tmax-1]*dtO/365; %converted to years, just to plot. Does not affect the computation




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time series
numberserie=10000;%2000;%2000;%10000;%if you change this you will change the actual values in the time series, rememebr!
numberevents=numberserie/2;
lag=exprnd(1.7*365,numberevents,1);lag(lag<1)=1;
duration=exprnd(0.01*365,numberevents,1);lag(lag<0.01)=0.01;
dtOserie=ones(numberevents*2,1);
dtOserie(1:2:end)=lag;
dtOserie(2:2:end)=duration;

%forger it, let's just do it constant
dtOserie=dtOserie*0+1*365;

time=cumsum(dtOserie)/365; %converted to years, just to plot. Does not affect the computation
time=[time(2:end);time(end)];

He=exprnd(0.5,numberevents,1);
% Hoseries=ones(numberevents*2,1);
% Hoseries(1:2:end)=2;%1;
% Hoseries(2:2:end)=4+He;

surge=ones(numberevents*2,1);
surge(1:2:end)=0;
surge(2:2:end)=0.5+He*0.6;%+(-0.5+rand(numberevents,1));
surge=surge*0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Swell wave direction
angleSWELLserie=NaN*ones(numberserie,1);
    randdir=rand(numberserie,1);
    dirsign=ones(numberserie,1);dirsign(randdir<=P.Pswelldir)=-1;
%rndhl=rand(numberserie,1);
%    a=find(rndhl>P.Pswellhighangle);angleSWELLserie(a)=dirsign(a).*(rand(length(a),1)*45);
%    a=find(rndhl<=P.Pswellhighangle);angleSWELLserie(a)=dirsign(a).*(45+rand(length(a),1)*45);
angleSWELLserie=dirsign.*(rand(numberserie,1)*45/2);

%SeaWave direction
angleWINDserie=rand(numberserie,1)*360; %every time step a random direction
%angleWINDserie=180+(2-rand(numberserie,1))*90;%*360; %every time step a random direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sediments_basindx50(P);
[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sediments_Barataria(P);
%[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry_3sediments_basindx50Georgia(P);
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,'r07TOPCo50w7Do1Dmud02kro1dmNOPONDS');
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,'bombar3NOPOND');
%load BasinNEWbig_surge05
%load BasinNEWbigfarsurge05mcontwind75
%load BwaverenewH2msicuro
%load nuovoR1
%load March_R3
%load wavetranslim1_anglespread45_R1
%load BasinRestoration
%Y1(215:218,100:130)=Y1(215:218,100:130)+5;
%Y1(215+4:218+4,100:130)=Y1(215+4:218+4,100:130)-5;

%load bombar07%NOPOND
%load gino3m
%load r07Dbase1R25long;%bomba
%load base_pozzex10
%load mudflta%NOPONDS
%load r07TOPCo50w7Do1Dmud02kro1dm
%load r07TOPCo50w7Do1Dmud02kro1dm%NOPONDS
%load r07TOPCo50w7Do1Dmud02
%load r07TOPCo60w7Do1
%load r07TOPCo50w7
%load r07TOPCo50w7NOPOND
%load r3TOPnoponds
%load r07TOP
%load r3TOP
%load r07LARGE
%load r3LARGE
%load r3LONG
%load r07
%load X
%load GINOCo25
%load GGGr3fox2
%load GGGr1
%load r3_D10_R05pondsdx10
%load r3_Do1_R05pondsdx10
%load CLEAN_R3
% Y2(310:360,100:140)=Y2(310:360,100:140)+1;
% Y2(310:460,1:50)=Y2(310:460,1:50)+1;

%Y2(610:660,340:410)=zb(610:660,340:410)+1-23+0.2;%1002.3;
%Y2(150+610:150+660,100:240)=Y2(150+610:150+660,100:240)-1.5;

%Y2(100+310:100+360,10:140)=Y2(100+310:100+360,10:140)-3;

%Y2(310:460,1:50)=Y2(310:460,1:50)+1;
%Y1(230+4:233+4,60:100)=Y1(230+4:233+4,60:100)-5;


%load Gang22R1
%load basin_noangleokbisBBBfaclong1bis
%load Basin
% Initck=2;
% zb=zb-Y1+Initck;
% Y1=Y1*0+Initck;
%load BhlimVEG___MUD
%zb(end-5:end,:)=zb(end-5:end,:)+50;

%Store value for mass balance check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumY1IN=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1IN=sum(sumY1IN(A==1));
sumY2IN=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2IN=sum(sumY2IN(A==1));
sumY3IN=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3IN=sum(sumY3IN(A==1));
FLX1=zeros(4,1);FLX2=zeros(4,1);FLX3=zeros(4,1);KBTOT=0;Y2OX=0;
FQsW_L=0;FQsW_R=0;
pondloss=0;

IO.Y1=Y1;IO.Y2=Y2;IO.Y3=Y3;
IO.flyr1=flyr1;IO.flyr2=flyr2;IO.flyr3=flyr3;
IO.flyrb1=flyrb1;IO.flyrb2=flyrb2;IO.flyrb3=flyrb3;
IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
IO.Active=Active;
fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
fIO.pondloss=pondloss;
fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
fIO.FQsW_R=FQsW_R;fIO.FQsW_L=FQsW_L;

zbedo=z-msl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



makevideo=0;
%v=VideoWriter('r07fox0co30D10_R5','Motion JPEG AVI');
%v=VideoWriter('kro1dm_r3LARGE_R5Co15SLOPE','Motion JPEG AVI');
%v=VideoWriter('r07R1Co40___Co5','Motion JPEG AVI');
v=VideoWriter('ALL','Motion JPEG AVI');


%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if makevideo==1;open(v);end 
s=0;step=0;tic;
for t=1:tmax;  %iteration over the tmax time stpes
    
    %P.RSLR=8/1000/365;
%     if t==500        
% IO.Y2(310:360,100:140)=IO.Y2(310:360,100:140)+1.5;
% IO.Y2(150+310:150+360,200:240)=IO.Y2(150+310:150+360,200:240)-1.5;
%     end
%if t>250;P.RSLR=8/1000/365;end

%Swellwave direction
angleSWELL=angleSWELLserie(t);
%SeaWave direction
angleWIND=angleWINDserie(t);%rand(1)*360; %every time step a random direction
%if angleWIND>90 & angleWIND<270;P.wind=6.5
%Lentgh of event
dtO=dtOserie(t);
%Storm surge
Hsurge=surge(t);


%[dtO P.Ho Hsurge]
if t==1;dto=0.00001;else;dto=dtO;end   

% if mod(t,4)==0
%     IO.Y2(310:360,100:140)=IO.Y2(310:360,100:140)+0.5;
%     IO.Y2(310:460,1:50)=IO.Y2(310:460,1:50)+0.5;
% end

dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=limitdeltaz+1;maxup=limitmaxup+1;
        while maxdeltaz>limitdeltaz | maxup>limitmaxup
        if firstattemp==1;else;dt=dt/2*min(limitdeltaz/maxdeltaz,limitmaxup/maxup);end;firstattemp=0;
        if t<=2;dt=min(0.5*365,dt);end
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,P,dx,dt,zb,IO,fIO,Hsurge,angleSWELL,angleWIND,t);
        step=step+1; %this is how many time you called the function mainevolution step
        end

    %the partial updating step was succefull! Keep going
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
if t>0%==149;%>0;%==249;
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


z=-zb+(Yb+plyr*P.dlyr)+(Y1+Y2+Y3);
Y=Y1+Y2+Y3;
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
flyrU1=max(0,Y1)./Ytot;flyrU1(Ytot==0)=1;
flyrU2=max(0,Y2)./Ytot;flyrU2(Ytot==0)=0;
flyrU3=max(0,Y3)./Ytot;flyrU3(Ytot==0)=0;
zbed=z-msl;


%STORE_t200_r07_ALL(:,:,contatore)=zbed;

%ax1 = subplot(4,4,contatore);%s+1 %set(IM,'alphadata',~(A==0));s+1
ax1 = subplot(1,1,1);%s+1 %set(IM,'alphadata',~(A==0));s+1
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
%IM=imagesc(y,x,zbed);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
IM=imagesc(x,y,zbed');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
%cmp=demcmap([-4 0.5],256); %-3 1 P.Trange/2
cmp=demcmap([-3 P.Trange/2],256); %-3 1 P.Trange/2
colormap(ax1,cmp)
%colormap('jet')
caxis([-3 P.Trange/2]);
%colorbar('hori') 
%xlim([0 1.5])
title(strcat(num2str(time(t)),' years ',num2str(step)))
%title(strcat(num2str(time(t)),' years   (RSLR=',num2str(P.RSLR*365*1000),' mm/yr)'))




% ax1 = subplot(2,1,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(x,y,S');axis equal;set(gca,'YDir','normal');%colormap('jet');
% %xlim([0 1.5])
 
% 
% ax1 = subplot(2,1,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(x,y,1000*SSC2');axis equal;set(gca,'YDir','normal');%colormap('jet');
% caxis([0 20])
% colormap('jet')

%  ax3 = subplot(1,3,3);
%  IM=imagesc(y,x,1-flyrU1);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% cptcmap('mycmap1', 'mapping', 'direct'); 
%  caxis([0 1])
% % % 
%  
% ax1 = subplot(1,3,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,U);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% colormap('jet')
% caxis([0 0.5]);
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



% %  ax1 = subplot(1,3,3); %set(IM,'alphadata',~(A==0));
% % %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% % %IM=imagesc(y,x,deltaPW+0*Hs);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% % IM=imagesc(y,x,SSC2*1000);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% % %cmp=demcmap([-3 P.Trange/2],256); %-3 1
% % %colormap(ax1,cmp)
% % caxis([0 50]);colormap('jet')
% % %colorbar('hori') 
% %  %colorbar
% %  
% % % 
% % %  ax1 = subplot(1,3,2); %set(IM,'alphadata',~(A==0));
% % % %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% % % %IM=imagesc(y,x,deltaPW+0*Hs);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% % % IM=imagesc(y,x,waveANGLE);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% % % %cmp=demcmap([-3 P.Trange/2],256); %-3 1
% % % %colormap(ax1,cmp)
% % % caxis([-90 90]);
% % % %colorbar('hori') 
% % %  %colorbar
% % %  
% % 
% % 
% %  ax1 = subplot(1,3,2); %set(IM,'alphadata',~(A==0));
% % %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% % IM=imagesc(y,x,Hsea);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% % %cmp=demcmap([-3 P.Trange/2],256); %-3 1
% % %colormap(ax1,cmp)
% % caxis([0 0.3]);
%colorbar('hori') 
 %colorbar
%  S=2*double(S)-1; S(B>0)=0;
%  S(zbed>P.Trange/2)=0;
%  imagesc(y,x,S);axis equal;set(gca,'YDir','normal')
%  caxis([-1 1])
 
 
%   ax1 = subplot(1,3,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Hsea);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% %cmp=demcmap([-3 P.Trange/2],256); %-3 1
% %colormap(ax1,cmp)
% caxis([0 0.3]);
% %colorbar('hori') 
%  %colorbar
%  
 
%    ax2 = subplot(2,2,3);
%   IM=imagesc(y,x,UR);set(IM,'alphadata',~(A==0));axis equal;set(gca,'YDir','normal');%colormap('jet');%colorbar('hori') 
%   caxis([0 0.3]);colormap('jet')
%  colormap(ax2,flipud(hot))
 

%  
 %max(1000*SSC(:))
%   ax2 = subplot(2,3,5);
%   IM=imagesc(y,x,1000*SSC2);set(IM,'alphadata',~(A==0));axis equal;set(gca,'YDir','normal');%colormap('jet');%colorbar('hori') 
%   caxis([0 10]);colormap('jet')
%  colormap(ax2,flipud(hot))
%  %colorbar
 
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






%  ax3 = subplot(2,3,4);
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
%TO PLOT USA QUESTO DIOCANE 
checksum=[[(sumY1IN-sumY1)+sumFLUX1/dx^2]+fIO.FQsW_L+fIO.FQsW_R    [(sumY2IN-sumY2)+sumFLUX2/dx^2]-pondloss+KBTOT-Y2OX    [(sumY3IN-sumY3)+sumFLUX3/dx^2]];
if abs(checksum(1))>0.1 |  abs(checksum(2))>0.1 | abs(checksum(3))>0.1 ;checksum,pause;end

if makevideo==1;V=getframe(figure(1));writeVideo(v,V);end
pause(0.1)
end

end



end; %end of panel plotting
end;%of of repetitie cycle
if makevideo==1;close(v);end %UNCOMMENT THIS TO CREATE A VIDEO

