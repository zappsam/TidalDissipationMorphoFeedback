function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry(P);
%shallower for Co=60mg/l
% %geometric parameters
% dx=50*2;
% N=200/2; %x
% M=200/2; %y%%%
%%
% %%
% dx=50/4;%
% N=(200-70)*4; %x
% M=100*4  *0.6; %y

% dx=50/5*2;%
% N=(200-86)*2*5/2; %x
% M=100*0.5*2*5/2; %y

% dx=5*2; 
% N=1000/2*1;%/2; %x
% M=600/2*1;%/4; %y
%
dx=20;%10*2;%%%20 m resolution to make tidal dissipation work 1/19/22
N=2200/2*1;%/2; %x  2200/2*1
M=400/2*1;%/4; %y    500/2*1
% N=2000/2*1;%/2; %x
% M=600/2*1;%/4; %y
%%%%%%%%%%%cell types
A=ones(N,M);
%sea b.c.
A(end,:)=2;
%A(end-1:end,:)=2;
%river b.c.
rivermouthfront=[];

% mouthW=1; %2
% A(1:2,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:2,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% S=A*0;S(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(S==1);clear S;

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;

%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping   0.5
%slope=0.2/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1;%sloping   0.5

%slope=1/1*0.2/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1;%sloping%ORIGINAL

%slope=1/1*0.1/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1; %change SLOPE %USED for tst58
%slope=1/1*0.05/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1; %tst60 slope
slope=1/1*0.0175/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1;
%slope=0;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1; 
%slope=.0000001;z=[0:N-1].^(slope)*dx;z=z'*ones(1,M)-1+1; %TRY and make
%curved slope                        %1.5 works to drop down elevation
%*********uncomment next 2 lines for curved slope*******************
% slopedepth=1;%1.2%approximate depth of depression, 1.2 works well to fill in 1000 yrs
% slope=.01;z=-(0.2.^([0:N-1].*slope))+slopedepth;z=z'*ones(1,M)-1+1;%z was previously multiplied by dx

%slope impacts how long it curves for


%z=z*0+5;
%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping
%zin=5;z=zin*ones(N,M);%flat

%%%make bowl depression%%%*******uncomment for bowl
% bowldepth=12;%0.6;  %approximate depth of depression (m)-->was as low as 0.15 for 10m grid
% bowl=((bowldepth/6/1000/10)*((1:length(y))-(length(y)./2)).^2); %  zbowl=z*(bowl').*dx;z=z'*ones(1,M)-1+1;
% zbowl=0*z;
% for row=1:length(x)
%     zbowl(row,:)=(z(row,:))-bowl;%.*(z(row,:));
% end
% z=zbowl;
%%%%%%%%%%%%%%%

%%%%%%%%%%make v shaped depression
VV=(1:length(y))-length(y)/2;
VVV=abs(VV);vshape=6*1.5*VVV./(length(y)/2);%multiplier is m depth of bowl(multiplier was 2 up to tst81)  %15multiplier + 9exponent-->drowning
zvshape=0*z;
for row=1:length(x)
    zvshape(row,:)=(z(row,:))-vshape/(((row/length(x))+1).^5);%.*(z(row,:));%multiplier was 1.3 through tst87
    %zvshape(row,:)=(z(row,:))-vshape/(((row/length(x))+1).^2);%.*(z(row,:));
end
%%upland 5km area
vshape=24*1.5*VVV./(length(y)/2);%multiplier is m depth of bowl(multiplier was 2 up to tst81)  %15multiplier + 9exponent-->drowning
for row=1:250
    zvshape(row,:)=(z(row,:))-vshape/(((row/length(x))+1).^6.5);%.*(z(row,:));
    %zvshape(row,:)=(z(row,:))-vshape/(((row/length(x))+1).^2);%.*(z(row,:));
end
z=zvshape;
%%%%%%%%%%
%drop down
z=z+0.3; %z=z+0.8; up to tst81  *1.4
%shift up
%z=z-0.5; %tst69
%z=z-0.2;
% slopebreak=30;
% z=[1:N];
% slope=1/1000;z(1:N-slopebreak)=[1:N-slopebreak]*slope*dx;
% slope=1/500;z(slopebreak+1:N)=+z(slopebreak)+[slopebreak+1:N]*slope*dx;
% z=z'*ones(1,M);%sloping

% %barriers
%pos=1;bwidth=1;%250/dx;
%A(end-pos:end-pos+bwidth,1:M/2-15)=0;
%A(end-pos:end-pos+bwidth,M/2+15:end)=0;
%z(end-pos:end-pos+bwidth,1:M/2-1)=-5.1;
%z(end-pos:end-pos+bwidth,M/2+1:end)=-5.1;

%z=z+1; %DROP EVERYTHING DOWN 1m**********REMOVE THIS

%z(z>1)=1;
%channel
%z(end-200:end,M/2-3:M/2+3)=2.5;%channel ***used up through tst 81

z(end-600:end,M/2-3:M/2+3)=3.5; % for 20m res 
%z(end-400:end,M/2-6:M/2+6)=2.5; % for 10m res 
%channel%z(:,M/2-2:M/2+2)=2;%channel

%z(1:end-600,M/2-3:M/2+3)=z(1:end-600,M/2-3:M/2+3)+1; % for 20m res **CHANNEL EXTENSION
% z(:,M/2)=2.5;%channel
% z(:,M/2-1)=2;z(:,M/2+1)=2;
% z(:,M/2-2)=1.2;z(:,M/2+2)=1.2;
% z(:,M/2-3)=0.5;z(:,M/2+3)=0.5;

%z(end-50:end,M/2-5:M/2+5)=3;%for tst58_3_by_5!!!!!
%^USED FOR tst700yr
%z(end-200:end,M/2-5:M/2+5)=5;%1/26/22 channel
%z(:,M/2-5:M/2+5)=50; %long channel
%z(1:150,1:2)=-2;

%THIS IS THE GOOD ONE
%z(end-pos-bwidth:end-pos,:)=4;%3;
%z(end-pos+1:end,:)=z(end-pos+1:end,:)+10;
%z(end-pos:end-pos+bwidth,1:M/2-4)=-2.1;
%z(end-pos:end-pos+bwidth,M/2+4:end)=-2.1;

%z(end-20:end,M/2-2:M/2+2)=3;

%UPLAND ADJUSTMENT
%z(1:50,1:95)=-12;z(50:100,1:90)=-12;z(100:150,1:80)=-12;z(150:200,1:70)=-12;z(200:250,1:60)=-12;
%z(1:50,105:end)=-12;z(50:100,110:end)=-12;z(100:150,120:end)=-12;z(150:200,130:end)=-12;z(200:250,140:end)=-12;
z(z<-4)=-4+(z(z<-4)+4)*0.4; %multiplier at end scales values (1 = no change), number less than 1 creates terrace
%z(z<-4)=1*(min(min(z))+4)+z(z<-4)
z(1:250,1:90)=-100;z(1:250,110:end)=-100;
%random
z=z+2*(rand(N,M)-0.5)*0.2;% *0.5   *0.2;%made moltiplier 0.05 from 0.2****

% z(:,M/2)=2.5;%channel
% z(:,M/2-1)=2;z(:,M/2+1)=2;
% z(:,M/2-2)=1.2;z(:,M/2+2)=1.2;
% z(:,M/2-3)=0.5;z(:,M/2+3)=0.5;

% z(100:end,M/2)=2.5;%channel
% z(100:end,M/2-1)=2;z(100:end,M/2+1)=2;
% z(100:end,M/2-2)=1.2;z(100:end,M/2+2)=1.2;
% z(100:end,M/2-3)=0.5;z(100:end,M/2+3)=0.5;
%river
%z(A==10)=P.hmouth;
%z(1:40,M/2-mouthW:M/2+mouthW)=P.hmouth;
%z(end-100:end,:)=3;


%z(end-10:end-2,:)=-0.3;
%z(end-10:end-2,100:200)=10;
%sea
z(end-1:end,:)=10;

%%%%%%%%%%%%%%%%
%z(end-50:end,:)=5;

% z(z<0)=-0.5;
% z(z>=0)=2;

% 
% %River
% mouthW=1; %2
% A(1:10,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:10,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% SPC=A*0;SPC(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(SPC==1);clear SPC;
% SPCLcell=struct;
% SPCLcell.rivermouthfront=rivermouthfront;
% %river
% z(A==10)=P.hmouth;
% z(1:52,M/2-mouthW:M/2+mouthW)=P.hmouth;



%load NEWnoASYM_t200_r07_ALL NEWnoASYM_t200_r07_ALL
%z=-squeeze(NEWnoASYM_t200_r07_ALL(:,:,8));
%load NEWnoASYM_t200_r07_NOWAVE NEWnoASYM_t200_r07_NOWAVE
%z=-squeeze(NEWnoASYM_t200_r07_NOWAVE(:,:,8));


% [N,M]=size(z);
% dx=50;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% A=z*0+1;
% A(end,:)=2;
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;

%A(1:200,:)=0;


% load mudfltaNOPONDS
% z=-z+msl;
% 
% %load zbedo;
% %load zbedor3LARGE; z=-zbedor3LARGE;




% load zbedor07; z=-zbedor07;
% [N,M]=size(z);
% G=-0.7/2-0.001*dx*[200:-1:0]'*ones(1,M)+0.1*(rand(200+1,M)-0.5);
% z=[G; z];
% [N,M]=size(z);
% A=z*0+1;
% A(end,:)=2;
% x=[0:N-1]'*dx/1000;y=[0:M-1]*dx/1000;

% load zbedor07; z=-zbedor07;
% [N,M]=size(z);
% %G=-3/2-0.001*dx*[200:-1:0]'*ones(1,M)+0.1*(rand(200+1,M)-0.5);
% %z=[G; z];
% %[N,M]=size(z);
% A=z*0+1;
% A(end,:)=2;
% x=[0:N-1]'*dx/1000;y=[0:M-1]*dx/1000;

%z=z(:,end:-1:1);

% load zbedo;z=-zbedo;
% dx=5;z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% [N,M]=size(z);
% A=z*0+1;
% % A(end,:)=2;




z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%


msl=0;
zbedo=-z-msl;
Active=zbedo<P.Trange/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%wave boundary conditions
%the lateral boundaries for waves
%postivie is left boundary; negative is right boudndary; 
%1 is no-gradient in wave. THIS IS ALSO no-gradient in wave-indcued lateral
%sediment transport
%2 is zero wave height . Also implies no sand transport at inlet
AW=A*0;%the internal points
%right boundary
AW(:,end)=-1;
%left boundary
AW(:,1)=1;
%AW(1:50,1)=1;
%%%%%%%%%%%%%%%%%%%%%
