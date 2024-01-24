function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,S]=initializegeometry(P);

%geometric parameters
dx=200/1;
N=140*1; %x
M=150*1; %y


%%%%%%%%%%%cell types
A=ones(N,M);
%sea b.c.
A(end,:)=2;
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


slope=0.6/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-5;%sloping
%zin=5;z=zin*ones(N,M);%flat

% slopebreak=30;
% z=[1:N];
% slope=1/1000;z(1:N-slopebreak)=[1:N-slopebreak]*slope*dx;
% slope=1/500;z(slopebreak+1:N)=+z(slopebreak)+[slopebreak+1:N]*slope*dx;
% z=z'*ones(1,M);%sloping

%barriers
pos=20*dx/dx;bwidth=1000/dx;
%z(end-pos:end-pos+bwidth,1:M/2-2)=-2.1;
%z(end-pos:end-pos+bwidth,M/2+2:end)=-2.1;
%z(end-pos:end-pos+bwidth,1:M/2-1)=-5.1;
%z(end-pos:end-pos+bwidth,M/2+1:end)=-5.1;

z(end-pos-bwidth:end-pos,:)=1;%3;
z(end-pos+1:end,:)=z(end-pos+1:end,:)+2;

%random
z=z+2*(rand(N,M)-0.5)*0.1;

%river
%z(A==10)=P.hmouth;
%z(1:40,M/2-mouthW:M/2+mouthW)=P.hmouth;

%sea
z(end-1:end,:)=50;

z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


S=0*A;  %the status of a cell; (for ponds dynamics). 1 is pond  0 is not a pond



msl=0;