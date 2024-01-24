function [dx,N,M,A,z,x,y,SPCLcell]=getmarshgeometry;
%geometric parameters
dx=50/1;
N=400*4; %x
M=50*1; %y


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


%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping   0.5
slope=0.01/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1;%sloping   0.5
%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping
%zin=5;z=zin*ones(N,M);%flat

% slopebreak=30;
% z=[1:N];
% slope=1/1000;z(1:N-slopebreak)=[1:N-slopebreak]*slope*dx;
% slope=1/500;z(slopebreak+1:N)=+z(slopebreak)+[slopebreak+1:N]*slope*dx;
% z=z'*ones(1,M);%sloping

%barriers
pos=35*dx/dx;bwidth=1000/dx;
%z(end-pos:end-pos+bwidth,1:M/2-2)=-2.1;
%z(end-pos:end-pos+bwidth,M/2+2:end)=-2.1;
%z(end-pos:end-pos+bwidth,1:M/2-1)=-5.1;
%z(end-pos:end-pos+bwidth,M/2+1:end)=-5.1;


%THIS IS THE GOOD ONE
z(end-pos-bwidth:end-pos,:)=1;%3;
z(end-pos+1:end,:)=z(end-pos+1:end,:)+5;

%z(end-20:end,M/2-2:M/2+2)=3;

slope=0.2/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)+2;%sloping   0.5

%z(end-100:end,24:26)=20;

%z(1:80,N/2-10:N/2+10)=0;

%random
z=z+2*(rand(N,M)-0.5)*0.5;

%river
%z(A==10)=P.hmouth;
%z(1:40,M/2-mouthW:M/2+mouthW)=P.hmouth;

%sea
%z(end-1:end,:)=50;

z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%



% dx=100;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';A=interp2([0:N-1],[0:M-1]',A',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;