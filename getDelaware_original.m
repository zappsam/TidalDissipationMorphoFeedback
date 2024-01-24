function [dx,N,M,A,z,x,y,SPCLcell]=getmarshgeometry;
%geometric parameters

dx=30;
load zDelaware2010s
z=zDelaware2010s(end:-1:1,:);
z=double(z);
z(z==1)=1;
z(z==0)=-0.5;

A=zDelaware2010s(end:-1:1,:);
A=double(A);
A(A==0)=1;
A(isnan(A))=0;

%figure;imagesc(z);pause

[N,M]=size(A);


%%%%%%%%%%%cell types

rivermouthfront=[];

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;


z(A==0)=NaN;

x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
%%%%%%%%%%%%%%%%%%



% dx=100;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';A=interp2([0:N-1],[0:M-1]',A',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;