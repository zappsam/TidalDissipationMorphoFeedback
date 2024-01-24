function [dx,N,M,A,z,x,y,SPCLcell]=getmarshgeometry;
%geometric parameters

dx=30;
load zDelaware2010dems
zDelaware2010dems(isnan(zDelaware2010dems))=0;
zDelaware2010dems=zDelaware2010dems*0-5;
[M,N]=size(zDelaware2010dems);
zDelaware2010dems=-ones(N,1)*[0:M]/M*5;
zDelaware2010dems=zDelaware2010dems';

 z=-zDelaware2010dems(end:-1:1,:);
 z=double(z);
% % z(z==1)=1;
% marsh=find(z==0);
% z(z==0)=-0.5;
% z(z>0 & z<1)=1;
% %z(z>0)=1;
% 
% %z=z+0.5;
% z(marsh)=-0.5;


A=zDelaware2010dems(end:-1:1,:);
A=double(A);
A=A*0+1;
A(isnan(A))=0;
% A(end-100:end,:)=0;
% A(:,1:5)=0;

%figure;imagesc(A);pause

[N,M]=size(A);


%%%%%%%%%%%cell types

rivermouthfront=[];

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;


z(A==0)=NaN;

x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
%%%%%%%%%%%%%%%%%%


%z=z(end:-1:1,:);
%A=A(end:-1:1,:);



% z=imrotate(z,45,'nearest');z=z(480:end,750:end-860);
% A=imrotate(A,45,'nearest');A=A(480:end,750:end-860);
% [N,M]=size(A);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
% 
% 
% 






% dx=100;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';A=interp2([0:N-1],[0:M-1]',A',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;