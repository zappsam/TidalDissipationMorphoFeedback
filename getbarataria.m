function [dx,N,M,A,z,x,y,SPCLcell]=getmarshgeometry;
%geometric parameters

%dx=200;
dx=100;
%load Abarataria;A=Abarataria(end:-1:1,:);
%load hbarataria;z=-double(hbarataria(end:-1:1,:));
%A(end-2:end,:)=0;

%load hterrebonne;z=-double(hterrebonne(end:-1:1,:));
%load hterrebonnedx100m;z=-double(hterrebonnedx100m(end:-1:1,:));
%z=z(30:end-100,100:end-100);

%load AterreH;z=-double(AterreH(end:-1:1,:));
load AterreH1930;z=-double(AterreH1930(end:-1:1,:));
z=z(1:end-0,100:end);

z=z(1:2:end,1:2:end);
A=z*0+1;
A(1,:)=2;

%figure;imagesc(A);pause
[N,M]=size(A);


%%%%%%%%%%%cell types

rivermouthfront=[];

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;

%landhighboudary
z(end,:)=-10;
%sea
z(1:2,:)=100;

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