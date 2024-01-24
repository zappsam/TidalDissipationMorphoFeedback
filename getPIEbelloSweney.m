function [dx,N,M,A,z,x,y,SPCLcell]=getPIE;
%geometric parameters

load zSW3;
zSW3=zSW3(1:end-10,:);
z=-zSW3;
z=double(z);




dx=2;

% dx=4;
% z=z(1:2:end,1:2:end);
% z(177,39)=-999;z(176,40)=-999;
% z(145:150,1:10)=-999;
% z(144,17)=-999;

%attach a piece ofm marhs left of top channel mouth in Sweney
z(end-25:end,132:141)=-1.35;

% figure;imagesc(z)
% pause



[N,M]=size(z);

%%%%%%%%%%%cell types
A=ones(N,M);

A(z==-999 | z<-1.9)=0;
A(end,141:155)=2;
%A(end,156:end)=0;

A(isnan(z))=0;

z(A==0)=-2;


% figure;imagesc(z);axis equal;set(gca,'YDir','normal');
% caxis([-3 3.1/2]);
% figure;imagesc(A);axis equal;set(gca,'YDir','normal');
% pause

%z(A==1)=z(A==1)-0.1;


%sea b.c.
%river b.c.
rivermouthfront=[];
SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;



%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;



%increase resolution
% N=N*2;M=M*2;dx=dx/2;
% xn=[0:N-1]*dx/1000;yn=[0:M-1]*dx/1000;
% 
% size(A)
% z=interp2(y',x,z,yn',xn);
% 
% A=ones(N,M);
% 
% A(:,end)=0;
% A(95:115,end)=2;
% 
% A(isnan(z))=0;
% 
% A(95*2:115*2,end)=2;
% 
% x=xn;y=yn;







%figure;imagesc(A);pause


%z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%


% load PIEcut
% d=PIEcut;
% d(144:150,1:17)=NaN;
% d(175:177,39:42)=NaN;
% PIEcutfill=d;
% save PIEcutfill PIEcutfill
% figure;
% ax1 = subplot(1,1,1); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,-d);axis equal;set(IM,'alphadata',~(d==0));set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-3 2.7/2],256); %-3 1
% colormap(ax1,cmp)
% caxis([-3 2.7/2]);