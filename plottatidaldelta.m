clear;close all;clc
load zbedi
% load zbed
% load zbedis
% load zbedis2
%load zbedG
load zbedG2
%load zbedG3
zbed=zbedG2;
load zbedD3D;zbedD3D=-zbedD3D;
A=zbed;
[N,M]=size(A);
Z=A*0;
Z(20:68,N/2-50:N/2+50)=999;
a=find(Z==999);
%zbedD3D(a)=NaN;



Z=A*0;
Z(83:end,:)=999;
b=find(Z==999 & isfinite(zbedD3D));
%zbedD3D(b)=NaN;

% ax1 = subplot(1,2,1);%s+1 %set(IM,'alphadata',~(A==0));s+1
% IM=imagesc(-zbedi);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-10 3],256); %-3 1 P.Trange/2 
% colormap(ax1,cmp) 
% caxis([-10 3]+1);%
ax1 = subplot(1,2,1);%s+1 %set(IM,'alphadata',~(A==0));s+1
IM=imagesc(-zbedD3D);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
%IM=imagesc(-zbedis);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-10 3],256); %-3 1 P.Trange/2 
colormap(ax1,cmp) 
caxis([-10 3]+1);%

ax1 = subplot(1,2,2);%s+1 %set(IM,'alphadata',~(A==0));s+1
IM=imagesc(-zbed);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-10 3],256); %-3 1 P.Trange/2 
colormap(ax1,cmp) 
caxis([-10 3]+1);%

 dif=zbedi-zbed;dif(dif<0)=0;
 sum(dif(a))
sum(zbedi(a))-sum(zbed(a))
% nansum(zbedi(a))-nansum(zbedis(a))
% nansum(zbedi(a))-nansum(zbedis2(a))

 dif=zbedi-zbedD3D;dif(dif<0)=0;
 sum(dif(a))
sum(zbedi(a))-sum(zbedD3D(a))

figure;
[N,X]=hist(zbedi(a),[-1:0.15:12]);N=N/sum(N);plot(X,N,'.-g')
hold on;[N,X]=hist(zbedD3D(a),[-1:0.15:12]);N=N/sum(N);plot(X,N,'.-')
hold on;[N,X]=hist(zbed(a),[-1:0.15:12]);N=N/sum(N);plot(X,N,'.-r')
xlim([-1 12])

legend('initial','Delft3D','CM2D')
