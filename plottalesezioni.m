clear;close all;clc

load WEST2010_measured 
load WEST2010_model
load WEST1900_model 
load WEST2500_R25co2_model
load WEST2500_R5co5_model
dx=2;

z1=WEST2010_measured;
z2=WEST2010_model;
%z2=WEST1900_model;
%z2=WEST2500_R25co2_model;%
%z2=WEST2500_R5co5_model;

ax1 = subplot(1,3,1); %set(IM,'alphadata',~(A==0));
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(-z1);axis equal;set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 3.1/2],256); %-3 1
colormap(ax1,cmp)
caxis([-3 3.1/2]);
%colorbar('hori') 
 %colorbar
 
 
 ax1 = subplot(1,3,2); %set(IM,'alphadata',~(A==0));
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(-z2);axis equal;set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 3.1/2],256); %-3 1
colormap(ax1,cmp)
caxis([-3 3.1/2]);
z2(z1==-2)=NaN;

jo=421;io=241;
xs=[0:50];i=io+xs;j=jo+floor(xs*0.6);
s=sub2ind(size(z1),i,j);
figure;
subplot(1,2,1);plot(xs,-z1(s),'.-',xs-3,-z2(s),'.--');xlim([0 45]);ylim([-2.5 1.5])

jo=355;io=420;
xs=[0:30];i=io+xs;j=jo+floor(xs*0.9);

s=sub2ind(size(z1),i,j);
subplot(1,2,2);plot(xs,-z1(s),'.-',xs-5,-z2(s),'.--');xlim([0 25]);ylim([-2.5 1.5])

% z1(s)=NaN;
% figure;imagesc(z1)




z=z2;

thwg=findthalweg(z(end:-1:1,:)+1.8,0.1,3000);
A=0*z+1;A(z<0)=0;
W=calculatewidth(A(end:-1:1,:),dx);

figure;imagesc(z)
figure;imagesc(W)
figure;imagesc(thwg)

th=find(thwg>0 & W>0);
th=th(1:5:end);

[N,M]=size(A);
for i=1:length(th)
    [I,J] = ind2sub(size(A),th(i));
    a=0;
    for x=-3:1:3;
        for y=-3:1:3
            a=max(a,z(min(N,max(1,I+x)),min(M,max(1,J+y))));%Dh(i)=z1(I,J);
        end
    end
    Dh(i)=a;
end
Wh=W(th);


[a,TXT,RAW]=xlsread('WestcreeksWD.xlsx');

% WD_model=[Dh; Wh'];
% save WD_model WD_model
load  WD_model

%WD_modelR5=[Dh; Wh'];
%save WD_modelR5 WD_modelR5
load  WD_modelR5

D=a(:,1);W=a(:,2);
Di=[0 4];Wi=Di*5.8;
figure;plot(W,D,'ob',Wi,Di,'-k',Wh,Dh,'.r',WD_model(2,:),WD_model(1,:),'.g',WD_modelR5(2,:),WD_modelR5(1,:),'.y');ylim([0 4]);xlim([0 110])%,W,D+1.3,'.'



