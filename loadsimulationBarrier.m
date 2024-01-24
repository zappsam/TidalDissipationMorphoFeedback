clear;close all;clc
load DATA;

dt=2;
figure;
sl=0;
for i=1:1: 3369;%500;
    if i<2000;slr=1/1000;else;slr=3/1000;end
    sl=sl+slr*dt;
z=DATA(i).z;%-sl;
% imagesc(-z)
% axis equal;set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-10 3],256); %-3 1 P.Trange/2
% colormap(cmp)
% caxis([-10 3]);
Z=prctile(-z',90);%median(-z')
hold off;plot(Z);

a=find(Z>0);
hold on
plot(a(end),Z(a(end)),'or')
p(i)=a(end);
ylim([-10 10]);xlim([0 300])


for j=1:150;
Z=-z(:,j);
a=find(Z>0);
pp(j)=a(end);
end
%p2(i)=mean(pp);
p2(i)=median(pp);


title(i*dt)
%pause(0.01)
end

figure
dx=0.200;
time=[1:i]*dt;
plot(time,p*dx,time,p2*dx,time,50-1/1000*1/(0.5/1000)*1/1000*time,'-k');ylim([30 70])