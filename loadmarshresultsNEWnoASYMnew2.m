clear;close all;clc

cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmap.cpt', cdata, ' ');

 cdata = [-1    50  255 255     0     255    255    255
          0    255 255 255     1     0    0    0];
 cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmapCLR.cpt', cdata, ' ');

cdata = [-1    205  165 0     0     0    0    0];
dlmwrite('mycmap1.cpt', cdata, ' ');


R=[2.5 5 7.5 10];
dx=10;





% load STOREr07ALL;ALL=STORE;
% load STOREr07UPLANDALL;ALLup=STORE;
% %B=squeeze(NOPOND(:,:,1))>0;
% %somma=sum(B(:)==1)
% Trange=0.7;
% L1=zeros(4,4);
% figure
% i=15
% 
% ax1 = subplot(2,2,1);
% G=squeeze(ALL(:,:,i));
% %G=[ones(201,300); G];
% G=[-2+0.7/2+0.001*dx*[200:-1:0]'*ones(1,300)+0.1*(rand(200+1,300)-0.5)*NaN ;G];
% 
% IM=imagesc(G');axis equal;set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
% colormap(ax1,cmp)
% caxis([-3 Trange/2]);
% 
% ax1 = subplot(2,2,2);
% Gup=squeeze(ALLup(:,:,i));
% IM=imagesc(Gup');axis equal;set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
% colormap(ax1,cmp)
% caxis([-3 Trange/2]);
% 
% ax1 = subplot(2,2,3);
% IM=imagesc(G'-Gup');axis equal;set(gca,'YDir','normal');%colormap('jet');
% %cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
% colormap('jet');xlim([200 700])
% colorbar
% caxis([-2 2])
% %caxis([-3 Trange/2]);
% 
% ax1 = subplot(2,2,4);
% A=G'*0;
% A(G'>0 & Gup'<0)=1;
% A(G'<0 & Gup'>0)=-11;
% IM=imagesc(A);axis equal;set(gca,'YDir','normal');%colormap('jet');
% %cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
% colormap('jet');xlim([200 700])
% colorbar
% caxis([-1 1])
% %caxis([-3 Trange/2]);








load STOREr07ALL;ALL=STORE;
load STOREr07ALLfox100;ALLup=STORE;
%B=squeeze(NOPOND(:,:,1))>0;
%somma=sum(B(:)==1)
Trange=0.7;
L1=zeros(4,4);
figure
i=15

ax1 = subplot(2,2,1);
G=squeeze(ALL(:,:,i));

IM=imagesc(G');axis equal;set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
colormap(ax1,cmp)
caxis([-3 Trange/2]);

ax1 = subplot(2,2,2);
Gup=squeeze(ALLup(:,:,i));
IM=imagesc(Gup');axis equal;set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
colormap(ax1,cmp)
caxis([-3 Trange/2]);

ax1 = subplot(2,2,3);
IM=imagesc(G'-Gup');axis equal;set(gca,'YDir','normal');%colormap('jet');
%cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
colormap('jet');
colorbar
caxis([-2 2])
%caxis([-3 Trange/2]);

ax1 = subplot(2,2,4);
A=G'*0;
A(G'>0 & Gup'<0)=1;
A(G'<0 & Gup'>0)=-11;
IM=imagesc(A);axis equal;set(gca,'YDir','normal');%colormap('jet');
%cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
colormap('jet');
colorbar
caxis([-1 1])
%caxis([-3 Trange/2]);


pause


% %load STOREr07ALL;ALL=STORE;
% load STOREr3UPLANDALL;ALL=STORE;
% %load STOREr3NOWAVE;ALL=STORE;
% %load STOREr3ALL;ALL=STORE;
% 
% %B=squeeze(NOPOND(:,:,1))>0;
% %somma=sum(B(:)==1)
% Trange=3;
% L1=zeros(4,4);
% figure
% for i=1:16
% ax1 = subplot(4,4,i);
% G=squeeze(ALL(:,:,i));
% IM=imagesc(G');axis equal;set(gca,'YDir','normal');%colormap('jet');
% %cmp=demcmap([-4 0.5],256); %-3 1 P.Trange/2
% cmp=demcmap([-3 Trange/2],256); %-3 1 P.Trange/2
% colormap(ax1,cmp)
% %colormap('jet')
% caxis([-3 Trange/2]);
% 
% end
% pause





figure
load STOREr07ALL;ALL=STORE;
load STOREr07NOWAVENOPOND;NOWAVENOPOND=STORE;
load STOREr07NOWAVENOPOND;REF=STORE;
A=squeeze(REF(:,:,1))>0;%the reference marsh
T=sum(A(:)==1);
A=squeeze(NOWAVENOPOND(:,:,1))>0 & squeeze(NOWAVENOPOND(:,:,1))<0.37;%the reference marsh
L1=zeros(4,4);
for i=1:16
subplot(4,4,i)
B=squeeze(ALL(:,:,i))>0 & squeeze(ALL(:,:,i))<0.37;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
%imagesc(B')
%colormap('jet')
L1(i)=(T-sum(B(:)==1))/T;
end

load STOREr07UPLANDALL;ALL=STORE;
load STOREr07UPLANDALL;NOWAVENOPOND=STORE;
load STOREr07NOWAVENOPOND;REF=STORE;
A=squeeze(REF(:,:,1))>0;%the reference marsh
T=sum(A(:)==1);
A=squeeze(NOWAVENOPOND(:,:,1))>0 & squeeze(NOWAVENOPOND(:,:,1))<0.37;%the reference marsh
L1up=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(ALL(:,:,i))>0 & squeeze(ALL(:,:,i))<0.37;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
%imagesc(B')
%colormap('jet')
L1up(i)=(T-sum(B(:)==1))/T;
end

for i=1:4
subplot(4,2,i*2-1)
plot(R,L1(i,:),R,L1up(i,:));
ylim([-0.4 1]);xlim([2.5 10])
end





load STOREr3ALL;ALL=STORE;
load STOREr3NOWAVENOPOND;NOWAVENOPOND=STORE;
load STOREr3NOWAVENOPOND;REF=STORE;
A=squeeze(REF(:,:,1))>0;%the reference marsh
T=sum(A(:)==1);
A=squeeze(NOWAVENOPOND(:,:,1))>0 & squeeze(NOWAVENOPOND(:,:,1))<1.6;%the reference marsh
L1=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(ALL(:,:,i))>0 & squeeze(ALL(:,:,i))<1.6;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
%imagesc(B')
%colormap('jet')
L1(i)=(T-sum(B(:)==1))/T;
end

load STOREr3UPLANDALL;ALL=STORE;
load STOREr3UPLANDALL;NOWAVENOPOND=STORE;
load STOREr3NOWAVENOPOND;REF=STORE;
A=squeeze(REF(:,:,1))>0;%the reference marsh
T=sum(A(:)==1);
A=squeeze(NOWAVENOPOND(:,:,1))>0 & squeeze(NOWAVENOPOND(:,:,1))<1.6;%the reference marsh
L1up=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(ALL(:,:,i))>0 & squeeze(ALL(:,:,i))<1.6;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
%imagesc(B')
%colormap('jet')
L1up(i)=(T-sum(B(:)==1))/T;
end

for i=1:4
subplot(4,2,i*2)
plot(R,L1(i,:),R,L1up(i,:));
ylim([-0.4 1]);xlim([2.5 10])
end














pause

%figure

load STOREr07ALL;ALL=STORE;
load STOREr07NOWAVE;NOWAVE=STORE;
load STOREr07NOWAVENOPOND;NOWAVENOPOND=STORE;
load STOREr07NOPOND;NOPOND=STORE;

A=squeeze(NOWAVENOPOND(:,:,1))>0;%the reference marsh
T=sum(A(:)==1);

L1=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(ALL(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L1(i)=(T-sum(B(:)==1))/T;
end
L1good=L1;

L2=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(NOPOND(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L2(i)=(T-sum(B(:)==1))/T; 
end

L3=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(NOWAVE(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L3(i)=(T-sum(B(:)==1))/T; 
end

L4=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(NOWAVENOPOND(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L4(i)=(T-sum(B(:)==1))/T; 
end

L5=zeros(4,4);
[Ddrn,idx]=bwdist(A==0);Ddrn=double(dx*Ddrn);%distance from natural channels
for i=1:16
%subplot(4,4,i)
B=squeeze(NOWAVENOPOND(:,:,i))>0;
B(Ddrn>50)=1;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
caxis([0 2])  
L5(i)=(T-sum(B(:)==1))/T; 
end

L2=L2-L4;
L3=L3-L4;
SUM=L2+L3+L4;
L2=L2./SUM.*L1;
L3=L3./SUM.*L1;
L4=L4./SUM.*L1;
L5=L5./SUM.*L1;

figure
for i=1:4
subplot(4,2,i*2-1)
X=[R; R; R; R];
Y=[ L4(i,:)-L5(i,:); L5(i,:);L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
area(X',Y');
ylim([0 1]);xlim([2.5 10])
end
legend('Drowning','channel','wave edge erosion','ponding')
xlabel('RSLR')
colormap('jet')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




load STOREr3ALL;ALL=STORE;
load STOREr3NOWAVE;NOWAVE=STORE;
load STOREr3NOWAVENOPOND;NOWAVENOPOND=STORE;
load STOREr3NOPOND;NOPOND=STORE;

A=squeeze(NOWAVENOPOND(:,:,1))>0;%the reference marsh
T=sum(A(:)==1);

L1=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(ALL(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
L1(i)=(T-sum(B(:)==1))/T;
end
L1good=L1;

L2=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(NOPOND(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L2(i)=(T-sum(B(:)==1))/T; 
end

L3=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(NOWAVE(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L3(i)=(T-sum(B(:)==1))/T; 
end

L4=zeros(4,4);
for i=1:16
%subplot(4,4,i)
B=squeeze(NOWAVENOPOND(:,:,i))>0;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
colormap('jet')
L4(i)=(T-sum(B(:)==1))/T; 
end

L5=zeros(4,4);
[Ddrn,idx]=bwdist(A==0);Ddrn=double(dx*Ddrn);%distance from natural channels
for i=1:16
%subplot(4,4,i)
B=squeeze(NOWAVENOPOND(:,:,i))>0;
B(Ddrn>50)=1;
I=(A-B);G=A*0+1;G(A==0)=0;G(I==1)=2;
L5(i)=(T-sum(B(:)==1))/T; 
end

L2=L2-L4;
L3=L3-L4;
SUM=L2+L3+L4;
L2=L2./SUM.*L1;
L3=L3./SUM.*L1;
L4=L4./SUM.*L1;
L5=L5./SUM.*L1;

for i=1:4
subplot(4,2,i*2)
X=[R; R; R; R];
Y=[ L4(i,:)-L5(i,:); L5(i,:);L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
area(X',Y');
ylim([0 1]);xlim([2.5 10])
end
legend('Drowning','channel','wave edge erosion','ponding')
xlabel('RSLR')
colormap('jet')









% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %figure
% A=squeeze(STORE_t200_r3_ALL(:,:,1))>0;
% T=sum(A(:)==1);
% 
% L1=zeros(4,4);
% L2=zeros(4,4);
% L3=zeros(4,4);
% L4=zeros(4,4);
% L5=zeros(4,4);
% 
% for i=1:16
% B=squeeze(STORE_t200_r3_ALL(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% L1(i)=(T-sum(B(:)==1))/T;
% end
% L1r3good=L1;
% 
% L2=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r3_NOPOND(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% L2(i)=(T-sum(B(:)==1))/T; 
% end
% 
% L3=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r3_NOWAVE(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% L3(i)=(T-sum(B(:)==1))/T; 
% end
% 
% L4=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r3_NOWAVENOPOND(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% L4(i)=(T-sum(B(:)==1))/T; 
% end
% 
% L5=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r3_NOWAVENOPONDNOCH(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% L5(i)=(T-sum(B(:)==1))/T; 
% end
% 
% 
% % R=[1:4];
% % %figure
% % for i=1:4
% % subplot(4,2,i*2-1)
% % plot(R,L1(i,:),R,L2(i,:),R,L3(i,:),R,L4(i,:))
% % ylim([0 1]);xlim([2.5 10])
% % end
% % legend('ponds+waves','no pond','no wave','no pond no wave')
% % xlabel('RSLR')
% 
% 
% 
% 
% L2=L2-L4;
% L3=L3-L4;
% SUM=L2+L3+L4;
% 
% L2=L2./SUM.*L1;
% L3=L3./SUM.*L1;
% L4=L4./SUM.*L1;
% L5=L5./SUM.*L1;
% 
% for i=1:4
% subplot(4,2,i*2)
% %plot(R,L1(i,:),R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
% X=[R; R; R; R];
% Y=[L5(i,:); L4(i,:)-L5(i,:); L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
% area(X',Y');
% ylim([0 1]);xlim([2.5 10])
% end
% legend('ALL','ALL','wave+drowning+widening','drowning+widening','summa')
% xlabel('RSLR')
% 
% 
% % figure
% % for i=1:4
% % subplot(4,2,i*2)
% % %plot(R,L1(i,:),R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
% % X=[R; R; R; R];
% % Y=[L5(i,:); L4(i,:)-L5(i,:); L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
% % area(X',Y');
% % ylim([0 1]);xlim([2.5 10])
% % end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% figure
% A=squeeze(STORE_t200_r07_ALLUPLAND(:,:,1)>0 & STORE_t200_r07_ALLUPLAND(:,:,1)<0.35);
% T=sum(A(:)==1)
% 
% L1=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r07_ALLUPLAND(:,:,i)>0 & STORE_t200_r07_ALLUPLAND(:,:,i)<0.35);
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% colormap('jet')
% caxis([0 2])   
% L1(i)=(T-sum(B(:)==1))/T;
% end
% 
% R=[1:4];
% for i=1:4
% subplot(4,2,i*2-1)
% plot(R,L1(i,:),R,L1good(i,:),'g')
% hold on;plot([R(1) R(end)],[0 0],'--k')
% ylim([-0.25 1])
% end
% legend('ponds+waves','waves','standard')
% xlabel('RSLR')
% 
% 
% L1=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r3_ALLUPLAND(:,:,i)>0 & STORE_t200_r07_ALLUPLAND(:,:,i)<0.35);
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% colormap('jet')
% caxis([0 2])   
% L1(i)=(T-sum(B(:)==1))/T;
% end
% 
% R=[1:4];
% for i=1:4
% subplot(4,2,i*2)
% plot(R,L1(i,:),R,L1r3good(i,:),'g')
% hold on;plot([R(1) R(end)],[0 0],'--k')
% ylim([-0.25 1])
% end
% legend('ponds+waves','waves','standard')
% xlabel('RSLR')




% Z=squeeze(NEWnoASYM_t200_r07_ALL(:,:,9));
% [S,AC,DIF]=findisolatedponds(Z,A,N,M,dx,zntwrk,zsill,distdr,minponddepth);
% sum(S(:)==1)/sum(Z(:)>0)

% L1=zeros(4,4);
% figure
% for i=1:16
% subplot(4,4,i)
% B=squeeze(ALL(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% IM=imagesc(G');
% set(IM,'alphadata',~(G'==0));
% colormap('jet')
% caxis([0 2])   
% L1(i)=(T-sum(B(:)==1))/T;
% %     Z=squeeze(NEWnoASYM_t200_r07_ALL(:,:,i));
% %     [S,AC,DIF]=findisolatedponds(Z,A,N,M,dx,zntwrk,zsill,distdr,minponddepth);
% %     sum(S(:)==1)/sum(Z(:)>0)
% end


% 