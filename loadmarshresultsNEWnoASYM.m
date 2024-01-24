clear;close all;clc

load SEDt
load SEDtx
load SEDtx2
load Ut
load Utx
load Utx2

figure;
subplot(4,2,1);plot(Utx(1:200,1:4));legend('50','25','10','5');%ylim([-5 10])
subplot(4,2,3);plot(Utx(1:200,4+1:4+4));%ylim([-5 10])
subplot(4,2,5);plot(Utx(1:200,8+1:8+4));%ylim([-5 10])
subplot(4,2,7);plot(Utx(1:200,12+1:12+4));%ylim([-5 10])

subplot(4,2,2);plot(Utx2(1:200,1:4));legend('50','25','10','5');%ylim([-5 10])
subplot(4,2,4);plot(Utx2(1:200,4+1:4+4));%ylim([-5 10])
subplot(4,2,6);plot(Utx2(1:200,8+1:8+4));%ylim([-5 10])
subplot(4,2,8);plot(Utx2(1:200,12+1:12+4));%ylim([-5 10])


N=300;M=500;dx=10;tmax=200;
figure;
subplot(4,2,1);plot(SEDtx(1:200,1:4)/(N*M*dx^2)/tmax*1000);legend('50','25','10','5');ylim([-5 2])
subplot(4,2,3);plot(SEDtx(1:200,4+1:4+4)/(N*M*dx^2)/tmax*1000);ylim([-5 2])
subplot(4,2,5);plot(SEDtx(1:200,8+1:8+4)/(N*M*dx^2)/tmax*1000);ylim([-5 2])
subplot(4,2,7);plot(SEDtx(1:200,12+1:12+4)/(N*M*dx^2)/tmax*1000);ylim([-5 2])

subplot(4,2,2);plot(SEDtx2(1:200,1:4)/(N*M*dx^2)/tmax*1000);legend('50','25','10','5');ylim([-5 2])
subplot(4,2,4);plot(SEDtx2(1:200,4+1:4+4)/(N*M*dx^2)/tmax*1000);ylim([-5 2])
subplot(4,2,6);plot(SEDtx2(1:200,8+1:8+4)/(N*M*dx^2)/tmax*1000);ylim([-5 2])
subplot(4,2,8);plot(SEDtx2(1:200,12+1:12+4)/(N*M*dx^2)/tmax*1000);ylim([-5 2])

figure
SEDtx=tsmovavg(SEDtx,'s',20,1);
G=785*diff(SEDtx)/dx./(Utx(2:200,:)*24*3600*365)*1000;
subplot(4,2,1);plot(G(:,1:4))
subplot(4,2,3);plot(G(:,4+1:4+4))
subplot(4,2,5);plot(G(:,8+1:8+4))
subplot(4,2,7);plot(G(:,12+1:12+4))

pause
load NEWnoASYM_t200_r07_ALL NEWnoASYM_t200_r07_ALL
load NEWnoASYM_t200_r07_NOWAVE NEWnoASYM_t200_r07_NOWAVE

Trange=0.7;
dx=1;
N=500;M=300;
zsill=Trange/2;  
zpondcr=-0.2;%P.Trange/4;%base of new pond formation with respect to MSL
minponddepth=0.1;%1; %minimum depth to define a pond after you identified the "lakes"
zntwrk=(Trange/2)*0.2;%P.Trange/2*0.9;%P.Trange/2-0.3;%0.3; %depth above msl that defines the channel network.  the smaller the harder to drain!
distdr=NaN;


A=squeeze(NEWnoASYM_t200_r07_ALL(:,:,1))>0;
T=sum(A(:)==1);


Z=squeeze(NEWnoASYM_t200_r07_ALL(:,:,9));
[S,AC,DIF]=findisolatedponds(Z,A,N,M,dx,zntwrk,zsill,distdr,minponddepth);
sum(S(:)==1)/sum(Z(:)>0)

L1=zeros(4,4);
figure
for i=1:4
subplot(4,4,i)
B=squeeze(NEWnoASYM_t200_r07_ALL(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
IM=imagesc(G');
set(IM,'alphadata',~(G'==0));
colormap('jet')
caxis([0 2])   
L1(i)=(T-sum(B(:)==1))/T;

    Z=squeeze(NEWnoASYM_t200_r07_ALL(:,:,i));
    [S,AC,DIF]=findisolatedponds(Z,A,N,M,dx,zntwrk,zsill,distdr,minponddepth);
    sum(S(:)==1)/sum(Z(:)>0)

end


% 
% pause
% 
% 
% 
% load STORE_t200_r07_ALL
% load STORE_t200_r07_NOPOND
% load STORE_t200_r07_NOPONDNOWAVE
% load STORE_t200_r07_NOPONDNOWAVENOCH
% load STORE_t200_r07_NOWAVE
% load STORE_t200_r07_ALLUPLAND
% 
% load STORE_t200_r3_ALL
% load STORE_t200_r3_ALLscarp2m
% load STORE_t200_r3_NOWAVE
% load STORE_t200_r3_NOPONDNOWAVE
% load STORE_t200_r3_NOPONDNOWAVENOCH
% load STORE_t200_r3_NOPOND
% load STORE_t200_r3_ALLUPLAND
% 
% a=STORE_t200_r3_NOPONDNOWAVE;
% STORE_t200_r3_NOPONDNOWAVE=STORE_t200_r3_NOPOND;
% STORE_t200_r3_NOPOND=a;
% 
% 
% A=squeeze(STORE_t200_r07_NOPONDNOWAVE(:,:,1))>0;
% T=sum(A(:)==1)
% 
% figure;imagesc(A);
% 
% 
% L1=zeros(4,4);
% figure
% for i=1:16
% subplot(4,4,i)
% B=squeeze(STORE_t200_r07_ALL(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% %IM=imagesc(G');
% %set(IM,'alphadata',~(G'==0));
% colormap('jet')
% caxis([0 2])   
% L1(i)=(T-sum(B(:)==1))/T;
% end
% L1good=L1;
% 
% L2=zeros(4,4);
% %figure
% for i=1:16
% subplot(4,4,i)
% B=squeeze(STORE_t200_r07_NOPOND(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% %IM=imagesc(G');
% %set(IM,'alphadata',~(G'==0));
% colormap('jet')
% caxis([0 2])  
% L2(i)=(T-sum(B(:)==1))/T; 
% end
% 
% L3=zeros(4,4);
% %figure
% for i=1:16
% subplot(4,4,i)
% B=squeeze(STORE_t200_r07_NOWAVE(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% %IM=imagesc(G');
% %set(IM,'alphadata',~(G'==0));
% colormap('jet')
% caxis([0 2])  
% L3(i)=(T-sum(B(:)==1))/T; 
% end
% 
% 
% L4=zeros(4,4);
% %figure
% for i=1:16
% subplot(4,4,i)
% B=squeeze(STORE_t200_r07_NOPONDNOWAVE(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% %IM=imagesc(G');
% %set(IM,'alphadata',~(G'==0));
% colormap('jet')
% caxis([0 2])  
% L4(i)=(T-sum(B(:)==1))/T; 
% end
% 
% L5=zeros(4,4);
% %figure
% for i=1:16
% subplot(4,4,i)
% B=squeeze(STORE_t200_r07_NOPONDNOWAVENOCH(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% %IM=imagesc(G');
% %set(IM,'alphadata',~(G'==0));
% colormap('jet')
% caxis([0 2])  
% L5(i)=(T-sum(B(:)==1))/T; 
% end
% 
% 
% R=[2.5 5 7.5 10];
% for i=1:4
% subplot(4,2,i*2-1)
% plot(R,L1(i,:),R,L2(i,:),R,L3(i,:),R,L4(i,:),R,L5(i,:))
% ylim([0 1]);xlim([2.5 10])
% end
% legend('ponds+waves','no pond','no wave','no pond no wave')
% xlabel('RSLR')
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
% subplot(4,2,i*2-1)
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
% B=squeeze(STORE_t200_r3_NOPONDNOWAVE(:,:,i))>0;
% I=(A-B);
% G=A*0+1;
% G(A==0)=0;
% G(I==1)=2;
% L4(i)=(T-sum(B(:)==1))/T; 
% end
% 
% L5=zeros(4,4);
% for i=1:16
% B=squeeze(STORE_t200_r3_NOPONDNOWAVENOCH(:,:,i))>0;
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