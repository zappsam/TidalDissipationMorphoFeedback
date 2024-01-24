clear;close all;clc
load STORE_t200_r07_ALL
load STORE_t200_r07_NOPOND
load STORE_t200_r07_NOPONDNOWAVE
load STORE_t200_r07_NOPONDNOWAVENOCH
load STORE_t200_r07_NOWAVE
load STORE_t200_r07_ALLUPLAND

load STORE_t200_r3_ALL
load STORE_t200_r3_ALLscarp2m
load STORE_t200_r3_NOWAVE
load STORE_t200_r3_NOPONDNOWAVE
load STORE_t200_r3_NOPONDNOWAVENOCH
load STORE_t200_r3_NOPOND
load STORE_t200_r3_ALLUPLAND

a=STORE_t200_r3_NOPONDNOWAVE;
STORE_t200_r3_NOPONDNOWAVE=STORE_t200_r3_NOPOND;
STORE_t200_r3_NOPOND=a;


A=squeeze(STORE_t200_r07_NOPONDNOWAVE(:,:,1))>0;
T=sum(A(:)==1)

figure;imagesc(A);


L1=zeros(4,4);
figure
for i=1:16
subplot(4,4,i)
B=squeeze(STORE_t200_r07_ALL(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
%IM=imagesc(G');
%set(IM,'alphadata',~(G'==0));
colormap('jet')
caxis([0 2])   
L1(i)=(T-sum(B(:)==1))/T;
end
L1good=L1;

L2=zeros(4,4);
%figure
for i=1:16
subplot(4,4,i)
B=squeeze(STORE_t200_r07_NOPOND(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
%IM=imagesc(G');
%set(IM,'alphadata',~(G'==0));
colormap('jet')
caxis([0 2])  
L2(i)=(T-sum(B(:)==1))/T; 
end

L3=zeros(4,4);
%figure
for i=1:16
subplot(4,4,i)
B=squeeze(STORE_t200_r07_NOWAVE(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
%IM=imagesc(G');
%set(IM,'alphadata',~(G'==0));
colormap('jet')
caxis([0 2])  
L3(i)=(T-sum(B(:)==1))/T; 
end


L4=zeros(4,4);
%figure
for i=1:16
subplot(4,4,i)
B=squeeze(STORE_t200_r07_NOPONDNOWAVE(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
%IM=imagesc(G');
%set(IM,'alphadata',~(G'==0));
colormap('jet')
caxis([0 2])  
L4(i)=(T-sum(B(:)==1))/T; 
end

L5=zeros(4,4);
%figure
for i=1:16
subplot(4,4,i)
B=squeeze(STORE_t200_r07_NOPONDNOWAVENOCH(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
%IM=imagesc(G');
%set(IM,'alphadata',~(G'==0));
colormap('jet')
caxis([0 2])  
L5(i)=(T-sum(B(:)==1))/T; 
end


R=[2.5 5 7.5 10];
for i=1:4
subplot(4,2,i*2-1)
plot(R,L1(i,:),R,L2(i,:),R,L3(i,:),R,L4(i,:),R,L5(i,:))
ylim([0 1]);xlim([2.5 10])
end
legend('ponds+waves','no pond','no wave','no pond no wave')
xlabel('RSLR')




L2=L2-L4;
L3=L3-L4;
SUM=L2+L3+L4;

L2=L2./SUM.*L1;
L3=L3./SUM.*L1;
L4=L4./SUM.*L1;
L5=L5./SUM.*L1;

for i=1:4
subplot(4,2,i*2-1)
%plot(R,L1(i,:),R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
X=[R; R; R; R];
Y=[L5(i,:); L4(i,:)-L5(i,:); L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
area(X',Y');
ylim([0 1]);xlim([2.5 10])
end
legend('ALL','ALL','wave+drowning+widening','drowning+widening','summa')
xlabel('RSLR')















%figure
A=squeeze(STORE_t200_r3_ALL(:,:,1))>0;
T=sum(A(:)==1);

L1=zeros(4,4);
L2=zeros(4,4);
L3=zeros(4,4);
L4=zeros(4,4);
L5=zeros(4,4);

for i=1:16
B=squeeze(STORE_t200_r3_ALL(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
L1(i)=(T-sum(B(:)==1))/T;
end
L1r3good=L1;

L2=zeros(4,4);
for i=1:16
B=squeeze(STORE_t200_r3_NOPOND(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
L2(i)=(T-sum(B(:)==1))/T; 
end

L3=zeros(4,4);
for i=1:16
B=squeeze(STORE_t200_r3_NOWAVE(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
L3(i)=(T-sum(B(:)==1))/T; 
end

L4=zeros(4,4);
for i=1:16
B=squeeze(STORE_t200_r3_NOPONDNOWAVE(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
L4(i)=(T-sum(B(:)==1))/T; 
end

L5=zeros(4,4);
for i=1:16
B=squeeze(STORE_t200_r3_NOPONDNOWAVENOCH(:,:,i))>0;
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
L5(i)=(T-sum(B(:)==1))/T; 
end


% R=[1:4];
% %figure
% for i=1:4
% subplot(4,2,i*2-1)
% plot(R,L1(i,:),R,L2(i,:),R,L3(i,:),R,L4(i,:))
% ylim([0 1]);xlim([2.5 10])
% end
% legend('ponds+waves','no pond','no wave','no pond no wave')
% xlabel('RSLR')




L2=L2-L4;
L3=L3-L4;
SUM=L2+L3+L4;

L2=L2./SUM.*L1;
L3=L3./SUM.*L1;
L4=L4./SUM.*L1;
L5=L5./SUM.*L1;

for i=1:4
subplot(4,2,i*2)
%plot(R,L1(i,:),R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
X=[R; R; R; R];
Y=[L5(i,:); L4(i,:)-L5(i,:); L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
area(X',Y');
ylim([0 1]);xlim([2.5 10])
end
legend('ALL','ALL','wave+drowning+widening','drowning+widening','summa')
xlabel('RSLR')


% figure
% for i=1:4
% subplot(4,2,i*2)
% %plot(R,L1(i,:),R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
% X=[R; R; R; R];
% Y=[L5(i,:); L4(i,:)-L5(i,:); L2(i,:); L3(i,:)];%%,R,L4(i,:)+L3(i,:)+L2(i,:),R,L4(i,:)+L2(i,:),R,L4(i,:),R,L2(i,:)+L3(i,:)+L4(i,:),'--k',R,L5(i,:))
% area(X',Y');
% ylim([0 1]);xlim([2.5 10])
% end






















figure
A=squeeze(STORE_t200_r07_ALLUPLAND(:,:,1)>0 & STORE_t200_r07_ALLUPLAND(:,:,1)<0.35);
T=sum(A(:)==1)

L1=zeros(4,4);
for i=1:16
B=squeeze(STORE_t200_r07_ALLUPLAND(:,:,i)>0 & STORE_t200_r07_ALLUPLAND(:,:,i)<0.35);
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
colormap('jet')
caxis([0 2])   
L1(i)=(T-sum(B(:)==1))/T;
end

R=[1:4];
for i=1:4
subplot(4,2,i*2-1)
plot(R,L1(i,:),R,L1good(i,:),'g')
hold on;plot([R(1) R(end)],[0 0],'--k')
ylim([-0.25 1])
end
legend('ponds+waves','waves','standard')
xlabel('RSLR')


L1=zeros(4,4);
for i=1:16
B=squeeze(STORE_t200_r3_ALLUPLAND(:,:,i)>0 & STORE_t200_r07_ALLUPLAND(:,:,i)<0.35);
I=(A-B);
G=A*0+1;
G(A==0)=0;
G(I==1)=2;
colormap('jet')
caxis([0 2])   
L1(i)=(T-sum(B(:)==1))/T;
end

R=[1:4];
for i=1:4
subplot(4,2,i*2)
plot(R,L1(i,:),R,L1r3good(i,:),'g')
hold on;plot([R(1) R(end)],[0 0],'--k')
ylim([-0.25 1])
end
legend('ponds+waves','waves','standard')
xlabel('RSLR')