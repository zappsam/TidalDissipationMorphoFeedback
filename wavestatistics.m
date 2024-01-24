clear;close all;clc
% a=dlmread('ST63189_v03.onlns');
% H=a(:,10);
% save H H
load H
%H=H(1:1000);

%trheshold
Ho=2;

S=H*0;
a=find(H>Ho);
S(a)=1;
d=S-[S(1); S(1:end-1)];
E=find(d==1);
Ee=find(d==-1)-1;
for i=1:length(E)
    He(i)=max(H(E(i):Ee(i)));
end
figure;
subplot(3,1,1);plot(S,'.')
subplot(3,1,2);plot(d,'.')
subplot(3,1,3);plot(H,'.')

He=He-Ho;

figure
duration=(Ee-E);
lag=[E(2:end); Ee(end)]-Ee;%diff(E);lag=[lag; lag(1)];
subplot(1,3,1);[y,x]=hist(He,50);semilogy(x,y','o')
subplot(1,3,2);[y,x]=hist(lag,50);semilogy(x,y','o')
subplot(1,3,3);[y,x]=hist(duration,50);semilogy(x,y,'o')



pH= expfit(He)
pL= expfit(lag)
pD= expfit(duration)

n=500;
l = exprnd(pL,n,1);
d= exprnd(pD,n,1);
h= exprnd(pH,n,1);

d=h*50;

S=[];
for i=1:n
    dur=floor(d(i));if mod(dur,2)==0;dur=dur+1;end
    v=h(i)*[[0:(dur-1)/2] [(dur-1)/2-1:-1:0]]'/((dur-1)/2);
    S=[S; 0*ones(floor(l(i)),1); v];
 end
S=S+Ho;
figure;
H(H<=Ho)=NaN;
S(S<=Ho)=NaN;
subplot(2,1,1);plot(H(1:length(S)));ylim([0 5]);xlim([0 length(S)])
subplot(2,1,2);plot(S);ylim([0 5]);xlim([0 length(S)])

% figure
% [y,x]=hist(He,50);y=y/sum(y);
% subplot(1,2,1);plot(x,y)
% subplot(1,2,2);plot(x,cumsum(y))
% 
% 
% % [parmhat,parmci] = lognfit(He)
% % X=[0.1:.01:5];
% % Y=lognpdf(x,parmhat(1),parmhat(2));Y=Y*(x(2)-x(1));
% % hold on;
% % %P=k/l*(x/l).^(k-1).*exp(-(x/l).^k)*(x(2)-x(1));
% % plot(x,Y,'.-r')
% % 
% 
% 
% [parmhat, parmci] = wblfit(He)
% l=parmhat(1);
% k=parmhat(2);
% 
% X=[0.1:.01:5];
% Y=wblpdf(x,parmhat(1),parmhat(2));Y=Y*(x(2)-x(1));
% P=k/l*(x/l).^(k-1).*exp(-(x/l).^k)*(x(2)-x(1));
% subplot(1,2,1);hold on;plot(x,Y,'.r',x,P)
% 
% Y=wblcdf(x,parmhat(1),parmhat(2));
% P=1-exp(-(x/l).^k);
% subplot(1,2,2);hold on;plot(x,Y,'.c',x,P)
% 
% yr=35;
% events=2459;
% epyr=events/35
% % 
% % RT=20;
% % Q=1-(1/(events*RT))

