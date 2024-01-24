clear;close all;clc

N=120;
M=120;

A=ones(N,M);h=A*0;

%A(3:30,3:60)=0;
dx=2;

r=25;
for i=-r:r
    for j=-r:r
        if i^2+j^2<=r^2
            h(N/2+i,M/2+j)=1;
        end
    end
end

% dx=2;
% h(1:50,50:55)=1;
% h(50:55,1:55)=1;
% h(1:50,50)=1;
% h(50,1:50)=1;

% N=N*4;M=M*4;
% A=ones(N,M);h=A*0;
% dx=2;
% h(1:50*4,50*4:55*4)=1;
% h(50*4:55*4,1:55*4)=1;

%h=h+0.00001*(rand(N,M)-0.5);
%h(:,1:N/3)=1;
    

periodic=0;
alpha=100;
ho=h;
h=smoothEcurvature(h,A,dx,periodic,alpha);
[CURV,HCURV,Cmax,Cmin,god] = surfatureGiulio(A,h,dx);
%X=[1:N]'*ones(M,1)'*dx;Y=ones(N,1)*[1:M]*dx;
%[CURV,HCURV,Cmax,Cmin,Pmax1,Pmin1] = surfature(X,Y,h);
%[CURV,HCURV,Cmax,Cmin,Pmax1,Pmin1] = surfatureNEW(X,Y,h,dx);

%CURV=CURV.*(h>0);
%CURV=smoothEcurvature(CURV,A,dx,periodic,alpha);
%HCURV=HCURV.*(h>0);
% Cmax=Cmax.*(h>0);
% Cmin=Cmin.*(h>0);

%Cmin=smoothEcurvature(Cmin,A,dx,periodic,alpha);
%Cmax=smoothEcurvature(Cmax,A,dx,periodic,alpha);

%cax=HCURV.*CURV;
CURV=sign(CURV).*sqrt(abs(CURV)).*ho*100;

%cax=(2*cax+circshift(cax,[0 1])+circshift(cax,[0 -1])+circshift(cax,[1 0])+circshift(cax,[-1 0]))/4;

%Cmin(h==0)=0.3;
figure
subplot(2,2,1);imagesc(CURV);colormap('jet');colorbar;%caxis([-0.015 0.015])
subplot(2,2,2);imagesc(HCURV);colorbar
subplot(2,2,3);imagesc(ho);colorbar
subplot(2,2,4);imagesc(god);colorbar
%subplot(2,2,4);imagesc(Cmin);colorbar% *(1+alpha/dx^2)


% figure
% imagesc(CURV);colormap('jet');colorbar
sum(CURV(:))

% figure
% subplot(2,2,1);imagesc((HCURV.^2 - CURV));colormap('jet');colorbar
% subplot(2,2,2);imagesc((Cmax+Cmin)/2);colorbar
% subplot(2,2,3);imagesc(Pmin1);colorbar
% subplot(2,2,4);imagesc(Pmax1);colorbar
%figure
%imagesc(CURV*(1+alpha/dx^2)^2-Cmin*(1+alpha/dx^2).*Cmax*(1+alpha/dx^2))

% X=[1:N]'*ones(M,1)'*dx;
% Y=ones(N,1)*[1:M]*dx;
% [K,H,Pmax,Pmin] = surfature(X,Y,h);
% 
% figure;imagesc(Cmin);colormap('jet');colorbar
%sum(cax(:))

%P=2*pi*r

%(sum(CURV(:)>0)-sum(CURV(:)<0))/P
%sum(HCURV(:))/P