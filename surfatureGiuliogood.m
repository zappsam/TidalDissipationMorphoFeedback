function [K] = surfatureGiulio(Z),
% SURFATURE -  COMPUTE GAUSSIAN AND MEAN CURVATURES OF A SURFACE
%   [K,H] = SURFATURE(X,Y,Z), WHERE X,Y,Z ARE 2D ARRAYS OF POINTS ON THE
%   SURFACE.  K AND H ARE THE GAUSSIAN AND MEAN CURVATURES, RESPECTIVELY.
%   SURFATURE RETURNS 2 ADDITIONAL ARGUEMENTS,
%   [K,H,Pmax,Pmin] = SURFATURE(...), WHERE Pmax AND Pmin ARE THE MINIMUM
%   AND MAXIMUM CURVATURES AT EACH POINT, RESPECTIVELY.

[s,t] = size(Z);
NN=s*t;

% First Derivatives
[Zu,Zv] = gradient(Z);

% Second Derivatives
[Zuu,Zuv] = gradient(Zu);
[Zuv,Zvv] = gradient(Zv);

% Reshape 2D Arrays into Vectors
Zu = Zu(:); 
Zv = Zv(:); 
Zuu = Zuu(:); 
Zuv = Zuv(:); 
Zvv = Zvv(:); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xuu=zeros(NN,1);Yuu=zeros(NN,1);Xuv=zeros(NN,1);
Yuv=zeros(NN,1);Xvv=zeros(NN,1);Yvv=zeros(NN,1);

Xu=zeros(NN,1);Xv=ones(NN,1);
Yu=ones(NN,1);Yv=zeros(NN,1);


XU          =   [Xu Yu Zu];
XV          =   [Xv Yv Zv];
XUU         =   [Xuu Yuu Zuu];
XUV         =   [Xuv Yuv Zuv];
XVV         =   [Xvv Yvv Zvv];

% First fundamental Coeffecients of the surface (E,F,G)
E           =   dot(XU,XU,2);%1+Zu.*Zu;%
F           =   dot(XU,XV,2);%Zu.*Zv;%
G           =   dot(XV,XV,2);%1+Zv.*Zv;%

m           =   cross(XU,XV,2);
p           =   sqrt(dot(m,m,2));
%p           =   sqrt(1+m(:,1).*m(:,1)+1+m(:,2).*m(:,2));
n           =   m./[p p p]; 

% Second fundamental Coeffecients of the surface (L,M,N)
L           =   dot(XUU,n,2);
M           =   dot(XUV,n,2);
N           =   dot(XVV,n,2);


% Gaussian Curvature
K = (L.*N - M.^2)./(E.*G - F.^2);
K = reshape(K,s,t);

%size(n)
% figure;
% a = reshape(n(:,3),s,t);
% subplot(2,1,1);imagesc(a)
% 
% a = reshape(Zu.*Zv,s,t);
% subplot(2,1,2);imagesc(a)
% pause

% % Mean Curvature
% H = (E.*N + G.*L - 2.*F.*M)./(2*(E.*G - F.^2));
% H = reshape(H,s,t);
% 
% % Principal Curvatures
% Pmax = H + sqrt(H.^2 - K);
% Pmin = H - sqrt(H.^2 - K);
