function [K,H,Pmax,Pmin,god] = surfatureGiulio(A,Z,dx),
% SURFATURE -  COMPUTE GAUSSIAN AND MEAN CURVATURES OF A SURFACE
%   [K,H] = SURFATURE(X,Y,Z), WHERE X,Y,Z ARE 2D ARRAYS OF POINTS ON THE
%   SURFACE.  K AND H ARE THE GAUSSIAN AND MEAN CURVATURES, RESPECTIVELY.
%   SURFATURE RETURNS 2 ADDITIONAL ARGUEMENTS,
%   [K,H,Pmax,Pmin] = SURFATURE(...), WHERE Pmax AND Pmin ARE THE MINIMUM
%   AND MAXIMUM CURVATURES AT EACH POINT, RESPECTIVELY.

[N,M] = size(Z);


p = find(A==1 | A==2);
NN=length(p);
G=0*A;G(p)=[1:NN];



% % First Derivatives
% [Zu,Zv] = gradient(Z);
% 
% % Second Derivatives
% [Zuu,Zuv] = gradient(Zu);
% [Zuv,Zvv] = gradient(Zv);
% 
% % Reshape 2D Arrays into Vectors
% Zu = Zu(:); 
% Zv = Zv(:); 
% Zuu = Zuu(:); 
% Zuv = Zuv(:); 
% Zvv = Zvv(:); 


Zu = zeros(NN,1);
Zv = zeros(NN,1);
Zuu = zeros(NN,1); 
Zuv = zeros(NN,1); 
Zvv = zeros(NN,1); 

%First Derivatives
for k = [N -1 1 -N]   
[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
    if k==1 | k==-1; Zv(G(p(a)))=Zv(G(p(a)))-sign(k)*(Z(p(a))-Z(q(a)))/2;
    else             Zu(G(p(a)))=Zu(G(p(a)))-sign(k)*(Z(p(a))-Z(q(a)))/2;
    end
end

%Second Derivatives
for k = [N -1 1 -N]   
[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
    if k==1 | k==-1; Zvv(G(p(a)))=Zvv(G(p(a)))-sign(k)*(Zv(G(p(a)))-Zv(G(q(a))))/2; Zuv(G(p(a)))=Zuv(G(p(a)))-sign(k)*(Zu(G(p(a)))-Zu(G(q(a))))/2;
    else             Zuu(G(p(a)))=Zuu(G(p(a)))-sign(k)*(Zu(G(p(a)))-Zu(G(q(a))))/2; %Z%%vu(p(a))=Zuv(p(a))-sign(k)*(Zv(p(a))-Zv(q(a)))/2;
    end
end



% 
% %Second Derivatives
% for k = [N -1 1 -N]   
% [a,q]=excludeboundarycell(k,N,M,p);
% a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
%     if k==1 | k==-1; Zvv(G(p(a)))=Zvv(G(p(a)))-(Z(p(a))-Z(q(a))); 
%     else             Zuu(G(p(a)))=Zuu(G(p(a)))-(Z(p(a))-Z(q(a))); 
%     end
% end
% 
% %Second derivative cross
% for k = [N+1 N-1 -N+1 -N-1]   
% [a,q]=excludeboundarycellDIAG(k,N,M,p);
% a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
%     if k==N+1;       Zuv(G(p(a)))=Zuv(G(p(a)))-Z(q(a));
%     elseif k==N-1;   Zuv(G(p(a)))=Zuv(G(p(a)))+Z(q(a));
%     elseif k==-N+1;  Zuv(G(p(a)))=Zuv(G(p(a)))+Z(q(a));
%     elseif k==-N-1;  Zuv(G(p(a)))=Zuv(G(p(a)))-Z(q(a));
%     end 
% end
% 




Zu = Zu/dx; 
Zv = Zv/dx; 
Zuu = Zuu/dx^2; 
Zuv = Zuv/dx^2; 
Zvv = Zvv/dx^2; 

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
E           =   dot(XU,XU,2);     %1+Zu.*Zu;%
F           =   dot(XU,XV,2);     %Zu.*Zv;%
G           =   dot(XV,XV,2);     %1+Zv.*Zv;%

m           =   cross(XU,XV,2);
p           =   sqrt(dot(m,m,2));               %p           =   sqrt(1+m(:,1).*m(:,1)+1+m(:,2).*m(:,2));
n           =   m./[p p p]; 

% Second fundamental Coeffecients of the surface (ELLE,EMME,ENNE)
ELLE           =   dot(XUU,n,2);
EMME           =   dot(XUV,n,2);
ENNE           =   dot(XVV,n,2);


% % Gaussian Curvature
% god=1./(E.*G - F.^2);
% K = (ELLE.*ENNE - EMME.^2).*god;
% 
% %Mean Curvature
% H = (E.*ENNE + G.*ELLE - 2.*F.*EMME)/2.*god;



god=1./(E.*G - F.^2);
K = (ELLE.*ENNE - EMME.^2).*god;

% Mean Curvature
H = (E.*ENNE + G.*ELLE - 2.*F.*EMME)./2.*god;



p = find(A==1 |A==2);
Ko=A*0;Ko(p)=K;K=Ko;
Ho=A*0;Ho(p)=H;H=Ho;
godo=A*0;godo(p)=god;god=godo;


%Principal Curvatures
Pmax = H + sqrt(H.^2 - K);
Pmin = H - sqrt(H.^2 - K);


%G=A*0;G(p)=god;
%figure;imagesc(G);colorbar;pause

% p=find(A==0);
% for k = [N -1 1 -N]   
% [a,q]=excludeboundarycell(k,N,M,p);
% a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells
% H(p(a))=0;H(q(a))=0;
% K(p(a))=0;K(q(a))=0;
% Pmax(p(a))=0;Pmax(q(a))=0;
% Pmin(p(a))=0;Pmin(q(a))=0;
% end
% %Second derivative cross
% for k = [N+1 N-1 -N+1 -N-1]   
% [a,q]=excludeboundarycellDIAG(k,N,M,p);
% a=a(A(q(a))==1);%exlcude the translated cell that are NOLAND cells
% H(p(a))=0;H(q(a))=0;
% K(p(a))=0;K(q(a))=0;
% Pmax(p(a))=0;Pmax(q(a))=0;
% Pmin(p(a))=0;Pmin(q(a))=0;
% end









%Kout = reshape(K,N,M);







% p=find(A==2);
% for k = [N -1 1 -N]   
% [a,q]=excludeboundarycell(k,N,M,p);
% a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
% H(p(a))=0;H(q(a))=0;
% K(p(a))=0;K(q(a))=0;
% Pmax(p(a))=0;Pmax(q(a))=0;
% Pmin(p(a))=0;Pmin(q(a))=0;
% end
% 
% p=find(A==0);
% for k = [N -1 1 -N]   
% [a,q]=excludeboundarycell(k,N,M,p);
% a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
% H(p(a))=0;H(q(a))=0;
% K(p(a))=0;K(q(a))=0;
% Pmax(p(a))=0;Pmax(q(a))=0;
% Pmin(p(a))=0;Pmin(q(a))=0;
% end
% 






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % First Derivatives
% [Zu,Zv] = gradient(Z);
% 
% % Second Derivatives
% [Zuu,Zuv] = gradient(Zu);
% [Zuv,Zvv] = gradient(Zv);
% 
% gg=Zu;
% 
% % Reshape 2D Arrays into Vectors
% Zu = Zu(:); 
% Zv = Zv(:); 
% Zuu = Zuu(:); 
% Zuv = Zuv(:); 
% Zvv = Zvv(:); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

