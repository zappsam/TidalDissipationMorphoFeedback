function F=diffusefetch(A,F,alpha,dx,URx,URy,MASK);

fac=0.001;

MASKI=double(MASK);
MASKI(MASKI==0)=fac;
%figure;imagesc(MASK');pause


[N,M]=size(A);
%D=A.*(alpha*3600*24)/(dx^2);

G=0*F;
p=find(A==1);%exclude the NOLAND CELLS (A==0)
NN=length(p);G(p)=[1:NN];rhs=F(p);[m,n]=size(G);i=[];j=[];s=[];

%figure;imagesc(A');pause

S=0*G;
[row col]=ind2sub(size(A),p);
for k = [m -1 1 -m]
%avoid to the the cells out of the domain (risk to make it periodic...)
if k==m;aa=find(col+1<=n);end;if k==-m;aa=find(col-1>0);end;if k==-1;aa=find(row-1>0);end;if k==1;aa=find(row+1<=m);end;

q=p+k;%the translated cell
a=aa(A(q(aa))==1 );%only inclued the cells in whcih you can creep to



value=0*A(p(a));%min(D(p(a)),D(q(a)));%value=(D(p(a))+D(q(a)))/2.*facNL;

% if (k==N);UR=URy(p(a));up=find(UR>0);Y=UR(up);end %East-west
% if (k==-N);UR=URy(p(a));up=find(UR<0);Y=-UR(up);end
% if (k==1);UR=URx(p(a));up=find(UR>0);Y=UR(up);end  %North-south
% if (k==-1);UR=URx(p(a));up=find(UR<0);Y=-UR(up);end
% value(up)=value(up)+alpha*Y/dx;

if (k==N);UR=URy(p(a)).*min(MASKI(p(a)),MASKI(q(a)));up=find(UR>0);Y=UR(up);end %East-west
if (k==-N);UR=URy(p(a)).*min(MASKI(p(a)),MASKI(q(a)));up=find(UR<0);Y=-UR(up);end
if (k==1);UR=URx(p(a)).*min(MASKI(p(a)),MASKI(q(a)));up=find(UR>0);Y=UR(up);end  %North-south
if (k==-1);UR=URx(p(a)).*min(MASKI(p(a)),MASKI(q(a)));up=find(UR<0);Y=-UR(up);end
value(up)=value(up)+alpha*Y/dx;

% if (abs(k)==N);UR=(URy(p(a))+URy(p(a)))/2;Y=UR;end %East-west
% if (abs(k)==1);UR=(URx(p(a))+URx(p(a)))/2;Y=UR;end  %North-south
% value=value+alpha*Y/dx;


S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell

% Fin=A(p(a))*0;Fin(A(p(a))==1)=1; % (A(q(a))==1 | A(q(a))==10)=1; to conserve the mass a the river mouth = no input
% Fout=A(p(a))*0;Fout(A(q(a))==1)=1; %needed not to affect the b.c. -> Do not choose 2 and 1p
% S(p(a))=S(p(a))+value.*Fin; %exit from that cell
% i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value.*Fout]; %gain from the neigborh cell

end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;F(G>0)=full(P(G(G>0)));

F(MASK==0)=F(MASK==0)/fac;


