function [U,Ux,Uy,q,P]=flowBasin(A,n_manning,h,dx,Q,Umouth,Uo);

Q=abs(Q);
%A(A==11)=1; %the river front cells are normal cells
%to impose no cross flux in inlet
imposenocrossflux=1;

%Uo=0.1;
manning=0*A+n_manning;
A(A==22)=1;  %this behaves as normal flow %but do not update A!
%consider the pond cells as A==0
A(A==3)=1; %the isoalted pond behaves as normal cell (btu different depth...) %but do not update A!

csi=h.^(1/3)./manning.^2./Uo;
%csi=1./manning.^2./Uo*24*3600;
D=csi.*h.^2/(dx^2);%.*eta;
G=0*A;a=find(A~=0);NN=length(a);G(a)=[1:NN];
rhs=zeros(NN,1);  
rhs(G(A==10))=Q/dx; %in m/s !!!

%A(A==10)=1;

[N,M] = size(G);i=[];j=[];s=[];

%boundary conditions imposed water level
a=find(A==2);
i=[i;G(a)]; j=[j;G(a)]; s=[s;ones(size(a))];rhs(G(a))=0;%water level zero

S=0*G;
%exclude the NOLAND CELLS (A==0)
p = find(A==1 | A==10);[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]

%avoid to the the cells out of the domain (risk to make it periodic...)
if k==N;a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells

DD=(D(p(a))+D(q(a)))/2;
%DD=max(D(p(a)),D(q(a)));
if imposenocrossflux==1;if (k==N | k==-N);cross=find(A(p(a))==10);DD(cross)=0;end;end

S(p(a))=S(p(a))+DD; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-DD]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)];

ds2 = sparse(i,j,s);%solve the matrix inversion

p=ds2\rhs;
P=G;P(G>0)=full(p(G(G>0)));
P(A==2 | A==21)=0;  %need when swtinching q and p
D=D./h*dx; %you need to put it later to have the rigth velocity on the cell

U1=0*A;Um1=0*A;UN=0*A;UmN=0*A;
Ux=0*A;Uy=0*A;

p = find(A==1 | A==10 | A==2);[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]
if k==N; a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1; a=find(row+1<=N);end;
q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
DD=(D(p(a))+D(q(a)))/2;
%DD=min(D(p(a)),D(q(a)));

if imposenocrossflux==1;if (k==N | k==-N);cross=find(A(p(a))==10);DD(cross)=0;end;end

% if (k==1 | k==-1);    Uy(p(a))=Uy(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD/2;%./h(p(a))*dx;
% elseif (k==N | k==-N);Ux(p(a))=Ux(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD/2;%./h(p(a))*dx;
% end
if (k==1); U1(p(a))=U1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;%./h(p(a))*dx;
elseif (k==-1); Um1(p(a))=Um1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;%./h(p(a))*dx;
elseif (k==N); UN(p(a))=UN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;%./h(p(a))*dx;
elseif (k==-N); UmN(p(a))=UmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;%./h(p(a))*dx;
end


end
%Ux(A==10)=2*Ux(A==10);Uy(A==10)=2*Uy(A==10);%because I had the factor 2

%Uy=max(abs(U1),abs(Um1));
%Ux=max(abs(UN),abs(UmN));
%Uy=min(abs(U1),abs(Um1));
%Ux=min(abs(UN),abs(UmN));

%ONLY THIS WORKS FOR THE RIVER!!!!
% Uy=(U1+Um1)/2; 
% Ux=(UN+UmN)/2;
% Ux(A==10)=Ux(A==10)*2;
% Uy(A==10)=Uy(A==10)*2;

Uy=max(abs(U1),abs(Um1)).*sign(U1+Um1);
Ux=max(abs(UN),abs(UmN)).*sign(UN+UmN);


%Uy(A==10)=max(abs(U1(A==10)),abs(Um1(A==10))).*sign(U1(A==10)+Um1(A==10));
%Ux(A==10)=max(abs(UN(A==10)),abs(UmN(A==10))).*sign(UN(A==10)+UmN(A==10));


%Ux(A==10)=Umouth;

Uy(A==10)=Q./h(A==10);  %this will impose more uniform discharge over the oulet
%Ux(A==10)=0;

U=sqrt(Ux.^2+Uy.^2);

q=U.*h;



% P1=[P(:,2:end) P(:,end)];P2=[P(:,1) P(:,1:end-1) ];P3=[P(2:end,:); P(end,:)];P4=[P(1,:); P(1:end-1,:)];
% D1=[D(:,2:end) D(:,end)];D2=[D(:,1) D(:,1:end-1) ];D3=[D(2:end,:); D(end,:)];D4=[D(1,:); D(1:end-1,:)];
% A1=[A(:,2:end) A(:,end)];A2=[A(:,1) A(:,1:end-1) ];A3=[A(2:end,:); A(end,:)];A4=[A(1,:); A(1:end-1,:)];
% 
% pp1=P(A1==0);pp2=P(A2==0);pp3=P(A3==0);pp4=P(A4==0);
% P1(A1==0)=pp1;P2(A2==0)=pp2;P3(A3==0)=pp3;P4(A4==0)=pp4;
% 
% pp1=D(A1==0);pp2=D(A2==0);pp3=D(A3==0);pp4=D(A4==0);
% D1(A1==0)=pp1;D2(A2==0)=pp2;D3(A3==0)=pp3;D4(A4==0)=pp4;
% 
% DD1=(D+D1)/2;DD2=(D+D2)/2;DD3=(D+D3)/2;DD4=(D+D4)/2;
% 
%   Ux=0.5*((P-P1).*DD1 + (P2-P).*DD2)./h*dx;
%   Uy=0.5*((P-P3).*DD3 + (P4-P).*DD4)./h*dx;
%   Ux(A==2 | A==10)=2*Ux(A==2 | A==10);
%   Uy(A==2 | A==10)=2*Uy(A==2 | A==10);
% %G=cat(3,(P-P1).*DD1,(P2-P).*DD2);Ux=nanmean(G,3);
% %G=cat(3,(P-P3).*DD3,(P4-P).*DD4);Uy=nanmean(G,3);



