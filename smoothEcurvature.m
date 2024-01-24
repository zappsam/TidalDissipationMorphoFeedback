function P=smoothEcurvature(z,A,dx,periodic,alpha,h);



G=0*z;a=find(A~=0);NN=length(a);G(a)=[1:NN];
%rhs=ones(NN,1).*min(DH,max(0,d(a)))/(T/2*3600*24); %in m/s!!!
rhs=z(a); %in m/s!!!

[N,M] = size(G);i=[];j=[];s=[];

S=0*G;
%exclude the NOLAND CELLS (A==0)
p = find(A~=0);[row col]=ind2sub(size(A),p);%
for k = [N -1 1 -N]

%the translated cells
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

%apazzo=a(A(q(a))==0);%exlcude the translated cell that are NOLAND cells

a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells

value=alpha/dx^2+0*z(p(a));
S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell

%value=alpha/dx^2+0*z(p(apazzo));
%S(p(apazzo))=S(p(apazzo))+value; %exit from that cell

end




%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];

ds2 = sparse(i,j,s);%solve the matrix inversion
p=ds2\rhs;
P=G;P(G>0)=full(p(G(G>0)));



