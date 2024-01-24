function [deltaY1,deltaY2,Pedge,Y2OX]=Edgeerosion(P,z,aw,fox,dt,dx,MASK,A,fracY1,fracY2,Y2OX);

%MASK 1 are the mudflat
%MASK 2 is where the are no waves
%fox=0;

%the non-periodic cell are walls
MASK(1,:)=0;
MASK(end,:)=0;

%you get the wave power from the 4 sourronding cells
Pedge=[P(:,1)*0 P(:,1:end-1)]+[P(1,:)*0; P(1:end-1,:)]+[P(:,2:end) P(:,end)*0]+[P(2:end,:); P(end,:)*0];
Pedge(MASK==1)=0;%the wave power in the mudflat becomes zero!
%find the marsh cells with some wave power around it
edg=find(A==1 & MASK==0 & Pedge>0);
%rng('shuffle');
r=rand(length(edg),1);% you only need rng at the beginnign of the loop
a=find(r<aw*Pedge(edg)*dt/dx);


%these cells will erode (the "high" cells)
edger=edg(a);

deltaY1=MASK*0;
deltaY2=MASK*0;

[N,M]=size(z);
for i=1:length(edger)
    
    %these are the adjacent cells
    [I,J] = ind2sub(size(z),edger(i));
    p=[sub2ind(size(z),mod(I+1-1,N)+1,J) sub2ind(size(z),mod(I-1-1,N)+1,J) sub2ind(size(z),I,mod(J+1-1,M)+1) sub2ind(size(z),I,mod(J-1-1,M)+1)];
    a=find(MASK(p)==1 & A(p)==1); %only choses the adjacent cells if they are mudlfat and if they are standard cells
 
 %pause(1)
    if a>0 %yes, there are adjacent cells   
            dz=([z(p(a))-z(edger(i))]);
            dz=mean(max(dz,0));

            %erosion of the boundary cell
            deltaY1(edger(i))=deltaY1(edger(i))+dz.*fracY1(edger(i));
            deltaY2(edger(i))=deltaY2(edger(i))+dz.*fracY2(edger(i));
    
            %redistirbution of the eroded sediment
            totcell=1+length(a);
    
            deltaY1(edger(i))=deltaY1(edger(i))-dz/totcell.*fracY1(edger(i));
            deltaY1(p(a))=deltaY1(p(a))-dz/totcell.*fracY1(edger(i));
    
            Y2OX=Y2OX+fox*dz.*fracY2(edger(i));
            deltaY2(edger(i))=deltaY2(edger(i))-dz/totcell.*fracY2(edger(i))*(1-fox);
            deltaY2(p(a))=deltaY2(p(a))-dz/totcell.*fracY2(edger(i))*(1-fox);
    end
   
end

       %ggg=find(A(p)==10);
    %if ggg>0 pause;end
    
%     g=find(isnan(dz));if g>0;
%         p
%         a
%         p(a)
%         A(p(a))
%         pause;end
    