function [zb]=seaboundaryNeumanbedelevationALLBOUNDARY(A,zb)

[N,M]=size(A);
p=find(A==2);%exclude the NOLAND CELLS (A==0)

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N] 

[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))==1);%only inclued the cells in whcih you can creep to



%figure;plot(Y2(A==2));%pause

[px,py]=ind2sub(size(A),p(a));
[qx,qy]=ind2sub(size(A),q(a));

    for i=1:length(a)
    zb(px(i),py(i))=zb(qx(i),qy(i));
    end
end
