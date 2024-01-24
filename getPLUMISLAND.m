function [d,A,M,m,n,dx,do]=getPLUMISLAND

dx=4;%5/2;%0/1;
load PIElarge
d=1.35-PIElarge(1:2:end,1:2:end);
d=double(d);
d(end-1:end,:)=NaN;
d(:,1:40)=NaN;
d(1:2,:)=NaN;
[m,n]=size(d);
A=ones(m,n);
M=A*0;
A(1:end,end)=2;
A(isnan(d))=0;
d(A==0)=NaN;
do=d;