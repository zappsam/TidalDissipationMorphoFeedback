function [dx,N,M,A,z,x,y,SPCLcell]=getPIE;
%geometric parameters

load ZBarn
z=-double(ZBarn);
z=z(1:end-200,1:1500-325);
z(z<-2)=-2;
%z=z*0+3;

dx=2;

figure;imagesc(ZBarn)
pause




[N,M]=size(z);

%%%%%%%%%%%cell types
A=ones(N,M);

A(z==-2)=0;
A(end,:)=0;
%figure;imagesc(A);pause

%river b.c.
rivermouthfront=[];
SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;

for k=1:20
    
    
A1=circshift(A,[0 1]);A2=circshift(A,[0 -1]);A3=circshift(A,[1 0]);A4=circshift(A,[-1 0]);
A1(:,1)=A(:,1);A2(:,end)=A(:,end);A3(1,:)=A(1,:);A4(end,:)=A(end,:);
A=min(A,min(A1,min(A2,min(A3,A4))));
    
CC = bwconncomp(A,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
A(A==1)=0;A(CC.PixelIdxList{idx}) = 1;

% 
% A1=circshift(A,[1 1]);A2=circshift(A,[1 -1]);A3=circshift(A,[-1 1]);A4=circshift(A,[-1 -1]);
% A1(:,1)=A(:,1);A2(:,end)=A(:,end);A3(1,:)=A(1,:);A4(end,:)=A(end,:);
% A=min(A,min(A1,min(A2,min(A3,A4))));

CC = bwconncomp(1-A,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
for i=1:length(CC.PixelIdxList)
    quanti=length(CC.PixelIdxList{i});
    if quanti<10000
    A(CC.PixelIdxList{i}) = 1;    
    end
end

end



%A(end-200:end,end-340:end)=0;
% for i=0:180    
% A(end-160-floor((180-i)/15):end,end-i)=0;
% end
% 
for i=0:200    
A(end-floor((200-i)*1.05):end,end-i)=0;
end
% 
for i=0:190    
A(end-floor((190-i)*1.05):end,end-600+i)=0;
end



for i=0:220    
A(end-floor((220-i)*1.05):end,end-940+i)=0;
end


for i=0:320    
A(1:floor(i*1.4),end-320+i)=0;
end


% GG=A*0;
% GG(:,end)=2;
% A(GG==2 & A==1)=2;


A(480:560,end)=2;

%A(z==-2)=0;
z(A==0)=nan;

%sum(z(:))


%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


A=A(1:1:952,:);
z=z(1:1:952,:);
[N,M]=size(z);
% figure;imagesc(z);
%figure;imagesc(A);pause


A=A(end:-1:1,:);
z=z(end:-1:1,:);


%CORREZIONE FOR MSL
z=z-0.1;

%z=z*0+3;


%increase resolution
% N=N*2;M=M*2;dx=dx/2;
% xn=[0:N-1]*dx/1000;yn=[0:M-1]*dx/1000;
% 
% size(A)
% z=interp2(y',x,z,yn',xn);
% 
% A=ones(N,M);
% 
% A(:,end)=0;
% A(95:115,end)=2;
% 
% A(isnan(z))=0;
% 
% A(95*2:115*2,end)=2;
% 
% x=xn;y=yn;







%figure;imagesc(A);pause


%z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%


% load PIEcut
% d=PIEcut;
% d(144:150,1:17)=NaN;
% d(175:177,39:42)=NaN;
% PIEcutfill=d;
% save PIEcutfill PIEcutfill
% figure;
% ax1 = subplot(1,1,1); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,-d);axis equal;set(IM,'alphadata',~(d==0));set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-3 2.7/2],256); %-3 1
% colormap(ax1,cmp)
% caxis([-3 2.7/2]);