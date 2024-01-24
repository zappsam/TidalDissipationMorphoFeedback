clear;close all
load zBarnsmall
Z=-zBarnsmall;
[N,M]=size(Z);
% DEM = GRIDobj([1:M]',[1:N],Z)
%  %DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%       DEMf = fillsinks(DEM);
% 
%       
%       toc
%       % display sink depth
%       DIFFDEM = DEMf-DEM;
% %       DIFFDEM.Z(DIFFDEM.Z==0) = nan;
% %       imageschs(DEM,DIFFDEM.Z);
%       
%       figure;imagesc(DIFFDEM.Z)
%       %figure;imagesc(-DEMf)
%       
%       
      ponddepth=0.1;
      ZZ=Z;ZZ(isnan(ZZ))=0;
      ZZ=imfill(ZZ);
      DIF = ZZ-Z;
      S=DIF>ponddepth;
      figure;
      subplot(2,1,1);imagesc(ZZ);caxis([1 1.5])
      subplot(2,1,2);imagesc(S);%caxis([1 1.5])