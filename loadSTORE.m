clear;close all;clc
load STORE

figure
a=squeeze(STORE(:,:,95)-STORE(:,:,100));
%subplot(2,2,1)
imagesc(a');
colormap('jet');caxis([-1 1])