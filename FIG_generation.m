
%load tst186v2_FIGoutputs.mat
TRANGEPLOT=(0.6+0.1)*2;
zbed=STORE(:,:,150);
figure
ax1 = subplot(1,1,1);
IM=imagesc(x,y,zbed');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 TRANGEPLOT/2]-0,256);
colormap(ax1,cmp)
caxis([-3 TRANGEPLOT/2]);
colorbar


figure
linecolors=copper(4);

plot(x,flip(squeeze(TrangePLOT(:,1,1))),'LineWidth',2)%,'color',linecolors(1,:))
hold on
plot(x,flip(squeeze(TrangePLOT(:,1,25))),'LineWidth',2)%,'color',linecolors(2,:))
plot(x,flip(squeeze(TrangePLOT(:,1,50))),'LineWidth',2)%,'color',linecolors(1,:))
plot(x,flip(squeeze(TrangePLOT(:,1,75))),'LineWidth',2)%,'color',linecolors(2,:))
plot(x,flip(squeeze(TrangePLOT(:,1,100))),'LineWidth',2)%,'color',linecolors(3,:))
plot(x,flip(squeeze(TrangePLOT(:,1,125))),'LineWidth',2)%,'color',linecolors(4,:))
plot(x,flip(squeeze(TrangePLOT(:,1,150))),'LineWidth',2)%,'color',linecolors(7,:))
hold off
xlim([0 17])
legend('t=1','t=500','t=1000','t=1500','t=2000','t=2500','t=3000')
%legend('t=1','t=1000','t=2000','t=3000')

h=figure
set(gcf,'Units', 'Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
fontsize=24;

subplot(1,3,1)
%time=20:20:3000;
plot(TOTMARSH,time,'LineWidth',2)
ylabel('years')
xlabel('total marsh area (km2)')
xlim([15 40])
axis ij
subplot(1,3,2)
plot(UVVRPLOT,time,'LineWidth',2)
xlabel('UVVR (-)')
xlim([0 1])
xticks([0 0.2 0.4 0.6 0.8 1])
axis ij
subplot(1,3,3)
plot(squeeze(TrangePLOT(600,1,:)),time,'LineWidth',2)
xlabel('tidal range at 10 km from seaward boundary (m)')
xlim([0 3])
%xticks([0 0.2 0.4 0.6 0.8 1])
axis ij


figure
plot(time,squeeze(TrangePLOT(1100,1,:)))
hold on
plot(time,squeeze(TrangePLOT(1000,1,:)))
plot(time,squeeze(TrangePLOT(900,1,:)))
plot(time,squeeze(TrangePLOT(800,1,:)))
plot(time,squeeze(TrangePLOT(700,1,:)))
plot(time,squeeze(TrangePLOT(600,1,:)))%
plot(time,squeeze(TrangePLOT(500,1,:)))
plot(time,squeeze(TrangePLOT(400,1,:)))
plot(time,squeeze(TrangePLOT(300,1,:)))
hold off
legend('boundary','0','2','4','6','8','10','12','14')