load('TidetimeseriesWIND.mat')

wcorr=45;%angle you want to rotate counterclockwise
newangle=(TidetimeseriesWIND(5,:)-wcorr)';
newangle(newangle<0)=360+newangle(newangle<0);
P.angleWIND=newangle;


figure
subplot(2,2,1)
histogram(TidetimeseriesWIND(1,:),40,'Normalization','probability')
xlabel('tidal range (m)')
% ylim('probability')
subplot(2,2,2)
histogram(TidetimeseriesWIND(2,:),40,'Normalization','probability')
xlabel('tidal period (days)')
% ylim('probability')
subplot(2,2,3)
histogram(TidetimeseriesWIND(3,:),40,'Normalization','probability')
xlabel('MSL anomaly (m)')
% ylim('probability')
% subplot(2,2,3)
% histogram(TidetimeseriesWIND(4,:),'Normalization','probability')
% xlabel('wind velocity (m/s)')

% subplot(2,2,4)
wind_rose_plot(newangle,TidetimeseriesWIND(4,:))