function wind_rose(wind_direction,wind_speed)
%WIND_ROSE Plot a wind rose
%   this plots a wind rose

%Al Mac (2022). wind_rose(wind_direction,wind_speed) (https://www.mathworks.com/matlabcentral/fileexchange/65174-wind_rose-wind_direction-wind_speed), MATLAB Central File Exchange. Retrieved February 8, 2022.


figure
pax = polaraxes;
polarhistogram(deg2rad(wind_direction(wind_speed<25)),deg2rad(0:10:360),'displayname','20 - 25 m/s','Normalization','probability')
hold on
polarhistogram(deg2rad(wind_direction(wind_speed<20)),deg2rad(0:10:360),'FaceColor','red','displayname','15 - 20 m/s','Normalization','probability')
polarhistogram(deg2rad(wind_direction(wind_speed<15)),deg2rad(0:10:360),'FaceColor','yellow','displayname','10 - 15 m/s','Normalization','probability')
polarhistogram(deg2rad(wind_direction(wind_speed<10)),deg2rad(0:10:360),'FaceColor','green','displayname','5 - 10 m/s','Normalization','probability')
polarhistogram(deg2rad(wind_direction(wind_speed<5)),deg2rad(0:10:360),'FaceColor','blue','displayname','0 - 5 m/s','Normalization','probability')
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
legend('Show')
%title('Wind Rose')
end
