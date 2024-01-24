function [Yb,Y1,Y2,zb,zs,plyr,flyr,flyrb]=initializestratigraphy(z,N,M,P)
Yb=ones(N,M)*P.tcko;%tickness of bed layer
plyr=ones(N,M)*P.levo; %initial level occupied 

zb=z+(Yb+plyr.*(P.dlyr))+P.YUi;%substrate (hard rock) elevation

%the heigth of the statigraphy column (without Y, the mobile layer)
zs=zb-(Yb+plyr.*(P.dlyr)); 

%the hight of the mobile layer
Y=zs-z;

%initial composition of the active layer
Y1=P.initialfU*Y;Y2=(1-P.initialfU)*Y;
%composition of the stratigraphy column
flyr=NaN*zeros(N,M,P.nlyr);
flyr(:,:,1:plyr)=P.initialf;
flyrb=ones(N,M)*P.initialf;%composition of bottom layer