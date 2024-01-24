clear; close all;clc
%the first point is alwys January 1st at 2400. the second point is jan 1at
%at 1 am
% 
% 


ko=0.002;fGOOS=log(10/ko)/log(3/ko);


a=dlmread('GOOSES2010.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2010,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2010,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS10 = struct('w',w,'d',d);


a=dlmread('GOOSES2011.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2011,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2011,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS11 = struct('w',w,'d',d);

a=dlmread('GOOSES2012.txt');yr=366;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2012,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2012,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS12 = struct('w',w,'d',d);


a=dlmread('GOOSES2013.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2013,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2013,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS13 = struct('w',w,'d',d);


a=dlmread('GOOSES2014.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2014,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2014,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS14 = struct('w',w,'d',d);


a=dlmread('GOOSES2015.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2015,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2015,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS15 = struct('w',w,'d',d);


a=dlmread('GOOSES2016.txt');yr=366;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2016,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2016,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS16 = struct('w',w,'d',d);

a=dlmread('GOOSES2017.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2017,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2017,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS17 = struct('w',w,'d',d);

a=dlmread('GOOSES2018.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2018,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2018,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS18 = struct('w',w,'d',d);

a=dlmread('GOOSES2019.txt');yr=365;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2019,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2019,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS19 = struct('w',w,'d',d);

a=dlmread('GOOSES2020.txt');yr=366;
w=ones(yr*24*6,1)*NaN;d=ones(yr*24*6,1)*NaN;
t = round((datenum(2020,a(:,2),a(:,3),a(:,4),a(:,5),0) -datenum(2020,1,1,0,0,0))*24*6+1);
w(t)=a(:,7);d(t)=a(:,6);
w=interp1([1:yr*24*6]/6,w,[1:yr*24]');
d=interp1([1:yr*24*6]/6,d,[1:yr*24]');
a=find(w==99);w(a)=NaN;a=find(d==999);d(a)=NaN;
WGOOS20 = struct('w',w,'d',d);


WGOOS=struct('y10',WGOOS10,'y11',WGOOS11,'y12',WGOOS12,'y13',WGOOS13,'y14',WGOOS14,'y15',WGOOS15,'y16',WGOOS16,'y17',WGOOS17,'y18',WGOOS18,'y19',WGOOS19,'y20',WGOOS20);
save WGOOS WGOOS


w=[];d=[];
% w=[w; WGOOS.y10.w*fGOOS];d=[d; WGOOS.y10.d];
% w=[w; WGOOS.y11.w*fGOOS];d=[d; WGOOS.y11.d];
% w=[w; WGOOS.y12.w*fGOOS];d=[d; WGOOS.y12.d];
% w=[w; WGOOS.y13.w*fGOOS];d=[d; WGOOS.y13.d];
% w=[w; WGOOS.y14.w*fGOOS];d=[d; WGOOS.y14.d];
% w=[w; WGOOS.y15.w*fGOOS];d=[d; WGOOS.y15.d];
w=[w; WGOOS.y16.w*fGOOS];d=[d; WGOOS.y16.d];
w=[w; WGOOS.y17.w*fGOOS];d=[d; WGOOS.y17.d];
w=[w; WGOOS.y18.w*fGOOS];d=[d; WGOOS.y18.d];
w=[w; WGOOS.y19.w*fGOOS];d=[d; WGOOS.y19.d];
w=[w; WGOOS.y20.w*fGOOS];d=[d; WGOOS.y20.d];


%w(w==0)=NaN;
weq=(nanmean(w.^2)).^(1/2)

% wind_rose(d+180,w,'ci',[0:4:12],'n',16,'di',[0:2.5:15],'dtype','meteo')