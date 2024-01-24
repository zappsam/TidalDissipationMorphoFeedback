function [E]=totalsedimenterosionMUD(U,fUpeak,MANN,VEG,B,fTide,UR,Uwave_sea,Uwave,Tp_sea,Tp_swell,fMFtide,fMFswell,fMFsea,fMFriver,taucr,me,h,lev);


taucro=U*0+taucr;
taucro(VEG==1)=0.8;
%a=find(lev<-1.55);taucro(a)=taucro(a)+0.5*(-lev(a)-1.55);

%figure;imagesc(taucro);pause
%taucr=0;
%taucr=VEG*0+taucr;
%taucr(VEG==1)=0.5;
%taucro=taucr;

%figure;imagesc(taucr);pause

Utide=U*fUpeak;%./max(0.1,fTide);
%Utide(fTide<1)=Utide(fTide<1)*2;
%tauC=1030*9.81/45^2.*Utide.^2;
tauC=1030*9.81*MANN.^2.*h.^(-1/3).*Utide.^2;

%CD=9.8/45^2;
%tauC=h.^(-1/3).*1030.*CD.*Utide.^2;
%tauC=1030.*CD.*Utide.^2;


%tauC(VEG==1)=tauC(VEG==1).*(facMann*B(VEG==1)).^2;
%tauC(VEG==1)=tauC(VEG==1).*facMann^2;


%figure;imagesc(tauC);pause


tauCRiver=1030*0.04^2*9.81*h.^(-1/3).*UR.^2;

ko=1/1000;
aw=Tp_sea.*Uwave_sea/(2*pi);
fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;
%fw=0.015;
tauWsea=(0.5*1030*fw.*Uwave_sea.^2).*fTide;

ko=1/1000;
aw=Tp_swell.*Uwave/(2*pi);
fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;
%fw=0.015;
tauWswell=(0.5*1030*fw.*Uwave.^2).*fTide;



%NONCONTA!!!
%E=me*max(0,tauC-taucr)./taucro*fMFtide;% +me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea +me*max(0,tauCRiver-taucr)./taucro*fMFriver;
%E=me*(sqrt(1+(tauC/taucr).^2)-1)*fMFtide;


E=me./taucro.*(sqrt(taucr.^2+tauC.^2)-taucr)*fMFtide;


%+me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea +me*max(0,tauCRiver-taucr)./taucro*fMFriver;
%%%
%E=h*NaN;
%Ec=me*max(0,tauC-taucr)./taucro*fMFtide;
%Ew=+me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea;

%E=E*0;
%E(h<=0.1)=0;
% Ec(h<=0.1)=0;
% Ew(h<=0.1)=0;

%Ew=Ew*0;
%Ec=Ec*0;