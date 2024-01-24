function [Um,TP,HS,F,kwave,PW]=SeaWaves(h,angle,hwSea_lim,range,wind,MASK,ndir,dx);
Fetchlim=600;%dx*2*10;
%dlo=hwSea_lim; %minimum water depth to calculate wave. below this you don't calculate it

F=calculatefetch(MASK,ndir,dx,angle);
Fo=F;


%diffuse the fetch field    
alphadiffusefetch=20;%%%QUESTO ERA 1 FINO AD APRILE 23 2018!!!!!
F=diffusefetch(MASK,F,alphadiffusefetch,dx); %F(F<=Fetchlim)=0;
F(Fo<=Fetchlim)=0;

Um=0*h;TP=0*h;HS=0*h;

a=find(Fo>Fetchlim);%h>dlo & 
%a=find(h>dlo);
D=h(a);Ff=F(a);

%[Hs,Tp]=YeV_correction(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
[Hs,Tp]=YeV(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
Tp(Tp==0)=1;

TP(a)=Tp;
HS(a)=Hs;

%do not diffuse in cells outside the MASK
%HS=diffusefetch(MASK,HS,alpha,dx);
%TP=diffusefetch(MASK,TP,alpha,dx);


D=h(a);
kwave=0*h;
kk=wavek(1./TP(a),D);%kk=wavekFAST(1./Tp,D);
kwave(a)=kk;
kwave(kwave==0)=1;

Um=pi*HS./(TP.*sinh(kwave.*h));

cg=(2*pi./kwave./TP)*0.5.*(1+2*kwave.*h./(sinh(2*kwave.*h)));
PW=cg*1030*9.8.*HS.^2/16;

Um(MASK==0)=0;
PW(MASK==0)=0;



%figure;imagesc(Um)
%pause