function [Um,TP,HS,F,kwave,PW]=SeaWaves(h,angle,hwSea_lim,range,wind,MASK,ndir,dx);

Um=0*h;TP=0*h;HS=0*h;





Fetchlim=max(50,dx*2);%dx*2;%600;%dx*2*10;
dlo=hwSea_lim; %minimum water depth to calculate wave. below this you don't calculate it

%The standard way
%F=calculatefetch(MASK,ndir,dx,angle);
% 
% %For Georgia
extrafetch=10000;%[m}
%extrafetch=0;%[m}
F=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,angle,extrafetch);
% F1=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,angle,extrafetch);
% F2=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,mod(angle+90,360),extrafetch);
% F3=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,mod(angle+180,360),extrafetch);
% F4=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,mod(angle+270,360),extrafetch);
% F=(F1+F2+F3+F4)/4;

%For the idealize basin
%extrafetch=10000;%[m}
%F=calculatefetchWITHEXTRAS(MASK,ndir,dx,angle,extrafetch);
% F1=calculatefetchWITHEXTRAS(MASK,ndir,dx,angle,extrafetch);
% F2=calculatefetchWITHEXTRAS(MASK,ndir,dx,mod(angle+90,360),extrafetch);
% F3=calculatefetchWITHEXTRAS(MASK,ndir,dx,mod(angle+180,360),extrafetch);
% F4=calculatefetchWITHEXTRAS(MASK,ndir,dx,mod(angle+270,360),extrafetch);
% F=(F1+F2+F3+F4)/4;

Fo=F;
            
%             %For all the modified ways. Creates a buffer on the side
%             %boundaries. Just used as a mask, the actual value is not
%             %importnat, just need to be larger than fetchlim.
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*
%             [N,M]=size(h);
%             Fo(2+floor(N*0.5):end-1,1:20)=9999;
%             Fo(2+floor(N*0.5):end-1,end-20:end)=9999;
%             %%%%%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*


F(Fo<=Fetchlim)=0;

%usa questo per isolared la mudflat
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASK(end-100:end,:)=1;
% F(end-100:end,:)=extrafetch;
% Fo(end-100:end,:)=extrafetch;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MASK(end-1:end,:)=1;
F(end-1:end,:)=extrafetch;
Fo(end-1:end,:)=extrafetch;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%diffuse the fetch field    
alphadiffusefetch=1;   %messo 10 for the VCR wave validation 10;%0;%%%QUESTO ERA 1 FINO AD APRILE 23 2018!!!!!
F=diffusefetch(MASK,F,alphadiffusefetch,dx); 
F(Fo<=Fetchlim | MASK==0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=find(Fo>Fetchlim & h>dlo & F>0 & MASK==1);%h>dlo & %a=find(Fo>dx*2);%h>dlo & %a=find(h>dlo);
D=h(a);Ff=F(a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%TRCUCCOZO TO AVOID depths too small
% hbedsheatresslim=0.5;
% h(h<hbedsheatresslim)=hbedsheatresslim;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[Hs,Tp]=YeV_correction(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
[Hs,Tp]=YeV(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
HS(a)=Hs;
TP(a)=Tp;TP(TP==0)=1;


%do not diffuse in cells outside the MASK
%HS=diffusefetch(MASK,HS,alpha,dx);
%TP=diffusefetch(MASK,TP,alpha,dx);


%hlimbedshearstress=0.5;
%h=max(hlimbedshearstress,h);% to reduce the bed shear stress for very small water depth


kwave=0*h;
kk=wavek(1./TP(a),h(a));%kk=wavekFAST(1./Tp,D);
kwave(a)=kk;
kwave(kwave==0)=1;

Um=pi*HS./(TP.*sinh(kwave.*h));

cg=(2*pi./kwave./TP)*0.5.*(1+2*kwave.*h./(sinh(2*kwave.*h)));
PW=cg*1030*9.8.*HS.^2/16;

Um(MASK==0)=0;
PW(MASK==0)=0;




%figure;imagesc(Um)
%pause