function [wsB Cv]=getvegetationWSandCd(B,Bmax,h);
hs=0.0609*(B*Bmax).^0.1876;
ns=250*(B*Bmax).^0.3032;
ds=0.0006*(B*Bmax).^0.3;
as=0.25*(B*Bmax).^0.5;
Cdoveg=0.1;Cv=1/2*as.*h.*Cdoveg;%vegetation drag coefficent

Uref=0.1;
neta=0.224*(Uref.*ds/10^-6).^0.718.*(100*10^-6./ds).^2.08;
wsB=Uref.*neta.*ds.*ns.*min(hs,h);wsB(B<=0)=0;