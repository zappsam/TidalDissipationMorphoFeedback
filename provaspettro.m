clear;close all;clc
Tp_swell=8;Ho=1;
omegai=2*pi/Tp_swell*[0.05:0.0005:10];Tperiodi=2*pi./omegai;
%omegai=2*pi/Tp_swell*[0.5 1 2 4 8];Tperiodi=2*pi./omegai;
[Ejonswap, domg] = JONSWAP(omegai , Ho, Tp_swell);

%Ejonswap=Ejonswap.*domg;


%Ejonswap=Ejonswap/sum(Ejonswap.*domg);
sum(Ejonswap.*domg)
figure
plot(omegai/2/pi,Ejonswap,'.-')


%omegai=2*pi/Tp_swell*[0.75 1 1.25 1.5 1.75 2];Tperiodi=2*pi./omegai;
%omegai=2*pi/Tp_swell*[0.05:0.01:10];
omegai=2*pi/Tp_swell*[1 1.5 2 2.5 3];Tperiodi=2*pi./omegai;
[Ejonswap, domg] = JONSWAP(omegai , Ho, Tp_swell);
for i=1:length(omegai)
    a= JONSWAP(omegai(i)+[-0.5:0.01:0.5]*domg(i), Ho, Tp_swell);Ejonswap(i)=mean(a);
end
%Ejonswap=Ejonswap.*domg;
%Ejonswap=sqrt(Ejonswap);
%Ejonswap=Ejonswap/sum(Ejonswap);
%Ejonswap=sqrt(Ejonswap);
sum(Ejonswap)

hold on
plot(omegai/2/pi,Ejonswap,'.-')

%Tmo1=2*pi./sum(omegai.*Ejonswap.*domg)
Tmo1=2*pi./sum(omegai.*Ejonswap)

x=[    0.0500
    0.0603
    0.0727
    0.0877
    0.1057
    0.1275
    0.1538
    0.1854
    0.2236
    0.2696
    0.3252
    0.3921
    0.4729
    0.5702
    0.6877
    0.8293
    1.0000];
y=[ -0.9900E+02
-0.9900E+02
-0.9900E+02
 0.5430E+03
 0.3570E+04
 0.1454E+05
 0.3417E+04
 0.1707E+04
 0.7670E+03
 0.3207E+03
 0.1297E+03
 0.5158E+02
 0.2036E+02
 0.8010E+01
 0.3146E+01
 0.1234E+01
 0.4841E+00
];

y2=[ -0.9900E+02
 -0.9900E+02
 -0.9900E+02
  0.2213E+01
  0.1798E+02
  0.1041E+03
  0.4377E+02
  0.4779E+02
  0.4836E+02
  0.4917E+02
  0.4480E+02
  0.2639E+02
  0.1089E+02
  0.3780E+01
  0.1249E+01
  0.3778E+00
  0.9458E-01
];

hold on;plot(x,y/4000,'.-g',x,y2/400,'.-c')

