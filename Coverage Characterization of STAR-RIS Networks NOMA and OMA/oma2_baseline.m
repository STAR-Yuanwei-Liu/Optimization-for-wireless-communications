function [g,h]=oma2_baseline(x,miu,gamma1,gamma2)
%r1=r(1),D2=r(2),p1=x(3),p2=x(4),beta1=x(5),beta2=x(6),r0=x(7),w1=x(8),w(2)=x(9)
%miu2=0.5;
Pmax=1;
% gamma1=4;
% gamma2=6;
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
g=[x(9)*(2^(gamma2/x(9))-1)/rou0*delta/0.0004*x(2)^2.2-x(4)
x(8)*(2^(gamma1/x(8))-1)/rou0*delta/0.00042*x(1)^2.2-x(3)
x(3)+x(4)-Pmax
miu*x(7)-x(1)
(1-miu)*x(7)-x(2)
x(8)+x(9)-1];
%abs(x(1)^2-x(2)^2)-2;];
h=[x(5)+x(6)-1
  ];