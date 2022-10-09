function [g,h]=noma2_fmin(x,miu,gamma1,gamma2)
%situation:h2<h1,D1^2*beta1>D2^2*beta2
%D1=x(1),D2=x(2),p1=x(3),p2=x(4),beta1=x(5),beta2=x(6),D0=x(7),w1=x(8),w(2)=x(9)
%miu2=0.5;
Pmax=1;
% gamma1=4;
% gamma2=6;
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
g=[(2^gamma1-1)*(x(4)*rou0+delta/0.0018*x(1)^2.2/x(5))-x(3)*rou0
(2^gamma2-1)/rou0*delta/0.0016*x(2)^2.2/x(6)-x(4)
%x(6)-2^gamma1*(2^gamma2-1)/(2^gamma1-1)*x(5);%线性条件
-x(6)/(x(2)^2.2)*0.0016+x(5)/(x(1)^2.2)*0.0018
x(3)+x(4)-Pmax
miu*x(7)-x(1)
(1-miu)*x(7)-x(2)
];
%abs(x(1)^2-x(2)^2)-2;];
h=[x(5)+x(6)-1
  ];