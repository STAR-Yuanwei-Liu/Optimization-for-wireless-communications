%Use to genarate channel
clear all;
clc;
la=[0,0]; %location of AP
lr=[50,0]; %location of ris
u1=[48,0];
u2=[52,0];
D1=norm(lr-u1);
D2=norm(lr-u2);
dar=norm(la-lr);
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
M=20; %RIS element
K=10; %Rician fading factor

for i=1:M
    glos(i)=exp(j*pi*(i-1));
    r1los(i)=exp(j*pi*(i-1));
    r2los(i)=exp(j*pi*(i-1));
end
an=angle(r1los);
% glos=ones(1,M);
% r1los=ones(1,M);
% r2los=ones(1,M);
var=1;
gnlos=sqrt(var/2)*(randn(1,M)+1i*randn(1,M)); %CSCG
r1nlos=sqrt(var/2)*(randn(1,M)+1i*randn(1,M));
r2nlos=sqrt(var/2)*(randn(1,M)+1i*randn(1,M));
g=sqrt(rou0/(dar^2.2))*(sqrt(K/(K+1))*glos+sqrt(1/(K+1))*gnlos);
r1=sqrt(rou0/(D1^2.2))*(sqrt(K/(K+1))*r1los+sqrt(1/(K+1))*r1nlos);
r2=sqrt(rou0/(D2^2.2))*(sqrt(K/(K+1))*r2los+sqrt(1/(K+1))*r2nlos);
r1avr=(sqrt(K/(K+1))*r1los+sqrt(1/(K+1))*r1nlos);
r2avr=(sqrt(K/(K+1))*r2los+sqrt(1/(K+1))*r2nlos);
q1=r1avr*diag(g); %r1 bar
p=angle(q1); 
q2=r2avr*diag(g);
temp1=norm(sum(q1))^2;
optq1=zeros(1,M);
optq2=zeros(1,M);
for j=1:M
    optq1(j)=norm(q1(j));
    optq2(j)=norm(q2(j));
end
temp1=sum(optq1)^2;
temp2=sum(optq2)^2;
h1=norm(sum(q1))^2*rou0/(D1^2.2); %|h1|^2
h2=norm(sum(q2))^2*rou0/(D2^2.2);
