%function[h_origin,opt_h]=channel_generate(M,channel,mode,user)
% R = raylrnd(1,1,10000);
% a=sum(R)/10000
% clear all;
% 
% mode=0;
% M=30;
% mode=0;
% user=6;
% channel=3;

if mode==1
M=M/2;
end
la=[0,0]; %location of AP
lr=[50,0]; %location of star-ris
%user=6;  %%revise
% channel=2; %%revise
% M=20;

% xu=[40,40,40,65,65,65,35,70,45,55]; %x-coordinates of users
% yu=[15,0,-20,25,0,-35,21,-22,2,-3]; %y-coordinates of users
% xu=[40,40,40,65,65,65,45,55]; %x-coordinates of users
% yu=[15,0,-20,25,0,-35,6,-8]; %y-coordinates of users
% xu=[40,40,40,40,55,55,55,55]; %x-coordinates of users
% yu=[15,5,-5,-15,4,2,-2,-4]; %y-coordinates of users

% xu=[40,40,40,60,60,60]; %x-coordinates of users
% yu=[5,0,-5,5,0,-5]; %y-coordinates of users
% xu=[45,39,41,65,53,56]; %x-coordinates of users
% yu=[6,2,-5,8,-6,-3]; %y-coordinates of users

% xu=[45,42,58,55];
% yu=[5,-8,8,-5];
% xu=[45,45,55,55];
% yu=[5,-5,5,-5];
if user==6
xu=[47,45,47,53,55,53];
yu=[4,0,-4,4,0,-4];
elseif user==4
% xu=[45,45,55,55];
% yu=[5,-5,5,-5];
xu=[47,47,53,53];
yu=[4,-4,4,-4];
% xu=[45,42,58,55];
% yu=[5,-8,8,-5];
elseif user==8
a=cos(60/180*pi);
b=sin(60/180*pi);
xu=[47,45,47,45-a,53,55,53,45+a];
yu=[4,0,-4,b,4,0,-4,-b];
end

for i=1:user
    disau(i)=sqrt((xu(i)-la(1))^2+(yu(i)-la(2))^2);
    disru(i)=sqrt((xu(i)-lr(1))^2+(yu(i)-lr(2))^2);
end
disar=norm(la-lr);
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
%M=40; %RIS element
K=2; %Rician fading 
var=1;
nlos=zeros(user,M);
los=ones(1,M);
h_origin=zeros(M,user,channel);


 %AP to RIS
for i=1:channel
    gnlos=sqrt(var/2)*(randn(1,M)+1i*randn(1,M));
    g=sqrt(rou0/(disar^2.2))*(sqrt(K/(K+1))*los+sqrt(1/(K+1))*gnlos); 
    for j=1:user
    nlos(j,:)=sqrt(var/2)*(randn(1,M)+1i*randn(1,M));   
    f(j,:)=sqrt(rou0/(disru(j)^2.8))*(sqrt(K/(K+1))*los+sqrt(1/(K+1))*nlos(j,:));
    h_origin(:,j,i)=diag(g)*(f(j,:)'); 
    end
end

for i=1:channel
    for j=1:user
        for k=1:M
            optq(k,j,i)=norm(h_origin(k,j,i));
        end
    end
end

for i=1:channel
    for j=1:user
        opt_h(j,i)=sum(optq(:,j,i))^2;  %先取norm，再求和，再平方
    end
end

%save('channel_60_6_asy','h_origin','opt_h');
