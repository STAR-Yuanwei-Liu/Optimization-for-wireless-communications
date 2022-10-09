function[h_origin,opt_h]=channel_generate2(dis_t,M,mode)
if mode==1
M=M/2;
end
user=6;
channel=3;
dis_r=10-dis_t;
la=[0,0]; %location of AP
lr=[50,0]; %location of star-ris

disru=[dis_r,dis_r,dis_r,dis_t,dis_t,dis_t];

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

