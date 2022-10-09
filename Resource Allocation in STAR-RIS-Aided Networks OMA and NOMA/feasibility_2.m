function[alpha_update,z]=feasibility_2(do,alpha,Hc,gamma,channel,Pmax)
% if mode==1
% M=M/2;
% end
% clear all;
% clc;
delta=1e-11; %delta^2=-110dB
%M=10;
max_c=2; %maximum users a channel contains
%channel=4;
%user=6;
% do=zeros(channel,max_c); %decoding order
% tr=zeros(channel,max_c);  % T or R user
% h=zeros(M,channel,max_c);
%load('initial','do','tr','h_origin','h');
% do=[2,1;2,1;2,1;2,1];
%p=ones(channel,max_c)*0.1;
%alpha=0.0036*ones(channel,1);
%load('itr_power','p','alpha_next');
%alpha=alpha_next;
% alpha=[[0.091085819654193,0.046501530707538,0.086090671650175,0.016994315589394]];
%load('channel','Q','h','rate');
%% CVX
%cvx_begin quiet
cvx_begin 
%cvx_precision best
variable x(channel,max_c);
expressions X(channel,max_c);
variable p(channel,max_c) 
variable z
%rate(j)=log2(1+p1*h(j)/delta*2)
for i=1:channel
    for j=1:max_c
        X(i,j)=log2(1+x(i,j));
    end
end
minimize z
subject to
z>=0;
for i=1:channel
    for j=1:max_c
        1+x(i,j)+z>=2^gamma;
        p(i,j)>=0;
    end
end

sum(sum(p))<=Pmax+z;
Hc=Hc*1e11;
for i=1:channel
    for j=1:max_c
        if do(i,j)==2   %strong user
%             if tr(i,j)==0  % t user 
                p(i,j)*Hc(i,j)+z>= x(i,j);
        else   %weak users
                p(i,j)* Hc(i,j)+z>= (Hc(i,j))*(alpha(i)/2*p(i,mod(j,2)+1)^2+1/2/alpha(i)*x(i,j)^2)+x(i,j);               
        end         
    end
end

cvx_end
%% update parameters

%update decoding order
for i=1:channel
if do(i,1)==2
    alpha_update(i)=x(i,2)/p(i,2);
else
    alpha_update(i)=x(i,1)/p(i,1);
end
end

for i=1:channel
    for j=1:max_c
        opt_rate(i,j)=log2(1+x(i,j));
    end
end



