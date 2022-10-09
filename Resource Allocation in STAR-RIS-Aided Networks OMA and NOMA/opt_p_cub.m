function[alpha_update,opt_rate,p,x]=opt_p_cub(do,alpha,Hc,gamma,channel,Pmax)
% if mode==1
% M=M/2;
% end
% clear all;
% clc;
delta=1e-11; %delta^2=-110dB
%M=10;
max_c=2; %maximum users a channel contains

%% CVX
%cvx_begin quiet
cvx_begin
%cvx_precision best
variable x(channel,max_c);
expressions X(channel,max_c);
variable p(channel,max_c) 
%rate(j)=log2(1+p1*h(j)/delta*2)
for i=1:channel
    for j=1:max_c
        X(i,j)=log2(1+x(i,j));
    end
end
maximize sum(sum(X))
subject to
for i=1:channel
    for j=1:max_c
        1+x(i,j)>=2^gamma;
        p(i,j)>=0;
    end
end

sum(sum(p))<=Pmax;
Hc=Hc*1e11;
for i=1:channel
    for j=1:max_c
        if do(i,j)==2   %strong user
%             if tr(i,j)==0  % t user 
                p(i,j)*Hc(i,j) >= x(i,j);
        else   %weak users
                p(i,j)* Hc(i,j)>= (Hc(i,j))*(alpha(i)/2*p(i,mod(j,2)+1)^2+1/2/alpha(i)*x(i,j)^2)+x(i,j);               
        end         
    end
end

cvx_end
%% update parameters

%update decoding order
for i=1:channel
if do(i,1)==2
    alpha_update(i)=x(i,2)/p(i,1);  %%!!!! 不是除以自己的p
else
    alpha_update(i)=x(i,1)/p(i,2);
end
end

for i=1:channel
    for j=1:max_c
        opt_rate(i,j)=log2(1+x(i,j));
    end
end



