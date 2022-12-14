function[alpha_update,z]=feasibility_initial_point(user,match_c,do,tr,h,h_origin,alpha,p,M,channel,mode,gamma)
if mode==1
M=M/2;
end
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
cvx_begin 
cvx_precision best
variable x(channel,max_c);
expressions X(channel,max_c);
variable R(M,M) hermitian semidefinite;
variable T(M,M) hermitian semidefinite;
variable z 
%rate(j)=log2(1+p1*h(j)/delta*2)
minimize z
subject to
z>=0;
for i=1:user
[m,n]=find(match_c==i);
1+x(m,n)+z>=2^gamma;
end
for i=1:channel
    for j=1:max_c
        H=h(:,i,j)*h(:,i,j)'*1e11;
        if do(i,j)==2   %strong user
            if tr(i,j)==0  % t user 
                p(i,j)*real(trace(T*H))+z >= x(i,j);
            else
                p(i,j)*real(trace(R*H))+z>= x(i,j);
            end
        else   %weak users
            if tr(i,j)==0
                p(i,j)*real(trace(T*H)) +z>= p(i,mod(j,2)+1)*(alpha(i)/2*(real(trace(T*H)))^2+1/2/alpha(i)*x(i,j)^2)+x(i,j);
            else
                p(i,j)*real(trace(R*H)) +z>= p(i,mod(j,2)+1)*(alpha(i)/2*(real(trace(R*H)))^2+1/2/alpha(i)*x(i,j)^2)+x(i,j);
            end                    
        end         
    end
end

if mode==1
diag(R)==1*ones(M,1);
diag(T)==1*ones(M,1);
else
diag(R)+diag(T)==1*ones(M,1);   
end
cvx_end
%% update parameters
for i=1:channel
    for j=1:max_c
        H=h(:,i,j)*h(:,i,j)';
        if do(i,j)==2   %strong user
            if tr(i,j)==0  % t user 
                Hc(i,j)=real(trace(T*H));
                opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/delta);
            else
                Hc(i,j)=real(trace(R*H));
                opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/delta);
            end
        else   %weak users
            if tr(i,j)==0
                Hc(i,j)=real(trace(T*H));
                opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/(p(i,mod(j,2)+1)*Hc(i,j)+delta));  %same as X(i,j)
                alpha_update(i)=x(i,j)/Hc(i,j)/1e11;
            else
                Hc(i,j)=real(trace(R*H));
                opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/(p(i,mod(j,2)+1)*Hc(i,j)+delta));
                alpha_update(i)=x(i,j)/Hc(i,j)/1e11;
            end                    
        end         
    end
end

%update decoding order
for i=1:channel
    for j=1:max_c
       if Hc(i,j)>=Hc(i,mod(j,2)+1)
           do_update(i,j)=2;
       else
           do_update(i,j)=1;
       end
    end
end

for i=1:channel
    for j=1:max_c
        opt_rate(i,j)=log2(1+x(i,j));
    end
end


%TR=[1,1,1,1;1,1,1,1;1,1,1,1;0,0,0,0;0,0,0,0;0,0,0,0]; %need to be revised
%TR=[1,1,1,1;1,1,1,1;1,1,1,1;0,0,0,0;0,0,0,0;0,0,0,0;1,1,1,1;0,0,0,0]; %need to be revised 8-user
TR=[1,1,1,1,1;1,1,1,1,1;1,1,1,1,1;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;1,1,1,1,1;0,0,0,0,0;1,1,1,1,1;0,0,0,0,0]; 
for i=1:user
    for j=1:channel
        temp_h=h_origin(:,i,j)*h_origin(:,i,j)';
        if TR(i,j)==0
            opt_H(i,j)=real(trace(T*temp_h));
        else
            opt_H(i,j)=real(trace(R*temp_h));
        end
    end
end
%save('itr_matrx','Hc');


