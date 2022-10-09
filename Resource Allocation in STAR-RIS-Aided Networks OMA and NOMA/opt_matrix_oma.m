function[Hc,opt_H,R,T,opt_rate]=opt_matrix_oma(match_c,tr,h,h_origin,p,M,channel,gamma)
% clear all;
% clc;
delta=1e-11; %delta^2=-110dB
max_c=2; %maximum users a channel contains
user=6;
%% CVX
cvx_begin 
variable x(channel,max_c);
%expressions X(channel,max_c);
variable R(M,M) hermitian semidefinite;
variable T(M,M) hermitian semidefinite;
%rate(j)=log2(1+p1*h(j)/delta*2)
for i=1:channel
    for j=1:max_c
        X(i,j)=1/2*log2(1+2*x(i,j));
    end
end
maximize sum(sum(X))
subject to
for i=1:user
[m,n]=find(match_c==i);
        %log2(1+x(m,n))>=gamma;
        1+x(m,n)>=2^(2*gamma)/2;
end

for i=1:channel
    for j=1:max_c
        H=h(:,i,j)*h(:,i,j)'*1e11;
            if tr(i,j)==0  % t user 
                p(i,j)*real(trace(T*H)) >= x(i,j);
            else
                p(i,j)*real(trace(R*H)) >= x(i,j);
            end        
    end
end
diag(R)+diag(T)==1*ones(M,1);

cvx_end
%% update parameters
for i=1:channel
    for j=1:max_c
        H=h(:,i,j)*h(:,i,j)';
            if tr(i,j)==0  % t user 
                Hc(i,j)=real(trace(T*H));     
            else   %weak users
                Hc(i,j)=real(trace(R*H));     
            end     
     opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/delta);                
     end         
end



TR=[1,1,1,1;1,1,1,1;1,1,1,1;0,0,0,0;0,0,0,0;0,0,0,0]; %need to be revised
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



