function[do_new,Hc,R,T]=do_max_combined(tr,h,channel,max_c,M)

cvx_begin 
variable R(M,M) hermitian semidefinite;
variable T(M,M) hermitian semidefinite;
cg=0; %combined channel gain
for i=1:channel
    for j=1:max_c
        H=h(:,i,j)*h(:,i,j)'*1e11;
            if tr(i,j)==0  % t user 
                cg=cg+real(trace(T*H)); 
            else
                cg=cg+real(trace(R*H));
            end                                     
    end
end
%subject to 

maximize cg
subject to
diag(R)+diag(T)==ones(M,1);
cvx_end



for i=1:channel
    for j=1:max_c
       H=h(:,i,j)*h(:,i,j)';
            if tr(i,j)==0  % t user 
                Hc(i,j)=real(trace(T*H)); 
            else
                Hc(i,j)=real(trace(R*H));
            end                                  
    end
end


for i=1:channel
if Hc(i,1)>=Hc(i,2)
    do_new(i,1)=2;
    do_new(i,2)=1;
else
    do_new(i,1)=1;
    do_new(i,2)=2;
end
end


[v1,d1]=eig(R);
[v2,d2]=eig(T);
er=v1(:,M)*sqrt(d1(M,M));
et=v2(:,M)*sqrt(d2(M,M));
for i=1:channel
    for j=1:max_c
        if tr(i,j)==0
        a0(i,j)=real(et'*h(:,i,j));
        b0(i,j)=imag(et'*h(:,i,j));
        else
        a0(i,j)=real(er'*h(:,i,j));
        b0(i,j)=imag(er'*h(:,i,j));
        end
    end
end





