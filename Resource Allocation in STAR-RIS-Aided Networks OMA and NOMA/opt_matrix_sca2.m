function[opt_rate,alpha,a0,b0,Hc]=opt_matrix_sca2(M,h,channel,tr,do,a,b,x,p)
a0=a*10^5.5;
b0=b*10^5.5;
max_c=2;
%h=h*1e11;
delta=1e-11;
cvx_begin
cvx_precision best
variable et(M,1) complex  %phase shift of T
variable er(M,1) complex  %phase shift of R
%variable t nonnegative
expression a(channel,max_c)
expression b(channel,max_c)
variable t(channel,1)
% minimize sum(sum(t))
maximize sum(t)


for i=1:channel
    t(i)>=0;   
    for j=1:max_c

        if tr(i,j)==0
        a(i,j)=real(et'*h(:,i,j))*10^5.5;
        b(i,j)=imag(et'*h(:,i,j))*10^5.5;
        else
        a(i,j)=real(er'*h(:,i,j))*10^5.5;
        b(i,j)=imag(er'*h(:,i,j))*10^5.5;
        end
    end
end
subject to
for i=1:channel
    if do(i,1)==2
%         1e11*(a0(i,1)^2+b0(i,1)^2+2*a0(i,1)*(a(i,1)-a0(i,1))+2*b0(i,1)*(b(i,1)-b0(i,1))-t(i,1))>=1e11*(a(i,2)^2+b(i,2)^2);   %do constraint
%         1e11*(a0(i,1)^2+b0(i,1)^2+2*a0(i,1)*(a(i,1)-a0(i,1))+2*b0(i,1)*(b(i,1)-b0(i,1))-t(i,1))>=1e11*(x(i,1)/p(i,1)*delta);
%         1e11*(a0(i,2)^2+b0(i,2)^2+2*a0(i,2)*(a(i,2)-a0(i,2))+2*b0(i,2)*(b(i,2)-b0(i,2))-t(i,2))>=1e11*(x(i,2)*p(i,1)/p(i,2)*(a(i,2)^2)+x(i,2)*p(i,1)/p(i,2)*(b(i,2)^2)+x(i,2)*delta);       
        a0(i,1)^2+b0(i,1)^2+2*a0(i,1)*(a(i,1)-a0(i,1))+2*b0(i,1)*(b(i,1)-b0(i,1))>=a(i,2)^2+b(i,2)^2;   %do constraint
        a0(i,1)^2+b0(i,1)^2+2*a0(i,1)*(a(i,1)-a0(i,1))+2*b0(i,1)*(b(i,1)-b0(i,1))-t(i)>=x(i,1)/p(i,1);
        a0(i,2)^2+b0(i,2)^2+2*a0(i,2)*(a(i,2)-a0(i,2))+2*b0(i,2)*(b(i,2)-b0(i,2))>=x(i,2)*p(i,1)/p(i,2)*(a(i,2)^2)+x(i,2)*p(i,1)/p(i,2)*(b(i,2)^2)+x(i,2)/p(i,2);       
    else
%         1e11*(a0(i,2)^2+b0(i,2)^2+2*a0(i,2)*(a(i,2)-a0(i,2))+2*b0(i,2)*(b(i,2)-b0(i,2))-t(i,2))>=1e11*(a(i,1)^2+b(i,1)^2); %do constraint
%         1e11*(a0(i,2)^2+b0(i,2)^2+2*a0(i,2)*(a(i,2)-a0(i,2))+2*b0(i,2)*(b(i,2)-b0(i,2))-t(i,2))>=1e11*(x(i,2)/p(i,2)*delta);  
%         1e11*(a0(i,1)^2+b0(i,1)^2+2*a0(i,1)*(a(i,1)-a0(i,1))+2*b0(i,1)*(b(i,1)-b0(i,1))-t(i,1))>=1e11*(x(i,1)*p(i,2)/p(i,1)*(a(i,1)^2+b(i,1)^2)+x(i,1)*delta);
        a0(i,2)^2+b0(i,2)^2+2*a0(i,2)*(a(i,2)-a0(i,2))+2*b0(i,2)*(b(i,2)-b0(i,2))>=a(i,1)^2+b(i,1)^2; %do constraint
        a0(i,2)^2+b0(i,2)^2+2*a0(i,2)*(a(i,2)-a0(i,2))+2*b0(i,2)*(b(i,2)-b0(i,2))-t(i)>=x(i,2)/p(i,2);  
        a0(i,1)^2+b0(i,1)^2+2*a0(i,1)*(a(i,1)-a0(i,1))+2*b0(i,1)*(b(i,1)-b0(i,1))>=x(i,1)*p(i,2)/p(i,1)*(a(i,1)^2+b(i,1)^2)+x(i,1)/p(i,1);
    end
end
for i=1:M
norm(et(i))^2+norm(er(i))^2<=1;
end
cvx_end


for i=1:channel
    for j=1:max_c
       Hc(i,j)=(a(i,j)^2+b(i,j)^2)/1e11;
    end
end

for i=1:channel
    for j=1:max_c
        if do(i,j)==2   %strong user
            x(i,j)=p(i,j)*Hc(i,j)/delta;
            opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/delta);
        else   %weak users
            x(i,j)=p(i,j)*Hc(i,j)/(p(i,mod(j,2)+1)*Hc(i,j)+delta);            
            opt_rate(i,j)=log2(1+p(i,j)*Hc(i,j)/(p(i,mod(j,2)+1)*Hc(i,j)+delta));                               
        end         
    end
end

for i=1:channel
if do(i,1)==2
alpha(i)=x(i,2)/p(i,1);    
else
alpha(i)=x(i,1)/p(i,2);      
end
end

for i=1:channel
    for j=1:max_c
        if tr(i,j)==0
        a(i,j)=real(et'*h(:,i,j));
        b(i,j)=imag(et'*h(:,i,j));
        else
        a(i,j)=real(er'*h(:,i,j));
        b(i,j)=imag(er'*h(:,i,j));
        end
    end
end

a0=a;
b0=b;