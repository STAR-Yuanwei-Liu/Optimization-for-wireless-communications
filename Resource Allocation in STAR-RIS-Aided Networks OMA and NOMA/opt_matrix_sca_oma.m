function[opt_rate,a0,b0,Hc]=opt_matrix_sca_oma(M,h,channel,tr,a,b,p,w,gamma)
% a0=a*10^5.5;
% b0=b*10^5.5;
a0=a;
b0=b;
max_c=2;
%h=h*1e11;
delta=1e-10;
cvx_begin

variable et(M,1) complex  %phase shift of T
variable er(M,1) complex  %phase shift of R

variable aa(channel,max_c)
variable bb(channel,max_c)
variable y(channel,max_c)
% minimize sum(sum(t))

for i=1:channel 
    for j=1:max_c
        r(i,j)=w(i,j)*log(1+p(i,j)*y(i,j)/w(i,j)/delta)/log(2);
        %r(i,j)=w(i,j)*log(1+p(i,j)*y(i,j)/w(i,j))/log(2);
    end
end
maximize sum(sum(r))

subject to
for i=1:channel 
    for j=1:max_c
%         if tr(i,j)==0
%         a(i,j)==real(et'*h(:,i,j))*10^5.5;
%         b(i,j)==imag(et'*h(:,i,j))*10^5.5;
%         else
%         a(i,j)==real(er'*h(:,i,j))*10^5.5;
%         b(i,j)==imag(er'*h(:,i,j))*10^5.5;
%         end
        if tr(i,j)==0
        aa(i,j)==real(et'*h(:,i,j));
        bb(i,j)==imag(et'*h(:,i,j));
        else
        aa(i,j)==real(er'*h(:,i,j));
        bb(i,j)==imag(er'*h(:,i,j));
        end
    end
end


for i=1:channel
for j=1:max_c
        r(i,j)>=gamma;
        a0(i,j)^2+b0(i,j)^2+2*a0(i,j)*(aa(i,j)-a0(i,j))+2*b0(i,j)*(bb(i,j)-b0(i,j))>=y(i,j);
        %a0(i,j)^2+b0(i,j)^2+2*a0(i,j)*(aa(i,j)-a0(i,j))>=y(i,j);
end
end
for m=1:M
%norm(et(i))^2+norm(er(i))^2<=1;
%pow_p(norm(et(i)),2)+pow_p(norm(er(i)),2)<=1;
square_abs(et(m))+square_abs(er(m))<=1;
%norm(et(i))+norm(er(i))<=1;
end
cvx_end


for i=1:channel
    for j=1:max_c
       Hc(i,j)=(aa(i,j)^2+bb(i,j)^2);
       %Hc(i,j)=(a(i,j)^2+b(i,j)^2)/1e11;
    end
end

for i=1:channel
    for j=1:max_c      
            opt_rate(i,j)=r(i,j);                                  
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



for i=1:channel 
    for j=1:max_c
        r(i,j)=w(i,j)*log(1+p(i,j)*Hc(i,j)/w(i,j)/1e-11)/log(2);
    end
end
opt_rate=r;