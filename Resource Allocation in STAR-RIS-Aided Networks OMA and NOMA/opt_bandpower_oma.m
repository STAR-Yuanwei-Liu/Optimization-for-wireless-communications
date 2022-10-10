function[opt_rate,p,w]=opt_bandpower_oma(Hc,channel,Pmax,gamma,mode)

delta=1e-11; %delta^2=-110dB
max_c=2; %maximum users a channel contains
Hc=Hc*1e11;
%% CVX
%cvx_begin quiet
cvx_begin
%cvx_precision best
variable p(channel,max_c) 
variable w(channel,max_c)
%rate(j)=log2(1+p1*h(j)/delta*2)

for i=1:channel
    for j=1:max_c
        r(i,j)=-rel_entr(w(i,j),w(i,j)+Hc(i,j)*p(i,j))*log2(exp(1));
    end
end
maximize sum(sum(r))
subject to
for i=1:channel
    for j=1:max_c
        r(i,j)>=gamma;
        p(i,j)>=0;
        w(i,j)>=0;
        if mode==2
        w(i,j)==0.5;
        end
        w(i,1)+w(i,2)<=1;
    end
end

sum(sum(p))<=Pmax;

cvx_end
%% update parameters


for i=1:channel
    for j=1:max_c
        opt_rate(i,j)=r(i,j);
    end
end




