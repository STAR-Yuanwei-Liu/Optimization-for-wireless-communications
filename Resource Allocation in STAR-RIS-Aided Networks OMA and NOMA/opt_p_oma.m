function[p,opt_rate,w]=opt_p_oma(user,match_c,Hc,Pmax,channel,mode,gamma)
delta=1e-11; %delta^2=-110dB

max_c=2; %maximum users a channel contains
Hc=Hc*1e11;
cvx_begin 
    variable p(channel,max_c) 
    variable w(channel,max_c) 
    expressions X(channel,max_c) 
    for i=1:channel
        for j=1:max_c
            X(i,j)=-rel_entr(w(i,j),w(i,j)+Hc(i,j)*p(i,j))*log2(exp(1));
        end
    end
maximize sum(sum(X))
    %find r1
    subject to
for i=1:channel
    w(i,1)+w(i,2)<=1;    
    for j=1:max_c
        w(i,j)>=0;
        if mode==2
        w(i,j)==0.5;
        end
        -rel_entr(w(i,j),w(i,j)+Hc(i,j)*p(i,j))*log2(exp(1))>=gamma;
    end
end
sum(sum(p))<=Pmax;
for i=1:channel
    for j=1:max_c
       p(i,j)>=0;
    end
end
cvx_end

Hc=Hc/1e11;
 for i=1:channel
        for j=1:max_c
            opt_rate(i,j)=w(i,j)*log2(1+p(i,j)*Hc(i,j)/delta/w(i,j));
        end
 end
