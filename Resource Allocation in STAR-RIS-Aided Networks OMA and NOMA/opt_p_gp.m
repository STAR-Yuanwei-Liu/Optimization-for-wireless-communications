function[p,alpha_next,opt_rate]=opt_p_gp(user,match_c,do,Hc,Pmax,channel,gamma)
delta=1e-11; %delta^2=-110dB
%M=10;
max_c=2; %maximum users a channel contains
%channel=4;
%user=6;
%Pmax=0.8;
% load('initial','do','tr','h_origin','opt_H','h','Hc');
% load('itr_matrx','Hc');
Mc=1e-11./Hc; %delta/Hc
%do=[2,1;2,1;2,1;2,1];
%w=[1,1.1;1.2,1.3;1.4,1.5;1.6,1.7];
cvx_begin gp 
%cvx_begin gp quiet
    variable r(channel,max_c) 
    variable r2
    z=1;
    x=0;
    y=0;
    for i=1:channel
        for j=1:max_c
            %z=z*r(i,j)^w(i,j);
            z=z*r(i,j);
        end
    end
    for i=1:channel
        if do(i,1)==2
            x=x+Mc(i,1)*(r(i,1)*r(i,2))+(Mc(i,2)-Mc(i,1))*r(i,2);
            y=y+Mc(i,2);
        else
            x=x+Mc(i,2)*(r(i,1)*r(i,2))+(Mc(i,1)-Mc(i,2))*r(i,1);
            y=y+Mc(i,1);
        end    
    end
    maximize log(z)/log(2)
    %find r1
    subject to
    x<=y+Pmax;
%     for i=1:channel
%         for j=1:max_c
%            r(i,j)>=2^gamma;
%         end
%     end
    for i=1:channel
        for j=1:max_c
           r(i,j)>=1;
        end
    end
for i=1:user
[m,n]=find(match_c==i);
r(m,n)>=2^gamma;
end
    
    
% r(1,1)>=2^gamma;
% r(1,2)>=2^gamma;
% r(2,1)*r(3,1)>=2^gamma;
% r(2,2)*r(4,2)>=2^gamma;
% r(3,2)>=2^gamma;
% r(4,1)>=2^gamma;

cvx_end

%% update parameters
%  for i=1:channel
%         for j=1:max_c
%            if do(i,j)==2
%            p(i,j)=(2^log2(r(i,j))-1)*Hc(i,j);
%            else
%            p_=(2^log2(r(i,mod(2,j)+1))-1)*Hc(i,mod(2,j)+1);
%            p(i,j)=(p_+Hc(i,j))*(2^log2(r(i,j))-1);
%            alpha_next(i)=(r(i,j)-1)/Hc(i,j);
%            end
%         end
%  end

 for i=1:channel
        for j=1:max_c
           if do(i,j)==2
           p(i,j)=(r(i,j)-1)*Mc(i,j);
           else
           p_=((r(i,mod(j,2)+1))-1)*Mc(i,mod(j,2)+1);
           p(i,j)=(p_+Mc(i,j))*(r(i,j)-1);
           alpha_next(i)=(r(i,j)-1)/Hc(i,j)/1e11;
           end
        end
 end
 for i=1:channel
        for j=1:max_c
            opt_rate(i,j)=log2(r(i,j));
        end
 end
%save('itr_power','p','alpha_next');
