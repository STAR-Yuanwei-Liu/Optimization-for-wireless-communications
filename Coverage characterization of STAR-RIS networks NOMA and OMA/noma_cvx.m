function[d1,d2]=noma_cvx(miu,gamma1,gamma2,K)

Pmax=1;
% gamma1=5;
% gamma2=5;
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
% if K==10
%     load('channel-k=10.mat')
% else
%     load('channel-k=2.mat')
% end
cvx_begin 
variable D1 nonnegative 
variable D2 nonnegative  
variable p1 nonnegative 
variable p2 nonnegative 
variable d2 nonnegative;
variable d1 nonnegative;
variable beta1 nonnegative 
variable beta2 nonnegative 
variable D0 nonnegative 
%variable tau(N) nonnegative
%expression v(100,1)

%maximize miu*D1+(1-miu)*D2
%maximize (miu*sqrt(D1)+(1-miu)*sqrt(D2))
%maximize pow_pos(D1,2)+pow_pos(D2,2)
maximize D0
subject to 
D1>=miu*D0;
D2>=(1-miu)*D0;
d2>=pow_p(D2,1.1);
d1>=pow_pos(D1,1.1);
(2^gamma2-1)*(p1*rou0+delta/temp2*quad_over_lin(d2,beta2))<=p2*rou0;
(2^gamma1-1)/rou0*delta/temp1*quad_over_lin(d1,beta1)<=p1;
%quad_over_lin(D1,beta1)-(D20^2/beta20+2*D20^1/beta20*(D2-D20)-D20^2/beta20^2*(beta2-beta20))<=0;
%D2 >= D1;
D1 >= 1;
D2 >= 1;
beta1 + beta2 <= 1;
p1 + p2 <= Pmax;
cvx_end
