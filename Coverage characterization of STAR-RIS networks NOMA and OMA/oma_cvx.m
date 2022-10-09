%function[d1,d2]=oma1_cvx(miu,gamma1,gamma2,K)
miu=0.4;
%miu=0.5;
Pmax=1;
gamma1=5;
gamma2=5;
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
K=10;
% if K==10
%     load('channel-k=10.mat')
% else
%     load('channel-k=2.mat')
% end
cvx_begin 
variable D1(1) nonnegative ;
variable D2(1) nonnegative  ;
variable D0 nonnegative 
variable p1 nonnegative ;
variable p2 nonnegative ;
variable d2 nonnegative;
variable d1 nonnegative;
variable tau1 nonnegative ;
variable tau2 nonnegative ;
% variable beta1 nonnegative ;
% variable beta2 nonnegative;
variable beta1(100,100) semidefinite ;
variable beta2(100,100) semidefinite;
%maximize D1
%maximize miu*D1+(1-miu)*D2
%maximize (miu*sqrt(D1)+(1-miu)*sqrt(D2))
maximize D0
%maximize pow_pos(D1,2)+pow_pos(D2,2)
%maximize (D1^2);
subject to 
%D = miu*sqrt(D1)+(1-miu)*sqrt(D2)
D1>=miu*D0;
D2>=(1-miu)*D0;
%0.5*(2^(2*gamma2)-1)/rou0*delta/0.0016*quad_over_lin(D2,beta2)<=p2;
d2>=pow_p(D2,1.1);
d1>=pow_pos(D1,1.1);
%log(D1)==4;
%0.5*(2^(2*gamma2)-1)/rou0*delta/0.0016*quad_over_lin(D2,p2)<=tau2;
%0.5*(2^(2*gamma1)-1)/rou0*delta/0.0017*quad_over_lin(D1,p1)<=tau1;
0.5*(2^(2*gamma2)-1)/rou0*delta*quad_over_lin(d2,p2)<=tau2;
0.5*(2^(2*gamma1)-1)/rou0*delta*quad_over_lin(d1,p1)<=tau1;

trace(beta1*(optq1'*optq1))>=tau1;
trace(beta2*(optq2'*optq2))>=tau2;
for i=1:100
beta1(i,i) + beta2(i,i) <= 1;
end
% beta1==semidefinite(10);
% beta2==semidefinite(10);
D1 >= 1;
D2 >= 1;
%p1>=0.4;
%p2>=0.4;
% beta1>=0.499;
% beta2>=0.499;

p1 + p2 <= Pmax;
%p1 <= Pmax;
%p2 <= Pmax;
cvx_end
