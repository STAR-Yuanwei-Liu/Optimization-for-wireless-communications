% Please refer to C. Wu, Y. Liu, X. Mu, X. Gu, and O. A. Dobre, °∞Coverage characterization of STAR-RIS networks: NOMA and OMA,°± IEEE Commun. Lett.,
%vol. 25, no. 9, pp. 3036®C3040, Sept. 2021.
% Date:2021/11/24
% Use to plot Fig.3
clear all;
clc;
%% initialization
Pmax=1;
rou0=0.001; %-30dB  2.4Ghz-40dB
delta=1e-11; %delta^2=-110dB
A = [];
b = [];
Aeq = [];
beq = [];
Aeq2=diag([0 0 0 0 1 1 0 0 0]);
beq2=[0 0 0 0 0.5 0.5 0 0 0];
Aeq3=diag([0 0 0 0 1 1 0 ]);
beq3=[0 0 0 0 0.5 0.5 0];
lb = [1 1 0 0 0 0 0];
lb2=[1 1 0 0 0 0 0 0 0] 
ub = [];
%x0=[2 2 0 0 0 0 0];
x0=[2 2 0.5 0.5 0.5 0.5 0.5];%≥ı÷µ
%x1=[2 2 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
x1=[2 2 0 0 0 0 0 0 0];%≥ı÷µ
miu=0.6;
rate=2:0.5:6;
N=length(rate);

%% main body

% tmp1=zeros(1,N);
% tmp2=zeros(1,N);
% 
% for i=1:N
% [x fval] = fmincon(@(x) obj_fmin(x,0.5),x0,A,b,Aeq,beq,lb,ub,@(x) oma1_fmin(x,miu,rate(i),5))
% tmp1(i)=x(2);
% tmp2(i)=x(7);
% end

tmp3=zeros(1,N);
tmp4=zeros(1,N);
for i=1:N
[x fval] = fmincon(@(x) obj_fmin(x,0.5),x1,A,b,Aeq,beq,lb2,ub,@(x) oma2_fmin(x,miu,rate(i),5))
tmp3(i)=x(2);
tmp4(i)=x(7);
end

% decoding order1
tmp5=zeros(1,N);
tmp6=zeros(1,N);
for i=1:N
[x fval] = fmincon(@(x) obj_fmin(x,0.5),x0,A,b,Aeq,beq,lb,ub,@(x) noma2_fmin(x,miu,rate(i),5))%% can use noma_cvx.m instead
tmp5(i)=x(1);
tmp6(i)=x(7);
end
%plot(r,tmp6,'s-');

% decoding order2
tmp7=zeros(1,N);
tmp8=zeros(1,N);
for i=1:N
[x fval] = fmincon(@(x) obj_fmin(x,0.5),x0,A,b,Aeq,beq,lb,ub,@(x) noma1_fmin(x,miu,rate(i),5))
tmp7(i)=x(1);
tmp8(i)=x(7);
end
%plot(r,tmp8,'*-');


% tmp9=zeros(1,N);
% tmp10=zeros(1,N);
% for i=1:N
% [x fval] = fmincon(@(x) obj_fmin(x,0.5),x0,A,b,Aeq,beq,lb,ub,@(x) oma1_baseline(x,miu,rate(i),5))
% tmp9(i)=x(2);
% tmp10(i)=x(7);
% end

tmp11=zeros(1,N);
tmp12=zeros(1,N);
for i=1:N
[x fval] = fmincon(@(x) obj_fmin(x,0.5),x1,A,b,Aeq,beq,lb2,ub,@(x) oma2_baseline(x,miu,rate(i),5))
tmp11(i)=x(2);
tmp12(i)=x(7);
end

% 
tmp13=zeros(1,N);
tmp14=zeros(1,N);
for i=1:N
[x fval] = fmincon(@(x) obj_fmin(x,0.5),x0,A,b,Aeq,beq,lb,ub,@(x) noma2_baseline(x,miu,rate(i),5))
tmp13(i)=x(1);
tmp14(i)=x(7);
end
%plot(r,tmp6,'s-');


tmp15=zeros(1,N);
tmp16=zeros(1,N);
for i=1:N
[x fval] = fmincon(@(x) obj_fmin(x,0.5),x0,A,b,Aeq,beq,lb,ub,@(x) noma1_baseline(x,miu,rate(i),5))
tmp15(i)=x(1);
tmp16(i)=x(7);
end
%plot(r,tmp8,'*-');

% combine docoding order
noma=zeros(1,N);
for i=1:N
noma(i)=max(tmp6(i),tmp8(i));
end

noma_bs=zeros(1,N);
for i=1:N
noma_bs(i)=max(tmp14(i),tmp16(i));
end
%plot(rate,miu*tmp2,'o-') %OMA-TYPE2
%plot(r,tmp2);

%% plot
plot(rate,noma,'*-'); %NOMA
hold on;
plot(rate,tmp4,'o-'); %OMA
hold on;

plot(rate,noma_bs,'*-.');
%plot(rate,miu*tmp10,'o-') %OMA-TYPE2-CR
%plot(r,tmp2);
hold on;
plot(rate,tmp12,'o-.');
hold on;


%legend('OMA-1','OMA-2','NOMA-part1','NOMA-part2');
%legend('OMA-I-STAR','OMA-II-STAR','NOMA-STAR','OMA-I-CR','OMA-II-CR','NOMA-CR');
legend('NOMA-STAR','OMA-STAR','NOMA-CR','OMA-CR');
xlabel('T user QoS requirements \gamma_t(bps/Hz)');
ylabel('Total coverage range D_0(m)');
box on;
