clear all;
clc;
beta_r=[0.2,0.3,0.4,0.5,0.6,0.7,0.8];
avr_mse_t=[0.6065,0.6914,0.808,0.97,1.21,1.614,2.424]*1e-14;
avr_mse_r=[2.4218,1.617,1.212,0.97,0.808,0.692,0.605]*1e-14;
avr_mse=[3.0584,2.3394,2.0505,1.97,2.049,2.337,3.06]*1e-14;
avr_mse_ts=[2,2,2,2,2,2,2]*1e-14;
plot(beta_r,avr_mse);
hold on;
plot(beta_r,avr_mse_t);
hold on;
plot(beta_r,avr_mse_r);
hold on;
plot(beta_r,avr_mse_ts);
xlim([0.2 0.8])
%ylim([0 4]*1e-3)
xlabel('\beta_r');
ylabel('MSE');
legend('Energy splitting sum MSE','T user MSE','R user MSE','Time switching sum MSE');


clear all;
clc;
pho=0.001;
d_as=50; %%AP-STAR
d_sut=5;   %STAR-T user
d_aut=55;  %AP-T user
d_sur=5;   %STAR-R user
d_aut=45;  %AP-R user
alpha_as=2.2;
alpha_su=2.8;
alpha_au=3.5;
L_as=pho*d_as^-alpha_as;
L_sut=pho*d_sut^-alpha_su;
L_aut=pho*d_aut^-alpha_au;
L_sur=pho*d_sut^-alpha_su;
L_aur=pho*d_aut^-alpha_au;
beta_t=0.5;
beta_r=0.5;
M=32;
T=2*(M+1);
SNR=140; %%dB P=1W -170dBm/Hz
iter=100;

for t=1:7
beta_t=0.1+0.1*t;
beta_r=1-beta_t;
sum_mse_es=0;
sum_mse_t=0;
sum_mse_r=0;
sum_mse_ts=0;
sum_nmse=0;
sum_nmse_t=0;
sum_nmse_r=0;
for k=1:iter
ht=sqrt(L_aut)*1/sqrt(2)*(randn(1)+1j*randn(1));
hr=sqrt(L_aur)*1/sqrt(2)*(randn(1)+1j*randn(1));
g=sqrt(L_as)*1/sqrt(2)*(randn(1,M)+1j*randn(1,M));
rt=sqrt(L_sut)*1/sqrt(2)*(randn(M,1)+1j*randn(M,1));
rr=sqrt(L_sur)*1/sqrt(2)*(randn(M,1)+1j*randn(M,1));
qt=diag(g)*rt;
qr=diag(g)*rr;

%% ideal 
x=[ht;hr;qt;qr];
sigma_2=10^(-SNR/10);
noise=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(T,1)+1j*randn(T,1));
H_temp=dftmtx(T);
H=H_temp(:,1:2*M+2);
H(:,3:2+M)=H(:,3:2+M)*sqrt(beta_t);
H(:,3+M:2*M+2)=H(:,3+M:2*M+2)*sqrt(beta_r);
y=H*x+noise;
error=pinv(H)*y-x; %%H'/T!=inv(H)
error_t=error(3:2+M);
error_r=error(3+M:2*M+2);
mse_matrix=error*error';
sum_mse_es=sum_mse_es+real(trace(mse_matrix));
sum_mse_t=sum_mse_t+error_t'*error_t;
sum_mse_r=sum_mse_r+error_r'*error_r;





%% Time swithing
H_ts_temp=dftmtx(T/2);
H_ts=H_ts_temp(:,1:M+1);
noise_t_ts=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(T/2,1)+1j*randn(T/2,1));
y_t=H_ts*[ht;qt]+noise_t_ts;
y_t_ts=pinv(H_ts)*y_t;
noise_r_ts=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(T/2,1)+1j*randn(T/2,1));
y_r=H_ts*[hr;qr]+noise_r_ts;
y_r_ts=pinv(H_ts)*y_r;
error_t=y_t_ts-[ht;qt];
error_r=y_r_ts-[hr;qr];
sum_mse_ts=sum_mse_ts+real(trace(error_t*error_t'))+real(trace(error_r*error_r'));

%% NMSE
sum_nmse=sum_nmse+[ht;hr;qt;qr]'*[ht;hr;qt;qr];
sum_nmse_t=sum_nmse_t+[ht;qt]'*[ht;qt];
sum_nmse_r=sum_nmse_r+[hr;qr]'*[hr;qr];
end
avr_mse_es(t)=sum_mse_es/sum_nmse;
avr_mse_t(t)=sum_mse_t/sum_nmse_t;
avr_mse_r(t)=sum_mse_r/sum_nmse_r;
avr_mse_ts(t)=sum_mse_ts/sum_nmse/2; %unified power
end



beta_t=[0.2,0.3,0.4,0.5,0.6,0.7,0.8];
plot(beta_t,avr_mse_es,'>-');
hold on;
plot(beta_t,avr_mse_t,'o--');
hold on;
plot(beta_t,avr_mse_r,'x--');
hold on;
plot(beta_t,avr_mse_ts,'s-');
hold on;
xlim([0.2 0.8])
%ylim([0 4]*1e-3)
xlabel('\beta_t');
ylabel('NMSE');
%grid on;
legend('ES','ES: T user','ES: R user','TS');


