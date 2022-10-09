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
avr_mse_es=zeros(1,6);
avr_mse_es_prac=zeros(1,6);
avr_mse_es_2p=zeros(1,6);
avr_mse_t=zeros(1,6);
avr_mse_r=zeros(1,6);
avr_mse_ts=zeros(1,6);
avr_mse_on_off=zeros(1,6);
T=68;  
iter=100;

for t=1:9
SNR=110+5*(t-1); %%dB P=1W -170dBm/Hz   

sum_mse_es=0;
sum_mse_t=0;
sum_mse_r=0;
sum_mse_es_prac=0;
sum_mse_es_2p=0;
sum_mse_ts=0;
sum_nmse=0;
sum_mse_on_off=0;

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


%error=H'/T*y-x; %%H'/T!=inv(H)
%% practical phase shift
% past method
% x1=ones(T,1);
% x2=x1;
% x2(2:2:T)=-1;
% c=randn(T,M)+randn(T,M)*1i;
% for i=1:T
%     for j=1:M
%     temp=randn(1);
%     c(i,j)=c(i,j)/norm(c(i,j));
%     if temp<0.5
%     d(i,j)=c(i,j)*1i*sqrt(beta_t); 
%     else
%     d(i,j)=c(i,j)*-1i*sqrt(beta_r);    
%     end
%     if mod(i,2)==0
%     d(i,j)=d(i,j)*(-1);
%     end
%     end
% end
% H_prac=[x1,x2,c,d];

% current mathod
H_temp=dftmtx(T);
H_prac=H_temp(:,1:2*M+2);
x_prac=[ht;qt;hr;qr];
H_prac(:,M+2)=[ones(T/2-1,1)*1i;ones(T/2+1,1)*-1i];
H_prac(:,2:M+1)=H_prac(:,2:M+1)*sqrt(beta_t);
H_prac(:,M+3:2*M+2)=H_prac(:,M+3:2*M+2)*sqrt(beta_r);

y_prac=H_prac*x_prac+noise;
error_prac=pinv(H_prac)*y_prac-x_prac; %%H'/T!=inv(H)
mse_matrix_prac=error_prac*error_prac';
sum_mse_es_prac=sum_mse_es_prac+real(trace(mse_matrix_prac));

%% Two phase

% noise1=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2,1)+1j*randn(2,1));
% y_1p=[1,1;1,-1]*[ht;hr]+noise1;
% y_1p_es=inv([1,1;1,-1])*y_1p;
% error_1p=inv([1,1;1,-1])*y_1p-[ht,hr]';
% mse_1p=error_1p'*error_1p;
% %[ht_es hr_es]=
% H_2p_temp=dftmtx(T-2)*sqrt(0.5);
% H_2p=H_2p_temp(:,1:2*M);
% noise2=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(T-2,1)+1j*randn(T-2,1));
% y_2p=H_2p*[qt;qr]+noise2+ones(T-2,1)*ht+ones(T-2,1)*hr*1i;
% y_2p_es=pinv(H_2p)*(y_2p-ones(T-2,1)*y_1p_es(1)-ones(T-2,1)*1i*y_1p_es(2));
% error_2p=y_2p_es-[qt;qr];
% mse_matrix_2p=error_2p*error_2p';
% sum_mse_es_2p=sum_mse_es_2p+real(trace(mse_matrix_2p));

% optimized
t1=(T-2*M)/2;
t2=T-t1;
noise1=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(t1,1)+1j*randn(t1,1));
x_1p=[ones(t1,1),[ones(t1/2,1);ones(t1/2,1)*-1]];
y_1p=x_1p*[ht;hr]+noise1;
y_1p_es=pinv(x_1p)*y_1p; %estimate
error_1p=pinv(x_1p)*y_1p-[ht,hr]';
mse_1p=error_1p'*error_1p;
H_2p_temp=dftmtx(t2)*sqrt(0.5);
H_2p=H_2p_temp(:,1:2*M);
noise2=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(t2,1)+1j*randn(t2,1));
y_2p=H_2p*[qt;qr]+noise2+ones(t2,1)*ht+ones(t2,1)*hr*1i;
y_2p_es=pinv(H_2p)*(y_2p-ones(t2,1)*y_1p_es(1)-ones(t2,1)*1i*y_1p_es(2));
error_2p=y_2p_es-[qt;qr];
mse_matrix_2p=error_2p*error_2p';
sum_mse_es_2p=sum_mse_es_2p+real(trace(mse_matrix_2p));


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



%% Time swithing On_OFF

noise_on_off=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2*M+2,1)+1j*randn(2*M+2,1));
sum_mse_on_off=sum_mse_on_off+real(trace(noise_on_off*noise_on_off'));


%% NMSE
sum_nmse=sum_nmse+[ht;hr;qt;qr]'*[ht;hr;qt;qr];
end
avr_mse_es(t)=sum_mse_es/sum_nmse;
avr_mse_t(t)=sum_mse_t/sum_nmse;
avr_mse_r(t)=sum_mse_r/sum_nmse;
avr_mse_es_prac(t)=sum_mse_es_prac/sum_nmse;
avr_mse_es_2p(t)=sum_mse_es_2p/sum_nmse;
avr_mse_ts(t)=sum_mse_ts/sum_nmse/2; %unified power
avr_mse_on_off(t)=sum_mse_on_off/sum_nmse; %unified power
end


time=0:5:40;
ideal=avr_mse_es;
ts=avr_mse_ts;
two_phase=avr_mse_es_2p;
prac=avr_mse_es_prac;
on_off=avr_mse_on_off;

%plot(time,on_off,'<-.');
semilogy(time,on_off,'<-.');
hold on;

% plot(time,ts,'s-');
semilogy(time,ts,'s-.');
hold on;


%plot(time,two_phase,'>--');
semilogy(time,two_phase,'>-');
hold on;

%plot(time,ideal,'o-');
semilogy(time,ideal,'o-');
hold on;

%plot(time,prac,'x-.');
semilogy(time,prac,'x-');
hold on;


set(gca,'XTick',time);

%ylim([0 4]*1e-3)
xlabel('Total transmit power (dBm)');
ylabel('NMSE');
grid on;
legend('TS, ON/OFF','TS','ES, two phase','ES, ideal phase shift','ES, practical phase shift');

