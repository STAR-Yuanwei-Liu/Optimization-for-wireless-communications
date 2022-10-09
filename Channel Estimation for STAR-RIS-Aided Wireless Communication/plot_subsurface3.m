clear all;
clc;
Nx=5;
Ny=16;
M0=Nx*Ny;
pho=0.001;
d_as=50; %%AP-STAR
d_sut=5;   %STAR-T user
d_aut=norm([0,0]-[54,3]);  %AP-T user
d_sur=5;   %STAR-R user
d_aur=norm([0,0]-[46,-3]);  %AP-R user
alpha_as=2.2;
alpha_su=2.8;
alpha_au=3.5;
L_as=pho*d_as^-alpha_as;
L_sut=pho*d_sut^-alpha_su;
L_aut=pho*d_aut^-alpha_au;
L_sur=pho*d_sur^-alpha_su;
L_aur=pho*d_aur^-alpha_au;

SNR=140;
kappa=10; %6 dB
aoa_theta_as=pi/2;
aod_theta_as=0;
aoa_phi_as=0;
aoa_theta_sut=atan(-4/3);
aod_theta_sut=pi/2-atan(-4/3);
aod_phi_sut=0;
aoa_theta_sur=atan(4/3);
aod_theta_sur=pi/2-atan(4/3);
aod_phi_sur=0;


avr_mse_es=zeros(1,8);
avr_mse_es_prac=zeros(1,8);
avr_mse_es_2p=zeros(1,8);
avr_mse_t=zeros(1,8);
avr_mse_r=zeros(1,8);
avr_mse_ts=zeros(1,8);
avr_mse_on_off=zeros(1,8);
iter=100;

for t=1:8
M=5*t; %%dB P=1W -170dBm/Hz   
sub=floor(M0/M);
left=M0-sub*(M-1);
sum_mse_es=0;
sum_mse_t=0;
sum_mse_r=0;
sum_mse_es_prac=0;
sum_mse_es_2p=0;
sum_mse_ts=0;
sum_nmse=0;
sum_mse_on_off=0;

for k=1:iter
g=zeros(M0,1);
rt=zeros(M0,1);
rr=zeros(M0,1);
for m= 0:Nx-1
    for n= 0:Ny-1
        %y(m*Nx+n+1) = exp( 1i* pi*( m*sin(phi)*sin(theta) + n*cos(theta)));
%         g(m*Ny+n+1,1) =sqrt(L_as)*(sqrt(kappa/(1+kappa))*exp(1i* pi*( m*sin(aoa_phi_as)*sin(aoa_theta_as)+n*cos(aoa_theta_as)))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
%         rt(m*Ny+n+1,1) =sqrt(L_sut)*(sqrt(kappa/(1+kappa))*exp(1i* pi*(m*sin(aod_phi_sut)*sin(aod_theta_sut)+n*cos(aod_theta_sut)))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
%         rr(m*Ny+n+1,1) =sqrt(L_sur)*(sqrt(kappa/(1+kappa))*exp(1i* pi*(m*sin(aod_phi_sur)*sin(aod_theta_sur)+n*cos(aod_theta_sur)))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
        g(m*Ny+n+1,1) =sqrt(L_as)*(sqrt(kappa/(1+kappa))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
        rt(m*Ny+n+1,1) =sqrt(L_sut)*(sqrt(kappa/(1+kappa))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
        rr(m*Ny+n+1,1) =sqrt(L_sur)*(sqrt(kappa/(1+kappa))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));    
    end
end
ht=sqrt(L_aut)*(sqrt(kappa/(1+kappa))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
hr=sqrt(L_aur)*(sqrt(kappa/(1+kappa))+sqrt(1/(1+kappa))*1/sqrt(2)*(randn(1)+1j*randn(1)));
qt=diag(g)*rt;
qr=diag(g)*rr;
qt_group=zeros(M,1);
qr_group=zeros(M,1);


for i=1:M-1
    for j=1:sub
    qt_group(i)=qt_group(i)+qt(j+sub*(i-1));
    qr_group(i)=qr_group(i)+qr(j+sub*(i-1));
    end
end
%if left~=M
for j=1:left
    qt_group(M)=qt_group(M)+qt(j+sub*(M-1));
    qr_group(M)=qr_group(M)+qr(j+sub*(M-1));
end



%% ideal 
x=[ht;qt_group;hr;qr_group];
sigma_2=10^(-SNR/10);
noise=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2*M+2,1)+1j*randn(2*M+2,1));
H_temp=dftmtx(2*M+2);
H=H_temp;
H(:,2:1+M)=H(:,2:1+M)*sqrt(0.5);
H(:,3+M:2*M+2)=H(:,3+M:2*M+2)*sqrt(0.5);

%end
y=H*x+noise;

y_g_ideal=pinv(H)*y;
%error=pinv(H)*y-x; %%H'/T!=inv(H)
error=y_g_ideal-[ht;qt_group;hr;qr_group];
sum_mse_es=sum_mse_es+error'*error;


%error=H'/T*y-x; %%H'/T!=inv(H)
%% practical phase shift
% current mathod

noise2=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2*M+2,1)+1j*randn(2*M+2,1));
H_prac=H;
for i=1:2*M+2
    H_prac(i,M+2)=(-1)^(i+1)*1i;
end
y_g_prac=pinv(H_prac)*(H_prac*x+noise2);

error_prac=y_g_prac-[ht;qt_group;hr;qr_group];
sum_mse_es_prac=sum_mse_es_prac+error_prac'*error_prac;

%% Two phase

noise1=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2,1)+1j*randn(2,1));
y_1p=[1,1;1,-1]*[ht;hr]+noise1;
y_1p_es=inv([1,1;1,-1])*y_1p;
error_1p=y_1p_es-[ht;hr];


H_2p=dftmtx(2*M)*sqrt(0.5);
noise2=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2*M,1)+1j*randn(2*M,1));
y_2p=H_2p*[qt_group;qr_group]+noise2+ones(2*M,1)*ht+ones(2*M,1)*hr*1i;
y_2p_es_g=pinv(H_2p)*(y_2p-ones(2*M,1)*y_1p_es(1)-ones(2*M,1)*1i*y_1p_es(2));

%error_2p=y_2p_es-[qt_group;qr_group];
error_2p=y_2p_es_g-[qt_group;qr_group];
sum_mse_es_2p=sum_mse_es_2p+error_2p'*error_2p+error_1p'*error_1p;



%% Time swithing
% H_ts_temp=dftmtx(T/2);
% H_ts=H_ts_temp(:,1:M+1);
% noise_t_ts=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(T/2,1)+1j*randn(T/2,1));
% y_t=H_ts*[ht;qt]+noise_t_ts;
% y_t_ts=pinv(H_ts)*y_t;
% noise_r_ts=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(T/2,1)+1j*randn(T/2,1));
% y_r=H_ts*[hr;qr]+noise_r_ts;
% y_r_ts=pinv(H_ts)*y_r;
% error_t=y_t_ts-[ht;qt];
% error_r=y_r_ts-[hr;qr];
% sum_mse_ts=sum_mse_ts+real(trace(error_t*error_t'))+real(trace(error_r*error_r'));

%grouping

H_ts=dftmtx(M+1);

noise_t_ts=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(M+1,1)+1j*randn(M+1,1));
y_t=H_ts*[ht;qt_group]+noise_t_ts;
y_g_ts=pinv(H_ts)*y_t;
% error_t=y_t_ts-[ht;qt_group];
%error_r=y_r_ts-[hr;qr_group];
error_ts=y_g_ts-[ht;qt_group];
sum_mse_ts=sum_mse_ts+2*(error_ts'*error_ts);





%% Time swithing On_OFF

noise_on_off=sqrt(10^(-SNR/10))*1/sqrt(2)*(randn(2*M,1)+1j*randn(2*M,1));
y_g_on_off=[qt_group;qr_group]+noise_on_off;

error_on_off=y_g_on_off-[qt_group;qr_group];

sum_mse_on_off=sum_mse_on_off+error_on_off'*error_on_off++error_1p'*error_1p;


%% NMSE
%sum_nmse=sum_nmse+[ht;hr;qt_group;qr_group]'*[ht;hr;qt_group;qr_group];
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


time=5:5:40;
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
xlabel('Number of STAR-RIS subsurfaces: M');
ylabel('NMSE');
grid on;
legend('TS, ON/OFF','TS','ES, two phase','ES, ideal phase shift','ES, practical phase shift');


