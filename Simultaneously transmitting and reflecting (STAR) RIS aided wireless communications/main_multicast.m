%Simulation codes for the multicast scenario
%Date: 01/10/2022
%Author: Xidong Mu

clear all;
tic;
% The number of users is K, the number of antennas at the BS is N, the
% number of reflection elements is M, the user is single-antenna
K=2;
N=4;
M=10;
PL=12;%noise power
K_BS=2;%The Rician factor of the BS-STAR channel, 3dB
a_BS=2.2;%The path-loss exponent of the BS-STAR channel
K_SU=2;%The Rician factor of the STAR-user channel, 3dB
a_SU=2.2;%The path-loss exponent of the STAR-use channel
% Initialize the locations of the BS, the STAR, and K users
location_BS=[-50;0;0];
location_STAR=[0;0;0];
%The users are uniformly distributed in a circle centered at the STAR with the radius of 3 m
radius(1)=3;
radius(2)=3;
maxOut=10;%Maximum number of outer iterations
maxIn=10;%Maximum number of inner iterations
c=10;
for k=1:K
    QoS(k)=10^(10/10);
end
for combination=1:1:5% ranging M from 10 to 50

    for realization=1:1 %channel realizations

        para = para_init();
        para.STAR_size = [combination*2,5]; % elements at STAR

        para.N = para.STAR_size(1)*para.STAR_size(2);
        %Generate channels
        [BS_array, STAR_array] = generate_arrays(para);
        [G, r] = generate_channel(para, BS_array, STAR_array);
        %The cascaded channel from the BS to the user with the aid of the STAR is
        H=[];
        for k=1:K
            H(:,:,k)=diag(r(:,k)',0)*G;
        end
        M=para.N; %the number of STAR elements
        N=4;%the number of BS antennas

        %Initialize the reflection and tranmission reflection matrix \Theta with random phase shifts from [0,2*pi), which needs to be optimized
        Theta=[];
        Q_r0=[];
        for k=1:K
            for m=1:M
                u(m,k)=exp(1i*pi*(2*rand(1,1)));
            end
            Theta(:,:,k)=diag(u(:,k),0);%M x M reflection matrix
            Q_r0(:,:,k)=u(:,k)*u(:,k)';
        end



        %%%%%%%%%%%%%%%% Energy splitting
        %Initialized a random passive beamforming and find feasible active beamforming
        Q_r_ES0=1/2*Q_r0;%Initialize eaqual energy splitting
        cvx_begin
        variable W_ES0(N,N) Hermitian semidefinite
        expressions p_ES0
        p_ES0=real(trace(W_ES0));
        minimize sum(p_ES0)
        subject to
        for k=1:K
            QoS(k)-real(trace(Q_r_ES0(:,:,k)*H(:,:,k)*W_ES0*H(:,:,k)'))<=0;
        end
        cvx_end
        if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
            consum_power_GES(realization,combination)=100;
            break;
        end
        %%%%%%%%%%%%%%%% joint active and passive beamforming optimization
        Q_r_ES=[];
        Q_r_ES=Q_r_ES0;%Initialize passive beamforming
        W_r_ES=W_ES0;%Initialize active beamforming
        for outer=1:maxOut
            if(outer==1)
                factor_ES=0;
            else
                factor_ES=0.0001*c^(outer-1);
            end
            for inter=1:maxIn
                uu_ES=[];
                for k=1:K
                    [E_ES,D_ES]=eig(Q_r_ES(:,:,k));
                    value_ES=diag(D_ES);
                    [value_m_ES,value_index_ES]=max(value_ES);
                    uu_ES(:,:,k)=E_ES(:,value_index_ES)*E_ES(:,value_index_ES)';
                end
                cvx_begin
                variable beta_ES(M,K)
                variable W_ES(N,N) Hermitian semidefinite
                variable Q_ES(M,M,K) Hermitian semidefinite
                expressions p_ES penlty_ES(K)
                for k=1:K
                    penlty_ES(k)=real(trace(Q_ES(:,:,k)))-real(trace(uu_ES(:,:,k)*(Q_ES(:,:,k)-Q_r_ES(:,:,k))))-norm(Q_r_ES(:,:,k),2);
                end
                p_ES=real(trace(W_ES));
                minimize sum(p_ES)+factor_ES*sum(penlty_ES)
                subject to
                for k=1:K
                    QoS(k)+1/2*pow_pos(norm(Q_ES(:,:,k)-H(:,:,k)*W_ES*H(:,:,k)','fro'),2)+1/2*pow_pos(norm(Q_r_ES(:,:,k),'fro'),2)-1*real(trace(Q_r_ES(:,:,k)'*Q_ES(:,:,k)))+1/2*pow_pos(norm(H(:,:,k)*W_r_ES*H(:,:,k)','fro'),2)-1*real(trace((H(:,:,k)'*H(:,:,k)*W_r_ES*H(:,:,k)'*H(:,:,k))'*W_ES))<=0;
                end
                for k=1:K
                    for m=1:M
                        Q_ES(m,m,k) == beta_ES(m,k);
                    end
                end
                for m=1:M
                    sum(beta_ES(m,:))==1;
                end
                beta_ES>=0;
                beta_ES<=1;
                cvx_end
                if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
                    consum_power_GES(realization,combination)=100;
                    break;
                end
                Q_r_ES=Q_ES;
                W_r_ES=W_ES;
                %max(penlty_ES)
                if(max(penlty_ES)<=10^(-7))
                    consum_power_ES(realization,combination)=sum(p_ES);
                    break;
                end
            end
            if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
                consum_power_GES(realization,combination)=100;
                break;
            end
            if(max(penlty_ES)<=10^(-7))
                consum_power_ES(realization,combination)=sum(p_ES);
                break;
            end
        end
        consum_power_GES(realization,combination)=sum(p_ES);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%% Mode switching
        for i=1:1
            Q_r_MS0=1/2*Q_r0;
            cvx_begin
            variable W_MS0(N,N) Hermitian semidefinite
            expressions p_MS0
            p_MS0=real(trace(W_MS0));
            minimize sum(p_MS0)
            subject to
            for k=1:K
                QoS(k)-real(trace(Q_r_MS0(:,:,k)*H(:,:,k)*W_MS0*H(:,:,k)'))<=0;
            end
            cvx_end
            if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
                consum_power_MS(realization,combination)=100;
                break;
            end
            %%%%%%%%%%%%%%% joint optimization rank one
            for k=1:K
                beta_MS_r(1:M,k)=1/2;
            end
            Q_r_MS=Q_r_MS0;
            W_r_MS=W_MS0;
            for outer=1:maxOut
                if(outer==1)
                    factor_MS=0;
                else
                    factor_MS=0.0001*c^(outer-1);
                end
                for inter=1:maxIn
                    uu_MS=[];
                    for k=1:K
                        [E_MS,D_MS]=eig(Q_r_MS(:,:,k));
                        value_MS=diag(D_MS);
                        [value_m_MS,value_index_MS]=max(value_MS);
                        uu_MS(:,:,k)=E_MS(:,value_index_MS)*E_MS(:,value_index_MS)';
                    end
                    cvx_begin
                    variable beta_MS(M,K)
                    variable W_MS(N,N) Hermitian semidefinite
                    variable Q_MS(M,M,K) Hermitian semidefinite
                    expressions p_MS penlty_MS(K) penlty_BMS(M)
                    for m=1:M
                        penlty_BMS(m)=(-2*beta_MS_r(m,1)+1)*beta_MS(m,1)+beta_MS_r(m,1)^2+(-2*beta_MS_r(m,2)+1)*beta_MS(m,2)+beta_MS_r(m,2)^2;
                    end
                    for k=1:K
                        penlty_MS(k)=real(trace(Q_MS(:,:,k)))-real(trace(uu_MS(:,:,k)*(Q_MS(:,:,k)-Q_r_MS(:,:,k))))-norm(Q_r_MS(:,:,k),2);
                    end
                    p_MS=real(trace(W_MS));
                    minimize sum(p_MS)+factor_MS*sum(penlty_MS)+factor_MS*sum(penlty_BMS)
                    subject to
                    for k=1:K
                        QoS(k)+1/2*pow_pos(norm(Q_MS(:,:,k)-H(:,:,k)*W_MS*H(:,:,k)','fro'),2)+1/2*pow_pos(norm(Q_r_MS(:,:,k),'fro'),2)-1*real(trace(Q_r_MS(:,:,k)'*Q_MS(:,:,k)))+1/2*pow_pos(norm(H(:,:,k)*W_r_MS*H(:,:,k)','fro'),2)-1*real(trace((H(:,:,k)'*H(:,:,k)*W_r_MS*H(:,:,k)'*H(:,:,k))'*W_MS))<=0;
                    end
                    for k=1:K
                        for m=1:M
                            Q_MS(m,m,k) == beta_MS(m,k);
                        end
                    end
                    for m=1:M
                        sum(beta_MS(m,:))==1;
                    end
                    beta_MS>=0;
                    beta_MS<=1;
                    cvx_end
                    if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
                        consum_power_MS(realization,combination)=100;
                        break;
                    end
                    beta_MS_r=beta_MS;
                    Q_r_MS=Q_MS;
                    W_r_MS=W_MS;
                    %max(penlty_MS)
                    %max(penlty_BMS)
                    if(max(max(penlty_MS),max(penlty_BMS))<=10^(-7))
                        consum_power_MS(realization,combination)=sum(p_MS);
                        break;
                    end
                end
                if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
                    consum_power_MS(realization,combination)=100;
                    break;
                end
                if(max(max(penlty_MS),max(penlty_BMS))<=10^(-7))
                    consum_power_MS(realization,combination)=sum(p_MS);
                    break;
                end
            end
            consum_power_MS(realization,combination)=sum(p_MS);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        %%%%%%%%%%%%%%%%%%%%%%Time Switching
        %%%%%%%%%%%%%%%%%% find the optimal individual W and Q
        Q_r_TS=Q_r0;
        Q_TS=[];
        for k=1:K
            factor_TS=0;
            for iter=1:10
                uu_TS=[];
                [E_TS,D_TS]=eig(Q_r_TS(:,:,k));
                value_TS=diag(D_TS);
                [value_m_TS,value_index_TS]=max(value_TS);
                uu_TS(:,:,k)=E_TS(:,value_index_TS)*E_TS(:,value_index_TS)';
                cvx_begin
                variable Q_TS_t(M,M) Hermitian semidefinite
                expression penlty_TS
                penlty_TS=real(trace(Q_TS_t))-real(trace(uu_TS(:,:,k)*(Q_TS_t-Q_r_TS(:,:,k))))-norm(Q_r_TS(:,:,k),2);
                minimize -real(trace(Q_TS_t*H(:,:,k)*H(:,:,k)'))+factor_TS*penlty_TS
                subject to
                for m=1:M
                    Q_TS_t(m,m) == 1;
                end
                cvx_end
                Q_r_TS(:,:,k)=Q_TS_t;
                factor_TS=(factor_TS+0.5);
                if(penlty_TS<=10^(-7))
                    break;
                end
            end
            Q_TS(:,:,k)=Q_TS_t;
            channel_gain_TS(1,k)=real(trace(Q_TS(:,:,k)*H(:,:,k)*H(:,:,k)'));
        end
        %%%%%%%%%%%%%%%%%% find the optimal power
        cvx_begin
        variable p_TS(K)
        variable lambda_TS(K)
        minimize sum(p_TS)
        subject to
        for k=1:K
            -rel_entr(lambda_TS(k),lambda_TS(k)+channel_gain_TS(1,k)*p_TS(k))/log(2) >=log2(1+QoS(k));
        end
        sum(lambda_TS)==1;
        lambda_TS>=0;
        cvx_end
        if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
            consum_power_TS(realization,combination)=100;
            break;
        end
        consum_power_TS(realization,combination)=sum(p_TS);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    end



end

toc;






