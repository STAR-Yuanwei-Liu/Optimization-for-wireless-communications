clear all;
tic;
% The number of users is K, the number of antennas at the BS is N, the
% number of reflection elements is M, the user is single-antenna
K=2;
N=1;
M=5;
PL=11;
K_BS=2;%BS-STAR Rician factor
a_BS=2.2;%BS-STAR path loss exponent
K_SU=2;%STAR-user Rician factor
a_SU=2.2;%STAR-user path loss exponent
a_BU=3.5;%BS-user path loss exponent
% Initialize the locations of the BS, the STAR, and K users
location_BS=[-50;0;0];
location_STAR=[0;0;0];
%The users are uniformly distributed in a circle centered at the STAR with the radius of 3 m
radius(1)=3;
radius(2)=3;

iteration=3;%Maximum number of iterations

for combination=1:1:5
    M=10*(combination);
    for realization=1:1
        %generate user's locations
        a(1)=1/2*pi+pi*rand(1,1);
        location_user(1,1)=radius(1)*cos(a(1));
        location_user(2,1)=radius(1)*sin(a(1));
        location_user(3,1)=0;
        a(2)=-1/2*pi+pi*rand(1,1);
        location_user(1,2)=radius(2)*cos(a(2));
        location_user(2,2)=radius(2)*sin(a(2));
        location_user(3,2)=0;
        %generate BS-STAR Rician channel
        distance_BS_STAR=norm(location_STAR-location_BS);
        path_loss_BS_STAR=sqrt(10^(-3)*distance_BS_STAR^(-a_BS));
        AoD_BS_STAR=atan((location_BS(1,1)-location_STAR(1,1))/(abs(location_BS(2,1)-location_STAR(2,1))));
        AoA_BS_STAR=1/2*pi-AoD_BS_STAR;
        array_response_BS_STAR_AOA=[];
        G_LoS=[];
        G_NLoS=[];
        g=[];
        array_response_STAR_user=[];
        r_LoS=[];
        r_NLoS=[];
        r=[];
        H=[];
        u=[];
        Theta=[];
        Q_r0=[];
        for m=1:M
            array_response_BS_STAR_AOA(m,1)=exp(1i*pi*(m-1)*sin(AoA_BS_STAR));
        end
        G_LoS=array_response_BS_STAR_AOA;
        G_NLoS=1/sqrt(2)*(randn(M,1)+1i*randn(M,1));
        g=path_loss_BS_STAR*(sqrt(K_BS/(1+K_BS))*G_LoS+sqrt(1/(1+K_BS))*G_NLoS);
        %generate STAR-user Rician channel and BS-user Rayleigh channel
        for k=1:K
            distance_STAR_user(k)=norm(location_user(:,k)-location_STAR);
            path_loss_STAR_user(k)=sqrt(10^(-3)*distance_STAR_user(k)^(-a_SU));
            AoD_STAR_user(k)=atan((location_user(2,k)-location_STAR(2,1))/(abs(location_user(1,1)-location_STAR(1,1))));
            distance_BS_user(k)=norm(location_user(:,k)-location_BS);
            path_loss_BS_user(k)=sqrt(10^(-3)*distance_BS_user(k)^(-a_BU));
        end
        for k=1:K
            for m=1:M
                array_response_STAR_user(m,k)=exp(1i*pi*(m-1)*sin(AoD_STAR_user(k)));
            end
            r_LoS(:,k)=array_response_STAR_user(:,k);
        end
        for k=1:K
            r_NLoS(:,k)=1/sqrt(2)*(randn(M,1)+1i*randn(M,1));
        end
        for k=1:K
            r(:,k)=path_loss_STAR_user(k)*(sqrt(K_SU/(1+K_SU))*r_LoS(:,k)+sqrt(1/(1+K_SU))*r_NLoS(:,k));%M x 1 matrix
            d(1,k)=path_loss_BS_user(k)*1/sqrt(2)*(randn(1,1)+1i*randn(1,1));
        end

        %The cascaded channel from the BS to the user with the aid of the STAR is
        z=[];
        for k=1:K
            z(:,k)=sqrt(10^(PL))*[r(:,k)'*diag(g,0) d(1,k)]';
        end

        %Initialize transmission and reflection coefficient
        for m=1:M
            ang(m,1)=rand(1,1)*2*pi;
            beta_r(m,1)=sqrt(1/2);
            beta_r(m,2)=sqrt(1/2);
            theta_r(m,1)=exp(ang(m,1)*1i);
            theta_r(m,2)=exp((ang(m,1)+1/2*pi)*1i);
        end
        theta_r(M+1,1)=1;
        theta_r(M+1,2)=1;
        beta_r(M+1,1)=1;
        beta_r(M+1,2)=1;

        %Initialize rate constraints
        rate(1)=2;
        rate(2)=2;
        for k=1:K
            gamma(k)=2^rate(k)-1;
            gammaO(k)=0.5*(2^(2*rate(k))-1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% element-wise optimization in NOMA with coupled phase shift
        for order=1:K
            if(order==1)
                %decoding order 1
                beta_1=beta_r;
                theta_1=theta_r;
                for inter=1:iteration
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%phase-shift optimization
                    for m=1:M
                        for k=1:K
                            f1(m,k)=norm((z(:,k)'*(beta_1(:,k).*theta_1(:,k))-z(m,k)'*(beta_1(m,k)*theta_1(m,k))),2)^2+norm(beta_1(m,k)*z(m,k)',2)^2;
                            f2(m,k)=z(:,k)'*(beta_1(:,k).*theta_1(:,k))-z(m,k)'*(beta_1(m,k)*theta_1(m,k));
                        end
                        for i=1:K
                            if(i==1)
                                cvx_begin
                                variable T(K)
                                variable theta1 complex
                                expressions gain1(K) p(K)
                                gain1(1)=f1(m,1)+2*real(z(m,1)'*(beta_1(m,1)*theta1)*f2(m,1)');
                                gain1(2)=f1(m,2)+2*real(z(m,2)'*(beta_1(m,2)*1i*theta1)*f2(m,2)');
                                p(1)=gamma(1)*inv_pos(T(1));
                                p(2)=gamma(2)*inv_pos(T(2))+gamma(1)*gamma(2)*inv_pos(T(1));
                                minimize sum(p)
                                subject to
                                for k=1:K
                                    gain1(k)>=T(k);
                                end
                                norm(theta1)<=1;
                                cvx_end
                                value(i)=sum(p);
                                norm(theta1)
                            else
                                cvx_begin
                                variable T(K)
                                variable theta2 complex
                                expressions gain2(K) p(K)
                                gain2(1)=f1(m,1)+2*real(z(m,1)'*(beta_1(m,1)*theta2)*f2(m,1)');
                                gain2(2)=f1(m,2)+2*real(z(m,2)'*(beta_1(m,2)*(-1i)*theta2)*f2(m,2)');
                                p(1)=gamma(1)*inv_pos(T(1));
                                p(2)=gamma(2)*inv_pos(T(2))+gamma(1)*gamma(2)*inv_pos(T(1));
                                minimize sum(p)
                                subject to
                                for k=1:K
                                    gain2(k)>=T(k);
                                end
                                norm(theta2)<=1;
                                cvx_end
                                value(i)=sum(p);
                            end
                        end
                        if(value(1)<value(2))
                            theta_1(m,1)=theta1;
                            theta_1(m,2)=1i*theta1;
                        else
                            theta_1(m,1)=theta2;
                            theta_1(m,2)=-1i*theta2;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%amplitude optimization
                    w=[];
                    for k=1:K
                        w(:,k)=sqrt(10^(PL))*[r(:,k)'*diag(theta_1(1:M,k),0)*diag(g,0) d(1,k)]';
                    end
                    for m=1:M
                        for k=1:K
                            l1(m,k)=norm((w(:,k)'*beta_1(:,k)-w(m,k)'*beta_1(m,k)),2)^2;%+norm(beta_11(m,k)*z(m,k)',2)^2;
                            l2(m,k)=w(:,k)'*(beta_1(:,k))-w(m,k)'*(beta_1(m,k));
                            l3(m,k)=2*real(w(m,k)'*l2(m,k)');
                        end
                        cvx_begin
                        variable T(K)
                        variable beta_t(K)
                        expressions gain_b(K) p(K)
                        for k=1:K
                            if(l3(m,k)>=0)
                                gain_b(k)=l1(m,k)+norm(w(m,k)',2)^2*beta_t(k)+l3(m,k)*sqrt(beta_t(k));
                            else
                                gain_b(k)=l1(m,k)+norm(w(m,k)',2)^2*beta_t(k)+l3(m,k)*(0.5*beta_1(m,k)^(-0.5)*(beta_t(k)-beta_1(m,k)));
                            end
                        end
                        p(1)=gamma(1)*inv_pos(T(1));
                        p(2)=gamma(2)*inv_pos(T(2))+gamma(1)*gamma(2)*inv_pos(T(1));
                        minimize sum(p)
                        subject to
                        for k=1:K
                            gain_b(k) >= T(k);
                        end
                        sum(beta_t)==1;
                        beta_t>=0;
                        cvx_end
                        for k=1:K
                            beta_1(m,k)=sqrt(beta_t(k));
                        end
                    end
                end
                result_temp(order)=sum(p);
            end
            if(order==2)
                %decoding order 2
                beta_2=beta_r;
                theta_2=theta_r;
                for inter=1:iteration
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%phase-shift optimization
                    for m=1:M
                        for k=1:K
                            f1(m,k)=norm((z(:,k)'*(beta_2(:,k).*theta_2(:,k))-z(m,k)'*(beta_2(m,k)*theta_2(m,k))),2)^2+norm(beta_2(m,k)*z(m,k)',2)^2;
                            f2(m,k)=z(:,k)'*(beta_2(:,k).*theta_2(:,k))-z(m,k)'*(beta_2(m,k)*theta_2(m,k));
                        end
                        for i=1:K
                            if(i==1)
                                cvx_begin
                                variable T(K)
                                variable theta1 complex
                                expressions gain1(K) p(K)
                                gain1(1)=f1(m,1)+2*real(z(m,1)'*(beta_2(m,1)*theta1)*f2(m,1)');
                                gain1(2)=f1(m,2)+2*real(z(m,2)'*(beta_2(m,2)*1i*theta1)*f2(m,2)');
                                p(2)=gamma(2)*inv_pos(T(2));
                                p(1)=gamma(1)*inv_pos(T(1))+gamma(1)*gamma(2)*inv_pos(T(2));
                                minimize sum(p)
                                subject to
                                for k=1:K
                                    gain1(k)>=T(k);
                                end
                                norm(theta1)<=1;
                                cvx_end
                                value(i)=sum(p);
                                norm(theta1)
                            else
                                cvx_begin
                                variable T(K)
                                variable theta2 complex
                                expressions gain2(K) p(K)
                                gain2(1)=f1(m,1)+2*real(z(m,1)'*(beta_2(m,1)*theta2)*f2(m,1)');
                                gain2(2)=f1(m,2)+2*real(z(m,2)'*(beta_2(m,2)*(-1i)*theta2)*f2(m,2)');
                                p(2)=gamma(2)*inv_pos(T(2));
                                p(1)=gamma(1)*inv_pos(T(1))+gamma(1)*gamma(2)*inv_pos(T(2));
                                minimize sum(p)
                                subject to
                                for k=1:K
                                    gain2(k)>=T(k);
                                end
                                norm(theta2)<=1;
                                cvx_end
                                norm(theta2)
                                value(i)=sum(p);
                            end
                        end
                        if(value(1)<value(2))
                            theta_2(m,1)=theta1;
                            theta_2(m,2)=1i*theta1;
                        else
                            theta_2(m,1)=theta2;
                            theta_2(m,2)=-1i*theta2;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%amplitude optimization
                    for k=1:K
                        w(:,k)=sqrt(10^(PL))*[r(:,k)'*diag(theta_2(1:M,k),0)*diag(g,0) d(1,k)]';
                    end
                    for m=1:M

                        for k=1:K
                            l1(m,k)=norm((w(:,k)'*beta_2(:,k)-w(m,k)'*beta_2(m,k)),2)^2;%+norm(beta_11(m,k)*z(m,k)',2)^2;
                            l2(m,k)=w(:,k)'*(beta_2(:,k))-w(m,k)'*(beta_2(m,k));
                            l3(m,k)=2*real(w(m,k)'*l2(m,k)');
                        end
                        cvx_begin
                        variable T(K)
                        variable beta_t(K)
                        expressions gain_b(K) p(K)
                        for k=1:K
                            if(l3(m,k)>=0)
                                gain_b(k)=l1(m,k)+norm(w(m,k)',2)^2*beta_t(k)+l3(m,k)*sqrt(beta_t(k));
                            else
                                gain_b(k)=l1(m,k)+norm(w(m,k)',2)^2*beta_t(k)+l3(m,k)*(0.5*beta_2(m,k)^(-0.5)*(beta_t(k)-beta_2(m,k)));
                            end
                        end
                        p(2)=gamma(2)*inv_pos(T(2));
                        p(1)=gamma(1)*inv_pos(T(1))+gamma(1)*gamma(2)*inv_pos(T(2));
                        minimize sum(p)
                        subject to
                        for k=1:K
                            gain_b(k) >= T(k);
                        end
                        sum(beta_t)==1;
                        beta_t>=0;
                        cvx_end
                        for k=1:K
                            beta_2(m,k)=sqrt(beta_t(k));
                        end
                    end
                end
                result_temp(order)=sum(p);
            end
        end
        result_NOMA(realization,combination)=min(result_temp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% element-wise optimization in OMA with coupled phase shift
        beta_o=beta_r;
        theta_o=theta_r;
        for inter=1:iteration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%phase-shift optimization
            for m=1:M
                for k=1:K
                    f1(m,k)=norm((z(:,k)'*(beta_o(:,k).*theta_o(:,k))-z(m,k)'*(beta_o(m,k)*theta_o(m,k))),2)^2+norm(beta_o(m,k)*z(m,k)',2)^2;
                    f2(m,k)=z(:,k)'*(beta_o(:,k).*theta_o(:,k))-z(m,k)'*(beta_o(m,k)*theta_o(m,k));
                end
                for i=1:K
                    if(i==1)
                        cvx_begin
                        variable p(K)
                        variable theta1 complex
                        expressions gain1(K)
                        gain1(1)=f1(m,1)+2*real(z(m,1)'*(beta_o(m,1)*theta1)*f2(m,1)');
                        gain1(2)=f1(m,2)+2*real(z(m,2)'*(beta_o(m,2)*1i*theta1)*f2(m,2)');
                        minimize sum(p)
                        subject to
                        gain1(1)>=gammaO(1)*inv_pos(p(1));
                        gain1(2)>=gammaO(2)*inv_pos(p(2));
                        p>=0;
                        norm(theta1)<=1;
                        cvx_end
                        norm(theta1)
                        value(i)=sum(p);
                    else
                        cvx_begin
                        variable p(K)
                        variable T(K)
                        variable theta2 complex
                        expressions gain2(K)
                        gain2(1)=f1(m,1)+2*real(z(m,1)'*(beta_o(m,1)*theta2)*f2(m,1)');
                        gain2(2)=f1(m,2)+2*real(z(m,2)'*(beta_o(m,2)*(-1i)*theta2)*f2(m,2)');
                        minimize sum(p)
                        subject to
                        gain2(1)>=gammaO(1)*inv_pos(p(1));
                        gain2(2)>=gammaO(2)*inv_pos(p(2));
                        p>=0;
                        norm(theta2)<=1;
                        cvx_end
                        value(i)=sum(p);
                    end
                end
                if(value(1)<value(2))
                    theta_o(m,1)=theta1;
                    theta_o(m,2)=1i*theta1;
                else
                    theta_o(m,1)=theta2;
                    theta_o(m,2)=-1i*theta2;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%amplitude optimization
            for k=1:K
                w(:,k)=sqrt(10^(PL))*[r(:,k)'*diag(theta_o(1:M,k),0)*diag(g,0) d(1,k)]';
            end
            for m=1:M

                for k=1:K
                    l1(m,k)=norm((w(:,k)'*beta_o(:,k)-w(m,k)'*beta_o(m,k)),2)^2;%+norm(beta_11(m,k)*z(m,k)',2)^2;
                    l2(m,k)=w(:,k)'*(beta_o(:,k))-w(m,k)'*(beta_o(m,k));
                    l3(m,k)=2*real(w(m,k)'*l2(m,k)');
                end
                cvx_begin
                variable p(K)
                variable beta_t(K)
                expressions gain_o(K)
                for k=1:K
                    if(l3(m,k)>=0)
                        gain_o(k)=l1(m,k)+norm(w(m,k)',2)^2*beta_t(k)+l3(m,k)*sqrt(beta_t(k));
                    else
                        gain_o(k)=l1(m,k)+norm(w(m,k)',2)^2*beta_t(k)+l3(m,k)*(0.5*beta_o(m,k)^(-0.5)*(beta_t(k)-beta_o(m,k)));
                    end
                end
                minimize sum(p)
                subject to
                gain_o(1)>=gammaO(1)*inv_pos(p(1));
                gain_o(2)>=gammaO(2)*inv_pos(p(2));
                sum(beta_t)==1;
                beta_t>=0;
                cvx_end
                for k=1:K
                    beta_o(m,k)=sqrt(beta_t(k));
                end
            end
        end
        if(strcmp(cvx_status,'Infeasible')||strcmp(cvx_status,'Failed'))
            result_OMA(realization,combination)=0;
        else
            result_OMA(realization,combination)=sum(p);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end
toc;


