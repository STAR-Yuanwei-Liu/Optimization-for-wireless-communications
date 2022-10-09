function [W, theta_t, theta_r, theta_tt, theta_rr, sum_rate_all, phas_diff_all] = algorithm_BCD(para, W, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho, G, h)
%The BCD algorithm for solve the augmented Lagrangian problem
%Date: 01/10/2022
%Author: Zhaolin Wang

sum_rate_pre = 0;
phas_diff_all = [];
sum_rate_all = [];

for i = 1:40
    % update auxiliary variables
    [omega, upsilon] = update_weights(para, W, theta_t, theta_r, G, h);
    % update beamformers
    [W] = update_W(para, omega, upsilon, theta_t, theta_r, G, h);
    % update original STAR coefficients
    [theta_t, theta_r] = update_theta(para, omega, upsilon, W, theta_tt, theta_rr, lambda_t, lambda_r, rho, G, h);
    % update auxiliary STAR coefficients
    [theta_tt, theta_rr] = update_theta_aux(para, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho);
    
    [gamma, ~] = SINR(para, W, theta_t, theta_r, G, h);
    sum_rate = sum(log2(1 + gamma));
    disp(['Inner loop - ' num2str(i) ', Sum rate - ' num2str(sum_rate)]);
    sum_rate_all = [sum_rate_all, sum_rate];
    phas_diff_all = [phas_diff_all abs(angle(theta_r) - angle(theta_t))];

    % check convergence
    reduction = abs(sum_rate - sum_rate_pre) / sum_rate;
    if reduction < 1e-3 % convergence
        break; 
    end
    sum_rate_pre = sum_rate;

end


end

%% update auxiliary variables
function [omega, upsilon] = update_weights(para, W, theta_t, theta_r, G, h)
    [gamma, upsilon] = SINR(para, W, theta_t, theta_r, G, h);
    omega = (1 + gamma)/log(2);
end

%% update beamformers
function [W] = update_W(para, omega, upsilon, theta_t, theta_r, G, h)
    Theta_t = diag(theta_t); ht = G' * Theta_t' * h; 
    Theta_r = diag(theta_r); hr = G' * Theta_r' * h; 
    cvx_begin quiet
        % optimization variables
        variable W(para.M, para.K) complex

        % constraint
        square_pos(norm(W,'fro')) <= para.Pt;
        
        % calculate objective function
        obj = 0;
        for k = 1:para.K
            if k <= para.K/2
                hk = ht(:,k); 
            else
                hk = hr(:,k);
            end   
            wk = W(:,k);
            ek = abs(upsilon(k))^2 * (square_pos(norm(sqrtm(hk*hk')*W,'fro')) + 1) - 2*real( conj(upsilon(k))*hk'*wk ) + 1;
            obj = obj + omega(k)*ek;
        end


        minimize(obj);
    cvx_end   
end

%% update original STAR coefficients
function [theta_t, theta_r] = update_theta(para, omega, upsilon, W, theta_tt, theta_rr, lambda_t, lambda_r, rho, G, h)

    cvx_begin quiet
        % optimization variables
        variable theta_t(para.N, 1) complex
        variable theta_r(para.N, 1) complex

        % constraint
        for i = 1:para.N 
            theta_t(i)*conj(theta_t(i)) + theta_r(i)*conj(theta_r(i)) <= 1;
        end
        
        % calculate objective function
        Theta_t = diag(theta_t); ht = G' * Theta_t' * h; 
        Theta_r = diag(theta_r); hr = G' * Theta_r' * h; 
        obj = 0;
        for k = 1:para.K
            if k <= para.K/2
                hk = ht(:,k); 
            else
                hk = hr(:,k);
            end   
            wk = W(:,k);
            ek = abs(upsilon(k))^2 * (quad_form(hk, W*W') + 1) - 2*real( conj(upsilon(k))*hk'*wk ) + 1;
            obj = obj + omega(k)*ek;
        end
        penalty = sum_square_abs(theta_tt - theta_t + rho*lambda_t) + sum_square_abs(theta_rr - theta_r + rho*lambda_r);
        obj = obj + 1/(2*rho) * penalty;
        minimize(obj);
    cvx_end  


end

%% update auxiliary STAR coefficients
function [theta_tt, theta_rr] = update_theta_aux(para, theta_t, theta_r, theta_tt, theta_rr, lambda_t, lambda_r, rho)

beta_tt = abs(theta_tt); 
beta_rr = abs(theta_rr); 

% update phase shift
[qtt, qrr] = update_phase(para, theta_t, theta_r, beta_tt, beta_rr, lambda_t, lambda_r, rho);

% update amplitude
[beta_tt, beta_rr] = update_amplitude(para, theta_t, theta_r, qtt, qrr, lambda_t, lambda_r, rho);

% combine amplitudes ad phase shifts
theta_tt = diag(beta_tt)*qtt;
theta_rr = diag(beta_rr)*qrr;

end

%% update phase shifts
function [qtt, qrr] = update_phase(para, theta_t, theta_r, beta_tt, beta_rr, lambda_t, lambda_r, rho)

qtt = zeros(para.N,1); qrr = zeros(para.N,1);

theta_t = -theta_t + rho*lambda_t; theta_t = theta_t'*diag(beta_tt);
theta_r = -theta_r + rho*lambda_r; theta_r = theta_r'*diag(beta_rr);

for n = 1:para.N
    
    phi_p = theta_t(n) + 1i*theta_r(n);
    phi_m = theta_t(n) - 1i*theta_r(n);

    qtt_1 = exp(1i*(pi - angle(phi_p))); qrr_1 =  1i*qtt_1;
    qtt_2 = exp(1i*(pi - angle(phi_m))); qrr_2 =  -1i*qtt_2;


    o1 = real(theta_t(n)*qtt_1) + real(theta_r(n)*qrr_1);
    o2 = real(theta_t(n)*qtt_2) + real(theta_r(n)*qrr_2);

    if o1 < o2
        qtt(n) = qtt_1; 
        qrr(n) = qrr_1;
    else
        qtt(n) = qtt_2; 
        qrr(n) = qrr_2;
    end
end

end

%% update amplitudes
function [beta_tt, beta_rr] = update_amplitude(para, theta_t, theta_r, qtt, qrr, lambda_t, lambda_r, rho)

beta_tt = zeros(para.N,1); beta_rr = zeros(para.N,1);
theta_t = -theta_t + rho*lambda_t; theta_t = theta_t'*diag(qtt);
theta_r = -theta_r + rho*lambda_r; theta_r = theta_r'*diag(qrr);

for n = 1:para.N
    a = abs(theta_t(n))*cos(angle(theta_t(n)));
    b = abs(theta_r(n))*cos(angle(theta_r(n)));
    
    phi = sign(b) * acos(a / sqrt(a^2 + b^2));

    if phi >= -pi && phi <-1/2*pi
        w = -1/2*pi - phi;
    elseif phi >= -1/2*pi && phi <= 1/4*pi
        w = 0;
    else
        w = 1/2*pi;
    end
    beta_tt(n) = sin(w);
    beta_rr(n) = cos(w);

end

end
