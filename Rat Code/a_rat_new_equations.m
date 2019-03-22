
% Sigma_tgf
f(10) = Sigma_tgf - ( 0.3408 + 3.449 / (3.88 + exp((Phi_mdsod - 3.859 * SF) / (-0.9617 * SF) )) );

% gamma_filsod
f(14) = gamma_filsod - ( 0.85 + 0.3 / (1 + exp((Phi_filsod - 18 * SF)/(138 * SF) )) );

% psi_al
f(20) = psi_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );

% lambda_dt
f(24) = lambda_dt - ( 0.82 + 0.39 / (1 + exp((Phi_dtsod - 1.7625 * SF) / (0.375 * SF) )) );

% Phi_win
f(27) = Phi_win - ( max( 0, 0.008 * SF / (1 + 86.1*1.822*(C_adh^(3*-1.607))) - 0.005 * SF ) );

% V_b
f(29) = V_b - ( 4.5479392962 * SF + 2.431217 * SF / (1 + exp(-(V_ecf - 18.11278 * SF) * (0.47437 / SF) )) );
% P_mf
f(30) = P_mf - ( ( (7.436 / SF) * V_b - 30.18) * epsilon_aum );

% P_ra
if     strcmp(gender,'male')
    f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8268 );
%     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8268 ) );
elseif strcmp(gender,'female')
    f(33) = P_ra - ( 0.2787 * exp(Phi_co * 0.2281 / SF) - 0.8245 );
%     f(33) = P_ra - ( max( 0, 0.2787 * exp(Phi_co * 0.2281) - 0.8245 ) );
end

% vas_f
f(35) = vas_f - ( (11.312 * exp(-Phi_co * 0.4799 / SF)) / 100000 );

% Phi_twreab
f(50) = Phi_twreab - ( 0.025 * SF - 0.001 * SF / (mu_al * mu_adh) + 0.8 * Phi_gfilt );
% mu_al
f(51) = mu_al - ( 0.17 + 0.94 / (1 + exp((0.48 - 0.87 * log10(C_al)) / 0.88)) );

% Phi_u
f(53) = Phi_u - ( max( 0.0003 * SF, Phi_gfilt - Phi_twreab ) );

% nu_mdsod
if     strcmp(gender,'male')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.731 * SF) / (0.6056 * SF) )) );
elseif strcmp(gender,'female')
    f(56) = nu_mdsod - ( 0.2262 + 28.04 / (11.56 + exp((Phi_mdsod - 1.637 * SF) / (0.6056 * SF) )) );
end

% C_al
if     strcmp(gender,'male')
    f(58) = C_al - ( max( 1, N_al * 395.3 ) );
elseif strcmp(gender,'female')
    f(58) = C_al - ( max( 1, N_al * 379.4 ) );
end


































