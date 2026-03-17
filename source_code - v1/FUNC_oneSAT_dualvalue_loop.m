function [optimal_sumrate2, optimal_sumpow2, optimal_obj2, Pow_Um, Rate_Um] = FUNC_oneSAT_dualvalue_loop(numUE, lambda_ref,...
    Beta_m, Rho_m, sigma_u, Pmax_overline, tau_ref2, Pow_cir_overline)

% %%%%%%%%%%%% Loop for dual value %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pow_Um = zeros(1,numUE); sum_pow_wf = 0;
for iUser = 1:numUE
    term_pow_wf = 1/(log(2)*tau_ref2) - Beta_m(iUser)*sigma_u/Rho_m(iUser);
    Pow_Um(iUser) = max(0,term_pow_wf);
    sum_pow_wf = sum_pow_wf + Pow_Um(iUser);
end

% %%%%%%% Case dual value = 0 %%%%%%%%%%%%%
if (sum_pow_wf <= Pmax_overline)
    sum_rate_wf = 0; Rate_Um = zeros(1,numUE);
    for iUser = 1:numUE
        Rate_Um(iUser) = log2(Beta_m(iUser) + Pow_Um(iUser)*Rho_m(iUser)/sigma_u);
        sum_rate_wf = sum_rate_wf + log2(Beta_m(iUser) + Pow_Um(iUser)*Rho_m(iUser)/sigma_u);
    end
    optimal_sumrate2 = sum_rate_wf;
    optimal_sumpow2 = sum_pow_wf;
    optimal_obj2 = optimal_sumrate2 - tau_ref2*(optimal_sumpow2 + Pow_cir_overline);
    
    % %%%%%%% Case dual value > 0 %%%%%%%%%%%%%
else
    lambda_sol = lambda_ref; epsilon_lambda = 1;
    lambda_underline = 0; lambda_overline = 1;
    iter_lambda = 0;
    
    % %%%%%%%%% Find lower bound of dual value %%%%%%%%%%%%%%%%%%%
    Pow_Um = zeros(1,numUE); sum_pow_wf = 0;
    for iUser = 1:numUE
        term_pow_wf = 1/(log(2)*(tau_ref2+lambda_sol)) - Beta_m(iUser)*sigma_u/Rho_m(iUser);
        Pow_Um(iUser) = max(0,term_pow_wf);
        sum_pow_wf = sum_pow_wf + Pow_Um(iUser);
    end
    while (sum_pow_wf > Pmax_overline)
        lambda_underline = lambda_sol;
        lambda_sol = 2*lambda_sol;
        Pow_Um = zeros(1,numUE); sum_pow_wf = 0;
        for iUser = 1:numUE
            term_pow_wf = 1/(log(2)*(tau_ref2+lambda_sol)) - Beta_m(iUser)*sigma_u/Rho_m(iUser);
            Pow_Um(iUser) = max(0,term_pow_wf);
            sum_pow_wf = sum_pow_wf + Pow_Um(iUser);
        end
    end
    
    % %%%%%%%%%%%% loop %%%%%%%%%%%%%%%%%%%%%
    while (epsilon_lambda >= 1e-9)
        iter_lambda = iter_lambda + 1;
        Pow_Um = zeros(1,numUE); sum_pow_wf = 0;
        for iUser = 1:numUE
            term_pow_wf = 1/(log(2)*(tau_ref2+lambda_sol)) - Beta_m(iUser)*sigma_u/Rho_m(iUser);
            Pow_Um(iUser) = max(0,term_pow_wf);
            sum_pow_wf = sum_pow_wf + Pow_Um(iUser);
        end
        
        if (sum_pow_wf > Pmax_overline)
            lambda_underline = lambda_sol;
        else
            lambda_overline = lambda_sol;
        end
        
        lambda_sol = (lambda_underline + lambda_overline)/2;
        epsilon_lambda = abs(lambda_overline - lambda_underline);
        loop_lambda_ref(iter_lambda) = lambda_sol;
        
        if(iter_lambda >= 500)
            epsilon_lambda = 0;
        end
    end
    
    sum_rate_wf = 0; Rate_Um = zeros(1,numUE);
    for iUser = 1:numUE
        Rate_Um(iUser) = log2(Beta_m(iUser) + Pow_Um(iUser)*Rho_m(iUser)/sigma_u);
        sum_rate_wf = sum_rate_wf + log2(Beta_m(iUser) + Pow_Um(iUser)*Rho_m(iUser)/sigma_u);
    end
    optimal_sumrate2 = sum_rate_wf;
    optimal_sumpow2 = sum_pow_wf;
    optimal_obj2 = optimal_sumrate2 - tau_ref2*(optimal_sumpow2 + Pow_cir_overline);
end

end

