% NumUE_in_SAT: 1xM
% IdClusUE: 1xU (the id cluster of each UE)
% Beta: MxU (numSAT X numUE)
% Rho: MxU (numSAT X numUE)
% Pmax_overline: 1xM
% Pow_cir_overline: 1xM
function [optimal_sumrate2, optimal_sumpow2, optimal_obj2, Pow, Rate] = FUNC_multiSAT_OPT(NumUE_in_SAT, IdClusUE, lambda_ref, Beta, Rho, sigma_u, ...
    Pmax_overline, Pow_cir_overline, tau_ref2, tau_underline2, tau_overline2)

numSAT = size(NumUE_in_SAT,2); numUE = sum(NumUE_in_SAT);

[optimal_sumrate2, optimal_sumpow2, optimal_obj2, Pow, Rate] = function_update_multiSAT(NumUE_in_SAT, IdClusUE, lambda_ref,...
    Beta, Rho, sigma_u, Pmax_overline, tau_ref2, Pow_cir_overline, numSAT, numUE);

while (optimal_obj2 > 0)
    tau_underline2 = tau_ref2;
    tau_ref2 = 2*tau_ref2;
    [optimal_sumrate2, optimal_sumpow2, optimal_obj2, Pow, Rate] = function_update_multiSAT(NumUE_in_SAT, IdClusUE, lambda_ref,...
    Beta, Rho, sigma_u, Pmax_overline, tau_ref2, Pow_cir_overline, numSAT, numUE);
end

% %%%%%%%%%%%% Loop for tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter_tau = 0; epsilon_tau = 1;
while (epsilon_tau >= 1e-9)
    iter_tau = iter_tau + 1;
    
    % %%%%%%%%%%%% Loop for dual value %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [optimal_sumrate2, optimal_sumpow2, optimal_obj2, Pow, Rate] = function_update_multiSAT(NumUE_in_SAT, IdClusUE, lambda_ref,...
    Beta, Rho, sigma_u, Pmax_overline, tau_ref2, Pow_cir_overline, numSAT, numUE);
    
    if (optimal_obj2 > 0)
        tau_underline2 = tau_ref2;
    else
        tau_overline2 = tau_ref2;
    end
    tau_ref2 = (tau_underline2 + tau_overline2)/2;
    epsilon_tau = abs(tau_overline2 - tau_underline2);
%     fprintf("The optimal value = %f \n", optimal_obj2);
%     fprintf("The optimal tau = %f \n", tau_ref2);
    
    loop_optimal_obj2(iter_tau) = optimal_obj2;
    loop_tau_ref(iter_tau) = tau_ref2;
    
    if (iter_tau >= 200)
        epsilon_tau = 0;
    end
end

end

function [optimal_sumrate2, optimal_sumpow2, optimal_obj2, Pow, Rate] = function_update_multiSAT(NumUE_in_SAT, IdClusUE, lambda_ref,...
    Beta, Rho, sigma_u, Pmax_overline, tau_ref2, Pow_cir_overline, numSAT, numUE)
optimal_sumrate2 = 0; optimal_sumpow2 = 0; optimal_obj2 = 0;
Pow = zeros(numSAT, numUE); Rate = zeros(1, numUE);
for m = 1:numSAT
    % users in SAT m
    Um = find(IdClusUE==m);
    Beta_m = Beta(m,Um); Rho_m = Rho(m,Um);
    [optimal_sumrate_m, optimal_sumpow_m, optimal_obj_m, Pow_Um_m, Rate_Um_m] = FUNC_oneSAT_dualvalue_loop(NumUE_in_SAT(m), lambda_ref,...
    Beta_m, Rho_m, sigma_u, Pmax_overline(3), tau_ref2, Pow_cir_overline(m));
    optimal_sumrate2 = optimal_sumrate2 + optimal_sumrate_m;
    optimal_sumpow2 = optimal_sumpow2 + optimal_sumpow_m;
    optimal_obj2 = optimal_obj2 + optimal_obj_m;
    Pow(m,Um) = Pow_Um_m; Rate(Um) = Rate_Um_m;
end

end

