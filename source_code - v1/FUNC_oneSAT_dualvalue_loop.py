import numpy as np


def FUNC_oneSAT_dualvalue_loop(numUE, lambda_ref, Beta_m, Rho_m, sigma_u, Pmax_overline, tau_ref2, Pow_cir_overline):
    if numUE == 0:
        return 0.0, 0.0, -tau_ref2 * Pow_cir_overline, np.zeros(0), np.zeros(0)

    pow_um = np.zeros(numUE, dtype=float)
    sum_pow_wf = 0.0
    for i_user in range(numUE):
        term_pow_wf = 1 / (np.log(2) * tau_ref2) - Beta_m[i_user] * sigma_u / Rho_m[i_user]
        pow_um[i_user] = max(0.0, term_pow_wf)
        sum_pow_wf += pow_um[i_user]

    if sum_pow_wf <= Pmax_overline:
        rate_um = np.zeros(numUE, dtype=float)
        sum_rate_wf = 0.0
        for i_user in range(numUE):
            rate_um[i_user] = np.log2(Beta_m[i_user] + pow_um[i_user] * Rho_m[i_user] / sigma_u)
            sum_rate_wf += rate_um[i_user]
        optimal_sumrate2 = sum_rate_wf
        optimal_sumpow2 = sum_pow_wf
        optimal_obj2 = optimal_sumrate2 - tau_ref2 * (optimal_sumpow2 + Pow_cir_overline)
    else:
        lambda_sol = lambda_ref
        epsilon_lambda = 1.0
        lambda_underline = 0.0
        lambda_overline = 1.0
        iter_lambda = 0

        pow_um = np.zeros(numUE, dtype=float)
        sum_pow_wf = 0.0
        for i_user in range(numUE):
            term_pow_wf = 1 / (np.log(2) * (tau_ref2 + lambda_sol)) - Beta_m[i_user] * sigma_u / Rho_m[i_user]
            pow_um[i_user] = max(0.0, term_pow_wf)
            sum_pow_wf += pow_um[i_user]
        while sum_pow_wf > Pmax_overline:
            lambda_underline = lambda_sol
            lambda_sol = 2 * lambda_sol
            pow_um = np.zeros(numUE, dtype=float)
            sum_pow_wf = 0.0
            for i_user in range(numUE):
                term_pow_wf = 1 / (np.log(2) * (tau_ref2 + lambda_sol)) - Beta_m[i_user] * sigma_u / Rho_m[i_user]
                pow_um[i_user] = max(0.0, term_pow_wf)
                sum_pow_wf += pow_um[i_user]

        while epsilon_lambda >= 1e-9:
            iter_lambda += 1
            pow_um = np.zeros(numUE, dtype=float)
            sum_pow_wf = 0.0
            for i_user in range(numUE):
                term_pow_wf = 1 / (np.log(2) * (tau_ref2 + lambda_sol)) - Beta_m[i_user] * sigma_u / Rho_m[i_user]
                pow_um[i_user] = max(0.0, term_pow_wf)
                sum_pow_wf += pow_um[i_user]
            if sum_pow_wf > Pmax_overline:
                lambda_underline = lambda_sol
            else:
                lambda_overline = lambda_sol
            lambda_sol = (lambda_underline + lambda_overline) / 2
            epsilon_lambda = abs(lambda_overline - lambda_underline)
            if iter_lambda >= 500:
                epsilon_lambda = 0.0

        rate_um = np.zeros(numUE, dtype=float)
        sum_rate_wf = 0.0
        for i_user in range(numUE):
            rate_um[i_user] = np.log2(Beta_m[i_user] + pow_um[i_user] * Rho_m[i_user] / sigma_u)
            sum_rate_wf += rate_um[i_user]
        optimal_sumrate2 = sum_rate_wf
        optimal_sumpow2 = sum_pow_wf
        optimal_obj2 = optimal_sumrate2 - tau_ref2 * (optimal_sumpow2 + Pow_cir_overline)

    return optimal_sumrate2, optimal_sumpow2, optimal_obj2, pow_um, rate_um
