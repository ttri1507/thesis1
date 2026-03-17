import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parent.parent))

from FUNC_oneSAT_dualvalue_loop import FUNC_oneSAT_dualvalue_loop


def FUNC_multiSAT_OPT(
    NumUE_in_SAT,
    IdClusUE,
    lambda_ref,
    Beta,
    Rho,
    sigma_u,
    Pmax_overline,
    Pow_cir_overline,
    tau_ref2,
    tau_underline2,
    tau_overline2,
):
    numSAT = len(NumUE_in_SAT)
    numUE = int(np.sum(NumUE_in_SAT))

    optimal_sumrate2, optimal_sumpow2, optimal_obj2, pow_, rate = function_update_multiSAT(
        NumUE_in_SAT,
        IdClusUE,
        lambda_ref,
        Beta,
        Rho,
        sigma_u,
        Pmax_overline,
        tau_ref2,
        Pow_cir_overline,
        numSAT,
        numUE,
    )

    while optimal_obj2 > 0:
        tau_underline2 = tau_ref2
        tau_ref2 = 2 * tau_ref2
        optimal_sumrate2, optimal_sumpow2, optimal_obj2, pow_, rate = function_update_multiSAT(
            NumUE_in_SAT,
            IdClusUE,
            lambda_ref,
            Beta,
            Rho,
            sigma_u,
            Pmax_overline,
            tau_ref2,
            Pow_cir_overline,
            numSAT,
            numUE,
        )

    iter_tau = 0
    epsilon_tau = 1.0
    loop_ee = []
    while epsilon_tau >= 1e-9:
        iter_tau += 1
        optimal_sumrate2, optimal_sumpow2, optimal_obj2, pow_, rate = function_update_multiSAT(
            NumUE_in_SAT,
            IdClusUE,
            lambda_ref,
            Beta,
            Rho,
            sigma_u,
            Pmax_overline,
            tau_ref2,
            Pow_cir_overline,
            numSAT,
            numUE,
        )
        if optimal_obj2 > 0:
            tau_underline2 = tau_ref2
        else:
            tau_overline2 = tau_ref2
        tau_ref2 = (tau_underline2 + tau_overline2) / 2
        epsilon_tau = abs(tau_overline2 - tau_underline2)
        loop_ee.append(optimal_sumrate2 / (np.sum(Pow_cir_overline) + optimal_sumpow2))
        if iter_tau >= 200:
            epsilon_tau = 0.0

    return optimal_sumrate2, optimal_sumpow2, optimal_obj2, pow_, rate, np.array(loop_ee, dtype=float)


def function_update_multiSAT(
    NumUE_in_SAT,
    IdClusUE,
    lambda_ref,
    Beta,
    Rho,
    sigma_u,
    Pmax_overline,
    tau_ref2,
    Pow_cir_overline,
    numSAT,
    numUE,
):
    optimal_sumrate2 = 0.0
    optimal_sumpow2 = 0.0
    optimal_obj2 = 0.0
    pow_ = np.zeros((numSAT, numUE), dtype=float)
    rate = np.zeros(numUE, dtype=float)
    pmax_fixed = Pmax_overline[2] if len(Pmax_overline) >= 3 else Pmax_overline[-1]

    for m in range(1, numSAT + 1):
        um = np.where(IdClusUE == m)[0]
        beta_m = Beta[m - 1, um]
        rho_m = Rho[m - 1, um]
        optimal_sumrate_m, optimal_sumpow_m, optimal_obj_m, pow_um_m, rate_um_m = FUNC_oneSAT_dualvalue_loop(
            int(NumUE_in_SAT[m - 1]),
            lambda_ref,
            beta_m,
            rho_m,
            sigma_u,
            pmax_fixed,
            tau_ref2,
            Pow_cir_overline[m - 1],
        )
        optimal_sumrate2 += optimal_sumrate_m
        optimal_sumpow2 += optimal_sumpow_m
        optimal_obj2 += optimal_obj_m
        pow_[m - 1, um] = pow_um_m
        rate[um] = rate_um_m
    return optimal_sumrate2, optimal_sumpow2, optimal_obj2, pow_, rate
