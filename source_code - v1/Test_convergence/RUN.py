import time

import matplotlib.pyplot as plt
import numpy as np

from FUNC_CoalitionGame import FUNC_CoalitionGame
from FUNC_Create_SAT_UE import FUNC_Create_SAT_UE
from FUNC_PATH_SAT import FUNC_PATH_SAT
from FUNC_multiSAT_OPT import FUNC_multiSAT_OPT
from FUNC_randomClustering import FUNC_randomClustering


def run():
    N = 100
    B = 25
    numSAT = 4
    numUE = 80
    Pmax = 10 ** (20 / 10)
    P0 = 10 ** (10 / 10)
    Bw = 500 * 1e6
    Noise_var = 1
    f = 20e9
    Area = np.array([2000e3, 2000e3], dtype=float)
    zm = 1000e3
    t_thres = 2000e3 / 3e8
    Rmin = 0.2

    model = {}
    sat, ue = FUNC_Create_SAT_UE(Area, numSAT, numUE, zm)
    model["SAT"] = sat
    model["UE"] = ue
    model["Channel"] = FUNC_PATH_SAT(model, N, Bw, f)

    theta = 0.5
    model_rc = FUNC_randomClustering(model.copy(), B, t_thres, theta, Area)
    id_clus_ue_rc = model_rc["IdClusUE"]
    numUE_in_SAT_rc = model_rc["NumUE_in_SAT"]
    start_gtc = time.time()
    model_gtc = FUNC_CoalitionGame(model.copy(), B, t_thres, theta, Area, id_clus_ue_rc, numUE_in_SAT_rc)
    _ = time.time() - start_gtc

    h_norm2 = np.sqrt(np.abs(np.sum(np.conj(model["Channel"]) * model["Channel"], axis=2)))
    precoding = np.conj(model["Channel"]) / np.repeat(h_norm2[:, :, np.newaxis], N, axis=2)
    rho = np.abs(np.sum(model["Channel"] * precoding, axis=2)) ** 2
    sigma_u = Noise_var
    pow_overline = sigma_u * (2 ** Rmin - 1) / rho
    beta = 1 + pow_overline * rho / sigma_u

    numUE_in_SAT = model_gtc["NumUE_in_SAT"]
    id_clus_ue = model_gtc["IdClusUE"]
    lambda_ref = 1
    tau_ref2 = 1
    tau_underline2 = 1e-10
    tau_overline2 = 1000
    sum_p_overline = np.zeros(numSAT, dtype=float)
    for m in range(1, numSAT + 1):
        um = np.where(id_clus_ue == m)[0]
        sum_p_overline[m - 1] = np.sum(pow_overline[m - 1, um])
    pmax_overline = Pmax - sum_p_overline
    pow_cir_overline = P0 + sum_p_overline
    sumrate_gtc, sumpow_gtc, _, _, _, loop_ee_gtc = FUNC_multiSAT_OPT(
        numUE_in_SAT,
        id_clus_ue,
        lambda_ref,
        beta,
        rho,
        sigma_u,
        pmax_overline,
        pow_cir_overline,
        tau_ref2,
        tau_underline2,
        tau_overline2,
    )
    _ = sumrate_gtc / (np.sum(pow_cir_overline) + sumpow_gtc)

    numUE_in_SAT = model_rc["NumUE_in_SAT"]
    id_clus_ue = model_rc["IdClusUE"]
    sum_p_overline = np.zeros(numSAT, dtype=float)
    for m in range(1, numSAT + 1):
        um = np.where(id_clus_ue == m)[0]
        sum_p_overline[m - 1] = np.sum(pow_overline[m - 1, um])
    pmax_overline = Pmax - sum_p_overline
    pow_cir_overline = P0 + sum_p_overline
    sumrate_rc, sumpow_rc, _, _, _, loop_ee_rc = FUNC_multiSAT_OPT(
        numUE_in_SAT,
        id_clus_ue,
        lambda_ref,
        beta,
        rho,
        sigma_u,
        pmax_overline,
        pow_cir_overline,
        tau_ref2,
        tau_underline2,
        tau_overline2,
    )
    _ = sumrate_rc / (np.sum(pow_cir_overline) + sumpow_rc)

    x = np.arange(1, 21)
    n = len(x)
    plt.semilogy(x, loop_ee_gtc[:n], "r-*", markersize=4, linewidth=1)
    plt.grid(True)
    plt.semilogy(x, loop_ee_rc[:n], "b--^", markersize=4, linewidth=1)
    plt.legend(["EE-GTC", "EE-RC"])
    plt.xlabel("Number of iterations")
    plt.ylabel("EE (bits/Joule/Hz)")
    plt.show()


if __name__ == "__main__":
    run()
