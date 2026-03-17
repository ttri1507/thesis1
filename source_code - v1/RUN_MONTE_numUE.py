import numpy as np
import matplotlib.pyplot as plt

from FUNC_CoalitionGame import FUNC_CoalitionGame
from FUNC_Compute_values import FUNC_Compute_values
from FUNC_Create_SAT_UE import FUNC_Create_SAT_UE
from FUNC_PATH_SAT import FUNC_PATH_SAT
from FUNC_multiSAT_OPT import FUNC_multiSAT_OPT
from FUNC_randomClustering import FUNC_randomClustering


def run():
    N = 100
    B = 25
    numSAT = 4
    Pmax = 10 ** (20 / 10)
    P0 = 10 ** (10 / 10)
    Bw = 500 * 1e6
    Noise_var = 1
    f = 20e9
    Area = np.array([2000e3, 2000e3], dtype=float)
    zm = 1000e3
    t_thres = 2000e3 / 3e8
    Rmin = 0.2

    numMonte = 100
    NumUE = np.array([50, 60, 70, 80, 90, 100], dtype=int)
    ee_total = np.zeros((4, len(NumUE)), dtype=float)

    for numScen in range(len(NumUE)):
        print(f"SCENARIO NUMBER {numScen + 1}")
        ee_monte = np.zeros((4, numMonte), dtype=float)
        numCurrent = 0
        while numCurrent < numMonte:
            model = {}
            sat, ue = FUNC_Create_SAT_UE(Area, numSAT, int(NumUE[numScen]), zm)
            model["SAT"] = sat
            model["UE"] = ue
            model["Channel"] = FUNC_PATH_SAT(model, N, Bw, f)

            theta = 0.5
            model_rc = FUNC_randomClustering(model.copy(), B, t_thres, theta, Area)
            id_clus_ue_rc = model_rc["IdClusUE"]
            numUE_in_SAT_rc = model_rc["NumUE_in_SAT"]
            model_gtc = FUNC_CoalitionGame(model.copy(), B, t_thres, theta, Area, id_clus_ue_rc, numUE_in_SAT_rc)

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
            sumrate_gtcpa, sumpow_gtcpa, _, _, _ = FUNC_multiSAT_OPT(
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
            ee_gtcpa = sumrate_gtcpa / (np.sum(pow_cir_overline) + sumpow_gtcpa)

            numUE_in_SAT = model_rc["NumUE_in_SAT"]
            id_clus_ue = model_rc["IdClusUE"]
            sum_p_overline = np.zeros(numSAT, dtype=float)
            for m in range(1, numSAT + 1):
                um = np.where(id_clus_ue == m)[0]
                sum_p_overline[m - 1] = np.sum(pow_overline[m - 1, um])
            pmax_overline = Pmax - sum_p_overline
            pow_cir_overline = P0 + sum_p_overline
            sumrate_rcpa, sumpow_rcpa, _, _, _ = FUNC_multiSAT_OPT(
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
            ee_rcpa = sumrate_rcpa / (np.sum(pow_cir_overline) + sumpow_rcpa)

            numUE_in_SAT = model_gtc["NumUE_in_SAT"]
            id_clus_ue = model_gtc["IdClusUE"]
            p_each_ue = Pmax / numUE_in_SAT
            pow_gtcep = np.zeros((numSAT, int(NumUE[numScen])), dtype=float)
            for u in range(int(NumUE[numScen])):
                m = id_clus_ue[u]
                pow_gtcep[m - 1, u] = p_each_ue[m - 1]
            sumrate_gtcep, sumpow_gtcep, _ = FUNC_Compute_values(numUE_in_SAT, id_clus_ue, rho, sigma_u, pow_gtcep)
            ee_gtcep = sumrate_gtcep / (sumpow_gtcep + numSAT * P0)

            numUE_in_SAT = model_rc["NumUE_in_SAT"]
            id_clus_ue = model_rc["IdClusUE"]
            p_each_ue = Pmax / numUE_in_SAT
            pow_rcep = np.zeros((numSAT, int(NumUE[numScen])), dtype=float)
            for u in range(int(NumUE[numScen])):
                m = id_clus_ue[u]
                pow_rcep[m - 1, u] = p_each_ue[m - 1]
            sumrate_rcep, sumpow_rcep, _ = FUNC_Compute_values(numUE_in_SAT, id_clus_ue, rho, sigma_u, pow_rcep)
            ee_rcep = sumrate_rcep / (sumpow_rcep + numSAT * P0)

            ee_monte[0, numCurrent] = ee_gtcpa
            ee_monte[1, numCurrent] = ee_rcpa
            ee_monte[2, numCurrent] = ee_gtcep
            ee_monte[3, numCurrent] = ee_rcep
            numCurrent += 1
        ee_total[:, numScen] = np.sum(ee_monte, axis=1) / numMonte

    x = NumUE
    plt.semilogy(x, ee_total[0, :], "r--^", markersize=4, linewidth=1)
    plt.grid(True)
    plt.semilogy(x, ee_total[1, :], "b--^", markersize=4, linewidth=1)
    plt.semilogy(x, ee_total[2, :], "r-*", markersize=4, linewidth=1)
    plt.semilogy(x, ee_total[3, :], "b-*", markersize=4, linewidth=1)
    plt.legend(["EE_GTCPA", "EE_RCPA", "EE_GTCEP", "EE_RCEP"])
    plt.xlabel("Number of UEs")
    plt.ylabel("EE (bits/Joule/Hz)")
    plt.show()


if __name__ == "__main__":
    run()
