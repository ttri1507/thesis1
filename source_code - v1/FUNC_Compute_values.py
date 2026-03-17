import numpy as np


def FUNC_Compute_values(NumUE_in_SAT, IdClusUE, Rho, sigma_u, Pow):
    numSAT = len(NumUE_in_SAT)
    numUE = int(np.sum(NumUE_in_SAT))
    rate = np.zeros(numUE, dtype=float)
    sumpow = float(np.sum(Pow))
    for m in range(1, numSAT + 1):
        um = np.where(IdClusUE == m)[0]
        for u in um:
            rate[u] = np.log2(1 + Pow[m - 1, u] * Rho[m - 1, u] / sigma_u)
    sumrate = float(np.sum(rate))
    return sumrate, sumpow, rate
