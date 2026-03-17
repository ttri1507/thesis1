import numpy as np


def FUNC_PATH_SAT(MODEL, N, Bw, f):
    sat = MODEL["SAT"]
    ue = MODEL["UE"]
    numSAT = sat.shape[1]
    numUE = ue.shape[1]
    distance = np.sqrt(
        (sat[0, :][:, None] - ue[0, :][None, :]) ** 2
        + (sat[1, :][:, None] - ue[1, :][None, :]) ** 2
        + (sat[2, :][:, None] - ue[2, :][None, :]) ** 2
    )
    h = np.zeros((numSAT, numUE, N), dtype=complex)
    for m in range(numSAT):
        for u in range(numUE):
            h[m, u, :] = func_path_1sat_1ss(N, distance[m, u], f, Bw)
    return h


def func_path_1sat_1ss(N, d_mu, f, Bw):
    lambda_ = 3e8 / f
    kB = 1.38e-23
    TR = 235.3
    GR = 10 ** (40.7 / 10)
    A_mu = 1
    phi = 2 * np.pi * np.random.rand(N)
    g_mu = np.exp(1j * phi)
    return (1 / np.sqrt(A_mu)) * (np.sqrt(GR) * g_mu) / (
        4 * np.pi * d_mu / lambda_ * np.sqrt(kB * TR * Bw)
    )
