import numpy as np


def FUNC_randomClustering(MODEL, B, t_thres, theta, Area):
    sat = MODEL["SAT"]
    ue = MODEL["UE"]
    numSAT = sat.shape[1]
    numUE = ue.shape[1]
    distance = np.sqrt(
        (sat[0, :][:, None] - ue[0, :][None, :]) ** 2
        + (sat[1, :][:, None] - ue[1, :][None, :]) ** 2
        + (sat[2, :][:, None] - ue[2, :][None, :]) ** 2
    )
    d_thres = t_thres * 3e8
    id_clus_ue = np.zeros(numUE, dtype=int)
    numUE_in_SAT = np.zeros(numSAT, dtype=int)
    numLoop = 0

    while True:
        numLoop += 1
        u_free = np.where(id_clus_ue == 0)[0]
        if u_free.size == 0:
            break

        for u in u_free:
            if id_clus_ue[u] == 0:
                near_sat = np.where(distance[:, u] < d_thres)[0] + 1
                free_sat = np.where(numUE_in_SAT < B)[0] + 1
                ok_sat = np.intersect1d(near_sat, free_sat)
                if ok_sat.size == 0:
                    ind_sat = int(np.argmin(distance[:, u])) + 1
                    users_in_sat = np.where(id_clus_ue == ind_sat)[0]
                    if users_in_sat.size > 0:
                        id_clus_ue[users_in_sat[0]] = 0
                    id_clus_ue[u] = ind_sat
                else:
                    ind_sat = int(np.random.choice(ok_sat))
                    numUE_in_SAT[ind_sat - 1] += 1
                    id_clus_ue[u] = ind_sat

        if numLoop == 1000:
            raise RuntimeError("Cannot perform clustering due to strict constraints")

    MODEL["IdClusUE"] = id_clus_ue
    MODEL["NumUE_in_SAT"] = numUE_in_SAT
    return MODEL
