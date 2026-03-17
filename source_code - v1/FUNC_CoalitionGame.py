import numpy as np


def FUNC_CoalitionGame(MODEL, B, t_thres, theta, Area, IdClusUE_RC, NumUE_in_SAT_RC):
    sat = MODEL["SAT"]
    ue = MODEL["UE"]
    chan = np.sqrt(np.abs(np.sum(np.conj(MODEL["Channel"]) * MODEL["Channel"], axis=2)))
    h_total = np.sum(chan)
    distance = np.sqrt(
        (sat[0, :][:, None] - ue[0, :][None, :]) ** 2
        + (sat[1, :][:, None] - ue[1, :][None, :]) ** 2
        + (sat[2, :][:, None] - ue[2, :][None, :]) ** 2
    )
    t_total = np.sum(distance) / 3e8
    d_thres = t_thres * 3e8
    c_not = distance > d_thres
    id_old = np.ones(ue.shape[1], dtype=int)
    id_clus_ue = np.array(IdClusUE_RC, dtype=int).copy()
    numUE_in_SAT = np.array(NumUE_in_SAT_RC, dtype=int).copy()

    numChange_old = ue.shape[1]
    numTimechange = 0
    numLoop = 0
    while not np.all(id_old == id_clus_ue):
        numLoop += 1
        id_old = id_clus_ue.copy()
        for u in range(ue.shape[1]):
            mp = id_clus_ue[u]
            vmp_old = func_utility(distance, chan, mp, id_clus_ue, theta, t_total, h_total)
            id_ok = id_clus_ue.copy()
            num_ok = numUE_in_SAT.copy()
            gap = 0.0
            m_ok_u = np.where(c_not[:, u] == 0)[0] + 1
            m_ok_u = m_ok_u[m_ok_u != mp]
            for m in m_ok_u:
                id_temp = id_clus_ue.copy()
                num_temp = numUE_in_SAT.copy()
                nUEm = num_temp[m - 1]
                vm_old = func_utility(distance, chan, m, id_clus_ue, theta, t_total, h_total)
                if nUEm == B:
                    um = np.where(id_clus_ue == m)[0]
                    for up in um:
                        if c_not[mp - 1, up]:
                            continue
                        id_temp = id_clus_ue.copy()
                        id_temp[u] = m
                        id_temp[up] = mp
                        vmp_new = func_utility(distance, chan, mp, id_temp, theta, t_total, h_total)
                        vm_new = func_utility(distance, chan, m, id_temp, theta, t_total, h_total)
                        gap_new = (vmp_old + vm_old) - (vmp_new + vm_new)
                        if gap_new > gap:
                            id_ok = id_temp
                            gap = gap_new
                else:
                    id_temp[u] = m
                    num_temp[mp - 1] -= 1
                    num_temp[m - 1] += 1
                    vmp_new = func_utility(distance, chan, mp, id_temp, theta, t_total, h_total)
                    vm_new = func_utility(distance, chan, m, id_temp, theta, t_total, h_total)
                    gap_new = (vmp_old + vm_old) - (vmp_new + vm_new)
                    if gap_new > gap:
                        id_ok = id_temp
                        num_ok = num_temp
                        gap = gap_new
            id_clus_ue = id_ok
            numUE_in_SAT = num_ok

        change = id_old == id_clus_ue
        numChange_new = np.where(change == 0)[0].size
        if (numChange_new - numChange_old) == 0:
            numTimechange += 1
        else:
            numTimechange = 0
        if numTimechange > 2 or numLoop == 100:
            break
        numChange_old = numChange_new

    MODEL["IdClusUE"] = id_clus_ue
    MODEL["NumUE_in_SAT"] = numUE_in_SAT
    return MODEL


def func_utility(distance, chan, cluster, id_clus_ue, theta, t_total, h_total):
    id_ue = np.where(id_clus_ue == cluster)[0]
    value = 0.0
    for u in id_ue:
        value += theta * distance[cluster - 1, u] / 3e8 / t_total - (1 - theta) * chan[cluster - 1, u] / h_total
    return value
