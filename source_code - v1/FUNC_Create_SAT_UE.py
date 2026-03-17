import numpy as np


def FUNC_Create_SAT_UE(Area, numSAT, numUE, zm):
    x_area = Area[0]
    y_area = Area[1]

    sat = np.zeros((3, 4), dtype=float)
    sat[0, :] = [x_area / 4, x_area * 3 / 4, x_area * 3 / 4, x_area / 4]
    sat[1, :] = [y_area / 4, y_area / 4, y_area * 3 / 4, y_area * 3 / 4]
    sat[2, :] = [zm, zm, zm, zm]

    ue = np.zeros((3, numUE), dtype=float)
    ue[0, :] = x_area * np.random.rand(numUE)
    ue[1, :] = y_area * np.random.rand(numUE)
    return sat, ue
