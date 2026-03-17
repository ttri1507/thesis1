import matplotlib.pyplot as plt


def FUNC_Plot_model_2D(MODEL, Area, mode):
    sat = MODEL["SAT"]
    ue = MODEL["UE"]
    x_area = Area[0]
    y_area = Area[1]
    plt.plot([0, x_area, x_area, 0, 0], [0, 0, y_area, y_area, 0], "k", markersize=8, linewidth=4)
    plt.grid(True)

    if mode == "normal":
        for m in range(sat.shape[1]):
            plt.plot(sat[0, m], sat[1, m], "r-^", markersize=10, linewidth=4)
        plt.plot(ue[0, :], ue[1, :], "b^", markersize=8, linewidth=1.5)
    elif mode == "cluster":
        plt.plot(ue[0, :], ue[1, :], "b^", markersize=8, linewidth=1.5)
        id_clus_ue = MODEL["IdClusUE"]
        for u in range(ue.shape[1]):
            m = id_clus_ue[u]
            plt.plot([ue[0, u], sat[0, m - 1]], [ue[1, u], sat[1, m - 1]], "k", linewidth=2)
        for m in range(sat.shape[1]):
            plt.plot(sat[0, m], sat[1, m], "r-^", markersize=10, linewidth=4)
