import numpy as np
import matplotlib.pyplot as plt
from py3dgespar import run_gespar3d_stats, run_gespar1d, compute_stats


def main():
    snr = np.inf
    max_t = 1000
    epoches = 30
    # ``run_gespar1d`` only supports even signal lengths.  The demo defaults to
    # an even cube dimension so the flattened signal length is even as well.
    dimlen = 6
    n = dimlen ** 3
    m = n
    ks = np.arange(2, 27, 2)
    gespar_1d = True
    if n % 2 == 1:
        print("WARNING: odd `dimlen` detected; skipping GESPAR-1D which requires an even length")
        gespar_1d = False
    verbose = False
    F = np.fft.fft(np.eye(n)) if gespar_1d else None

    hist_all = np.zeros((len(ks), 6))
    hist_all_1d = np.zeros_like(hist_all)

    for idx, k in enumerate(ks):
        nmse = np.zeros(epoches)
        supp = np.zeros(epoches)
        t_arr = np.zeros(epoches)
        iter_arr = np.zeros(epoches)

        nmse1 = np.zeros(epoches)
        supp1 = np.zeros(epoches)
        t_arr1 = np.zeros(epoches)
        iter_arr1 = np.zeros(epoches)

        for e in range(epoches):
            x = np.zeros(n, dtype=complex)
            locs = np.random.permutation(n // 2)
            x[locs[:k]] = (1 + np.random.rand(k)) * (-1) ** np.random.randint(1, 3, k)

            nmse[e], supp[e], t_arr[e], iter_arr[e] = run_gespar3d_stats(
                x, dimlen, k, m, max_t, snr, verbose
            )

            if gespar_1d:
                nmse1[e], supp1[e], t_arr1[e], iter_arr1[e] = run_gespar1d(
                    x, F, n, k, max_t, snr, verbose
                )

        hist_all[idx, :] = compute_stats(t_arr, iter_arr, nmse)
        if gespar_1d:
            hist_all_1d[idx, :] = compute_stats(t_arr1, iter_arr1, nmse1)

        print("method: 3d-GESPAR")
        print(
            f"k={k}, recovery rate={hist_all[idx,0]:.2f}\n"
            f"When succeed, mean t={hist_all[idx,1]:.4f} secs, std t={hist_all[idx,2]:.4f},"
            f" mean iteration={hist_all[idx,3]:.2f}, std iteration={hist_all[idx,4]:.2f}"
        )
        if gespar_1d:
            print("method: GESPAR-1D")
            print(
                f"k={k}, recovery rate={hist_all_1d[idx,0]:.2f}\n"
                f"When succeed, mean t={hist_all_1d[idx,1]:.4f} secs, std t={hist_all_1d[idx,2]:.4f},"
                f" mean iteration={hist_all_1d[idx,3]:.2f}, std iteration={hist_all_1d[idx,4]:.2f}"
            )

    if gespar_1d:
        plt.figure(1)
        plt.plot(ks, hist_all[:, 0], "-o", label="3D-GSPR")
        plt.plot(ks, hist_all_1d[:, 0], "-o", label="GESPAR")
        plt.xlabel("Sparsity")
        plt.ylabel("Recovery rate")
        plt.ylim([0, 1])
        plt.legend(loc="lower right")

        plt.figure(2)
        plt.plot(ks, hist_all[:, 3], "-o", label="3D-GSPR")
        plt.plot(ks, hist_all_1d[:, 3], "-o", label="GESPAR")
        plt.xlabel("Sparsity")
        plt.ylabel("Mean # iterations")
        plt.legend(loc="upper left")

        plt.figure(3)
        plt.plot(ks, hist_all[:, 4], "-o", label="3D-GSPR")
        plt.plot(ks, hist_all_1d[:, 4], "-o", label="GESPAR")
        plt.xlabel("Sparsity")
        plt.ylabel("Std # iteration")
        plt.legend(loc="upper left")

        plt.figure(4)
        plt.plot(ks, hist_all[:, 1], "-o", label="3D-GSPR")
        plt.plot(ks, hist_all_1d[:, 1], "-o", label="GESPAR")
        plt.xlabel("Sparsity")
        plt.ylabel("Mean running time[sec]")
        plt.legend(loc="upper left")

        plt.figure(5)
        plt.plot(ks, hist_all[:, 5], "-o", label="3D-GSPR")
        plt.plot(ks, hist_all_1d[:, 5], "-o", label="GESPAR")
        plt.xlabel("Sparsity")
        plt.ylabel("NMSE")
        plt.legend(loc="upper left")
    else:
        plt.figure(1)
        plt.plot(ks, hist_all[:, 0], "-o", label="3D-GSPR")
        plt.xlabel("Sparsity")
        plt.ylabel("Recovery rate")
        plt.ylim([0, 1])
        plt.legend(loc="lower right")

        plt.figure(2)
        plt.plot(ks, hist_all[:, 3], "-o", label="3D-GSPR")
        plt.xlabel("Sparsity")
        plt.ylabel("Mean # iterations")
        plt.legend(loc="upper left")

        plt.figure(3)
        plt.plot(ks, hist_all[:, 4], "-o", label="3D-GSPR")
        plt.xlabel("Sparsity")
        plt.ylabel("Std # iteration")
        plt.legend(loc="upper left")

        plt.figure(4)
        plt.plot(ks, hist_all[:, 1], "-o", label="3D-GSPR")
        plt.xlabel("Sparsity")
        plt.ylabel("Mean running time[sec]")
        plt.legend(loc="upper left")

        plt.figure(5)
        plt.plot(ks, hist_all[:, 5], "-o", label="3D-GSPR")
        plt.xlabel("Sparsity")
        plt.ylabel("NMSE")
        plt.legend(loc="upper left")

    plt.show()


if __name__ == "__main__":
    main()
