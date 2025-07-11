import argparse
import numpy as np
import matplotlib.pyplot as plt
from .gespar3d import run_gespar3d


def plot_volume(true_x, est_x, dimlen):
    coords = np.indices((dimlen, dimlen, dimlen)).reshape(3, -1)
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    mask1 = np.abs(true_x) > 1e-8
    ax1.scatter(coords[0][mask1], coords[1][mask1], coords[2][mask1],
                c=np.abs(true_x)[mask1], cmap="viridis")
    ax1.set_title("True")
    mask2 = np.abs(est_x) > 1e-8
    ax2.scatter(coords[0][mask2], coords[1][mask2], coords[2][mask2],
                c=np.abs(est_x)[mask2], cmap="viridis")
    ax2.set_title("Recovered")
    plt.tight_layout()
    plt.show()


def plot_flattened(true_x, est_x):
    """Compare flattened amplitudes in a single 1D plot."""
    true_flat = np.abs(true_x).ravel()
    est_flat = np.abs(est_x).ravel()
    idx = np.arange(len(true_flat))
    fig, ax = plt.subplots()
    ax.scatter(idx, true_flat, color="blue", label="True", marker=".")
    ax.scatter(idx, est_flat, facecolors="none", edgecolors="red",
               label="Recovered", marker="o")
    ax.set_xlabel("Index")
    ax.set_ylabel("Magnitude")
    ax.legend()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Demo for the Python 3D-GESPAR port")
    parser.add_argument("--plot", action="store_true",
                        help="show a 1D comparison plot of true vs recovered")
    args = parser.parse_args()
    dimlen = 6
    n = dimlen ** 3
    k = 3
    x = np.zeros(n, dtype=complex)
    idx = np.random.choice(n, k, replace=False)
    x[idx] = np.random.randn(k) + 1j * np.random.randn(k)
    nmse, supp_match, recovered = run_gespar3d(x, dimlen, k, n, max_t=50, snr=30, verbose=False)
    print("NMSE:", nmse)
    print("Support match:", supp_match)
    if args.plot:
        plot_flattened(x, recovered)


if __name__ == "__main__":
    main()
