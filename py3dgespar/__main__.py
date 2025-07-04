import numpy as np
from gespar3d import run_gespar3d


def main():
    dimlen = 4
    n = dimlen ** 3
    k = 3
    x = np.zeros(n, dtype=complex)
    idx = np.random.choice(n, k, replace=False)
    x[idx] = np.random.randn(k) + 1j * np.random.randn(k)
    nmse, supp_match = run_gespar3d(x, dimlen, k, n, max_t=50, snr=30, verbose=False)
    print("NMSE:", nmse)
    print("Support match:", supp_match)


if __name__ == "__main__":
    main()
