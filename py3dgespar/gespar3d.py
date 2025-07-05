import numpy as np
import time


def objective_fun(w, c, x):
    """Consistency error cost function."""
    return np.sum(w * ((np.abs(np.fft.fftn(x)) ** 2 - c) ** 2))


def grad_f(w, c, x):
    """Gradient of the objective function."""
    z = np.fft.fftn(x)
    return np.fft.ifftn(w * z * (np.abs(z) ** 2 - c)) * (x.size * 4)


def wg_cost_3d(c, x, w):
    dimlen = int(round(x.size ** (1 / 3)))
    y = np.abs(np.fft.fftn(x.reshape(dimlen, dimlen, dimlen))) ** 2
    return (np.sum(w * ((y.ravel() - c) ** 2))) ** (1 / 3)


def wg_grad_3d(c, x):
    dimlen = int(round(x.size ** (1 / 3)))
    z = np.fft.fftn(x.reshape(dimlen, dimlen, dimlen))
    c = c.reshape(dimlen, dimlen, dimlen)
    out = np.fft.ifftn((np.abs(z) ** 2 - c) * z)
    return out.ravel()


def gn_3d(support, c, n, x0, iterations, w):
    dimlen = int(round(n ** (1 / 3)))
    x = np.zeros(n, dtype=complex)
    x[support] = x0
    s = 0.5
    W = np.fft.fft(np.eye(dimlen))
    MM = []
    for ind_s in support:
        # ``support`` uses zero-based indexing whereas the MATLAB code assumes
        # one-based indices. Adjust the index accordingly before computing the
        # 3D coordinates used to form the DFT columns.
        ind = ind_s + 1
        alpha = int(np.ceil(ind / dimlen ** 2))
        ind -= (alpha - 1) * dimlen ** 2
        beta = int(np.ceil(ind / dimlen))
        gamma = ind % dimlen
        if gamma == 0:
            gamma = dimlen
        # Convert back to zero-based column indices for NumPy
        MM.append(
            np.kron(
                W[:, alpha - 1],
                np.kron(W[:, beta - 1], W[:, gamma - 1])
            )
        )
    MM = np.stack(MM, axis=1)
    err_vec = []
    for _ in range(iterations):
        s = min(2 * s, 1)
        z = np.fft.fftn(x.reshape(dimlen, dimlen, dimlen)).ravel()
        B = (np.real(z)[:, None] * np.real(MM) + np.imag(z)[:, None] * np.imag(MM)) * np.sqrt(w)[:, None]
        b = np.sqrt(w) * (c + np.abs(z) ** 2)
        x_old = x.copy()
        f_old = wg_cost_3d(c, x_old, w)
        x = np.zeros(n, dtype=complex)
        x[support] = np.linalg.lstsq(2 * B, b, rcond=None)[0]
        x_new = x.copy()
        while wg_cost_3d(c, x_old + s * (x_new - x_old), w) > f_old:
            s *= 0.5
        x = x_old + s * (x_new - x_old)
        err_vec.append(f_old)
        if np.linalg.norm(x - x_old) < 1e-4:
            break
    return x, np.array(err_vec)



def best_match_3d(x1, x2):
    dimlen = int(round(len(x1) ** (1 / 3)))
    x1 = x1.reshape(dimlen, dimlen, dimlen)
    x2 = x2.reshape(dimlen, dimlen, dimlen)
    min_err = np.inf
    x_best = x1
    phase_vec = np.linspace(0, np.pi, 2)
    for kk in range(dimlen):
        for jj in range(dimlen):
            for qq in range(dimlen):
                for phase in phase_vec:
                    for flip in [False, True]:
                        shifted = np.roll(x1, (kk, jj, qq), axis=(0, 1, 2)) * np.exp(1j * phase)
                        if flip:
                            shifted = np.flip(np.flip(np.flip(shifted, 0), 1), 2)
                        dis = x2 - shifted
                        err = np.linalg.norm(dis.ravel())
                        if err < min_err:
                            x_best = shifted
                            min_err = err
    return x_best.ravel()


def greedysparse_rec_3d(c, k, measurement_set, n, tind, max_t, verbose=False):
    p = np.random.permutation(n)
    supp = p[:k]
    iterations = 1000
    c = c[measurement_set]
    w = (1 + (np.random.rand(len(measurement_set)) < 0.5)).astype(float)
    # ``n`` corresponds to the total signal length. The original MATLAB code
    # doubled ``n`` when calling ``GN_3d`` since it assumed an even cube size
    # and used ``numel(x)/2``.  When ``n`` is already the actual length (for
    # odd cube sizes), doubling it results in a reshape error.  Use the given
    # ``n`` directly so both even and odd dimensions are supported.
    x_k, _ = gn_3d(
        supp,
        c,
        n,
        np.random.randn(k) + 1j * np.random.rand(k),
        iterations,
        w,
    )
    f_min = wg_cost_3d(c, x_k, w)
    while True:
        supp = supp[np.argsort(np.abs(x_k[supp]))]
        f_grad = wg_grad_3d(c, x_k)
        off_supp = np.setdiff1d(np.arange(n), supp)
        off_supp = off_supp[np.argsort(-np.abs(f_grad[off_supp]))]
        improved = False
        for i in supp:
            for j in off_supp[:1]:
                supp_temp = supp.copy()
                supp_temp[supp_temp == i] = j
                tind += 1
                x_temp, _ = gn_3d(supp_temp, c, n, x_k[supp_temp], iterations, w)
                f_temp = wg_cost_3d(c, x_temp, w)
                if f_temp < f_min:
                    if verbose:
                        print(f"it: ?, T: {tind}, Replaced {i} with {j}   f= {f_temp:3.3f}")
                    x_k = x_temp
                    supp = supp_temp
                    improved = True
                    f_min = f_temp
                    if f_temp < 1e-3:
                        return f_min, x_k, tind
                    break
            if improved:
                break
        if not improved or tind > max_t:
            if verbose:
                print("no possible improvement - trying new initial guess")
            return f_min, x_k, tind


def run_gespar3d(x, dimlen, k, m, max_t, snr, verbose=False):
    c = np.abs(np.fft.fftn(x.reshape(dimlen, dimlen, dimlen))) ** 2
    cn = c + np.random.normal(scale=np.sqrt(np.mean(c) / (10 ** (snr / 10))), size=c.shape)
    measurement_set = np.arange(m)
    f_min = np.inf
    x_best = np.zeros_like(x)
    t_ind = 0
    n = x.size
    while t_ind <= max_t:
        f_val, x_n, t_ind = greedysparse_rec_3d(
            cn.ravel(),
            k,
            measurement_set,
            n,
            t_ind,
            max_t,
            verbose,
        )
        if f_val < f_min:
            f_min = f_val
            x_best = x_n
            if f_min < 1e-4:
                break
    matched = best_match_3d(x_best, x)
    nmse = np.linalg.norm(x - matched) / np.linalg.norm(x)
    supp_match = len(np.intersect1d(np.nonzero(x)[0], np.nonzero(matched)[0]))
    return nmse, supp_match, matched


def run_gespar3d_stats(x, dimlen, k, m, max_t, snr, verbose=False):
    """Version of ``run_gespar3d`` that also reports runtime and iterations."""
    start = time.time()
    c = np.abs(np.fft.fftn(x.reshape(dimlen, dimlen, dimlen))) ** 2
    cn = c + np.random.normal(scale=np.sqrt(np.mean(c) / (10 ** (snr / 10))),
                              size=c.shape)
    measurement_set = np.arange(m)
    f_min = np.inf
    x_best = np.zeros_like(x)
    t_ind = 0
    n = x.size
    while t_ind <= max_t:
        f_val, x_n, t_ind = greedysparse_rec_3d(
            cn.ravel(),
            k,
            measurement_set,
            n,
            t_ind,
            max_t,
            verbose,
        )
        if f_val < f_min:
            f_min = f_val
            x_best = x_n
            if f_min < 1e-4:
                break
    runtime = time.time() - start
    matched = best_match_3d(x_best, x)
    nmse = np.linalg.norm(x - matched) / np.linalg.norm(x)
    supp_match = len(np.intersect1d(np.nonzero(x)[0], np.nonzero(matched)[0]))
    if f_min >= 1e-4:
        runtime = -runtime
        t_ind = -t_ind
    return nmse, supp_match, runtime, t_ind
