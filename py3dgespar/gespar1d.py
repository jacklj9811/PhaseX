import numpy as np
import time
from .gespar3d import objective_fun, grad_f


def best_match(x1, x2):
    """Find the best match of ``x1`` that aligns with ``x2`` via DFT ambiguities."""
    n = len(x1)
    min_err = np.inf
    x_best = x1
    for shift in range(n):
        for sign_ind in (1, -1):
            for flip in (False, True):
                shifted = np.roll(x1, shift) * sign_ind
                if flip:
                    shifted = np.flip(shifted)
                err = np.linalg.norm(x2 - shifted)
                if err < min_err:
                    x_best = shifted
                    min_err = err
    return x_best


def dgn(w, S, c, n, x0, iterations, F):
    x = np.zeros(2 * n, dtype=complex)
    x[S] = x0
    s = 0.5
    for _ in range(iterations):
        s = min(2 * s, 1)
        y = np.fft.fft(x)
        B = (np.real(y)[:, None] * np.real(F[:, S]) +
             np.imag(y)[:, None] * np.imag(F[:, S])) * np.sqrt(w)[:, None]
        b = np.sqrt(w) * (c + np.abs(y) ** 2)
        x_old = x.copy()
        f_old = objective_fun(w, c, x_old)
        x = np.zeros(2 * n, dtype=complex)
        x[S] = np.linalg.lstsq(2 * B, b, rcond=None)[0]
        x_new = x.copy()
        while objective_fun(w, c, x_old + s * (x_new - x_old)) > f_old:
            s *= 0.5
        x = x_old + s * (x_new - x_old)
        if np.linalg.norm(x - x_old) < 1e-4:
            break
    return x


def gespar_1df(c, n, k, iterations, verbose, F, ac, noisy,
               replacements_so_far, total_replacements,
               thresh, random_weights):
    if noisy:
        no_ac_supp = True
    else:
        no_ac_supp = False
    ac = ac.copy()
    ac[np.abs(ac) < 1e-8] = 0
    nz = np.nonzero(ac)[0]
    max_ac = nz.max() + 1 if nz.size else 1
    ac_off_supp = np.where(ac == 0)[0] + 1
    ac_off_supp_max = max_ac + 1 - np.where(ac == 0)[0] - 1
    ac_off_supp_max = ac_off_supp_max[(ac_off_supp_max >= 1) & (ac_off_supp_max <= n)]
    ac_off_supp = np.unique(np.concatenate([ac_off_supp, ac_off_supp_max]))
    if no_ac_supp:
        ac_off_supp = np.array([], dtype=int)
    p = np.setdiff1d(np.arange(2, n + 1), ac_off_supp)
    np.random.shuffle(p)
    replacements = 0
    w = 1 + (np.random.rand(2 * n) < 0.5) if random_weights else np.ones(2 * n)
    p = p[(p != 1) & (p != max_ac)]
    p = p[p != 0]
    k = min(k, len(p) + 2)
    supp = np.concatenate(([1, max_ac], p[:k - 2]))
    if no_ac_supp:
        supp = np.concatenate(([1], p[:k - 1]))
    supp0 = supp - 1
    x_k = dgn(w, supp0, c, n, np.random.randn(k), iterations, F)
    replacements += 1
    f_min = objective_fun(w, c, x_k)
    x_n = x_k
    while True:
        idx = np.argsort(np.abs(x_k[supp0]))
        supp0 = supp0[idx]
        f_grad = grad_f(w, c, x_k)
        off_supp = np.setdiff1d(np.arange(2 * n), supp0)
        off_supp = np.setdiff1d(off_supp, ac_off_supp - 1)
        off_supp = off_supp[np.argsort(-np.abs(f_grad[off_supp]))]
        improved = False
        for i in supp0[:1]:
            if no_ac_supp:
                if i == 0:
                    continue
            else:
                if i == 0 or i == max_ac - 1:
                    continue
            for j in off_supp[:1]:
                supp_temp = supp0.copy()
                supp_temp[supp_temp == i] = j
                x_temp = dgn(w, supp_temp, c, n, x_k[supp_temp], iterations, F)
                f_temp = objective_fun(w, c, x_temp)
                replacements += 1
                if f_temp < f_min:
                    if verbose:
                        print(
                            f"replacement: {replacements - 1} Replacing {i + 1} with {j + 1}   f= {f_temp:3.3f}")
                    x_k = x_temp
                    x_n = x_k.copy()
                    supp0 = supp_temp
                    f_min = f_temp
                    improved = True
                    if f_temp < thresh:
                        return x_n, f_min, replacements
                    break
                else:
                    x_n = x_temp
                if replacements_so_far + replacements + 1 > total_replacements:
                    return x_n, f_min, replacements
            if improved:
                break
        if not improved:
            return x_n, f_min, replacements


def run_gespar1d(x, F, n, k, maxT, snr, verbose=False):
    c = np.abs(np.fft.fft(x)) ** 2
    noise_std = 0.0 if np.isinf(snr) else np.sqrt(np.mean(c) / (10 ** (snr / 10)))
    cn = c + np.random.normal(scale=noise_std, size=c.shape)
    ac = np.fft.ifft(cn)[: n // 2]
    iterations = 100
    start = time.time()
    it_so_far = 0
    f_value_min = np.inf
    success = False
    while it_so_far < maxT:
        noisy = 1
        f_thresh = 1e-3
        random_weights = 1
        x_n, f_value, its = gespar_1df(
            cn, n // 2, k, iterations, verbose, F, ac, noisy,
            it_so_far, maxT, f_thresh, random_weights
        )
        it_so_far += its
        if f_value < 1e-4:
            runtime = time.time() - start
            success = True
            f_value_min = f_value
            x_best = x_n
            break
        if f_value < f_value_min:
            f_value_min = f_value
            x_best = x_n
    if not success:
        runtime = -(time.time() - start)
        iter_count = -it_so_far
    else:
        iter_count = it_so_far
    x_best = best_match(x_best, x)
    nmse = np.linalg.norm(x - x_best) / np.linalg.norm(x)
    supp_match = len(np.intersect1d(np.nonzero(x)[0], np.nonzero(x_best)[0]))
    return nmse, supp_match, runtime, iter_count
