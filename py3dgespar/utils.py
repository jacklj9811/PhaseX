import numpy as np


def compute_stats(t_array, iter_array, nmse_array):
    """Calculate summary statistics for an experiment."""
    recovery_rate = np.sum(t_array > 0) / len(t_array)
    mean_t = np.mean(np.abs(t_array))
    std_t = np.std(np.abs(t_array))
    mean_iter = np.mean(np.abs(iter_array))
    std_iter = np.std(np.abs(iter_array))
    nmse = np.mean(nmse_array)
    return np.array([recovery_rate, mean_t, std_t, mean_iter, std_iter, nmse])
