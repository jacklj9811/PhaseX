function out = compute_stats(t_array, iter_array, nmse_array)
%COMPUTE_STATS Calculate summary statistics for the experiment

recovery_rate = sum(t_array > 0) / numel(t_array);
mean_t = mean(abs(t_array));
std_t = std(abs(t_array));
mean_iter = mean(abs(iter_array));
std_iter = std(abs(iter_array));
nmse = mean(nmse_array);

out = [recovery_rate, mean_t, std_t, mean_iter, std_iter, nmse];
end
