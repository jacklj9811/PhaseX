% Demo for 3D-GESPAR with optional 1D GESPAR comparison
close all; clear;
warning('off','all');

%% Random seed for reproducibility
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

%% Parameters
snr    = inf;   % measurement SNR
maxT   = 1000;  % maximal number of replacements
epoches = 30;   % number of trials per sparsity level
dimlen = 6;     % cube dimension
n      = dimlen^3;
m      = n;     % number of measurement vectors
ks     = 2:2:26;
gespar_1d = 1;  % set to 1 to run the 1D variant
verbose   = 0;
F = fft(eye(n));

hist_all    = zeros(numel(ks),6);
hist_all_1d = hist_all;

for idx = 1:numel(ks)
    k = ks(idx);

    nmse      = zeros(epoches,1);
    supp      = nmse;
    t_arr     = nmse;
    iter_arr  = nmse;

    nmse1     = nmse;
    supp1     = nmse;
    t_arr1    = nmse;
    iter_arr1 = nmse;

    for e = 1:epoches
        %% Generate sparse signal
        x = zeros(n,1);
        locs = randperm(n/2); % oversampling
        x(locs(1:k)) = (1+rand(k,1)).*(-1).^randi(2,[k,1]);

        [nmse(e), supp(e), t_arr(e), iter_arr(e)] = ...
            run_gespar3d(x, dimlen, k, m, maxT, snr, verbose);

        if gespar_1d
            [nmse1(e), supp1(e), t_arr1(e), iter_arr1(e)] = ...
                run_gespar1d(x, F, n, k, maxT, snr, verbose);
        end
    end

    hist_all(idx,:) = compute_stats(t_arr, iter_arr, nmse);
    if gespar_1d
        hist_all_1d(idx,:) = compute_stats(t_arr1, iter_arr1, nmse1);
    end

    fprintf('method: 3d-GESPAR\n');
    fprintf(['k=%d, recovery rate=%.2f\nWhen succeed, mean t=%.4f secs, ' ...
             'std t=%.4f, mean iteration=%.2f, std iteration=%.2f\n'], ...
            k, hist_all(idx,1), hist_all(idx,2), hist_all(idx,3), ...
            hist_all(idx,4), hist_all(idx,5));

    if gespar_1d
        fprintf('method: GESPAR-1D\n');
        fprintf(['k=%d, recovery rate=%.2f\nWhen succeed, mean t=%.4f secs, ' ...
                 'std t=%.4f, mean iteration=%.2f, std iteration=%.2f\n'], ...
                k, hist_all_1d(idx,1), hist_all_1d(idx,2), hist_all_1d(idx,3), ...
                hist_all_1d(idx,4), hist_all_1d(idx,5));
    end
end

%% Plot statistics
if gespar_1d
    figure(1); plot(ks,hist_all(:,1),'-o',ks,hist_all_1d(:,1),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Recovery rate'); ylim([0 1]);
    legend({'3D-GSPR','GESPAR'},'Location','southwest');

    figure(2); plot(ks,hist_all(:,4),'-o',ks,hist_all_1d(:,4),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Mean # iterations');
    legend({'3D-GSPR','GESPAR'},'Location','northwest');

    figure(3); plot(ks,hist_all(:,5),'-o',ks,hist_all_1d(:,5),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Std # iteration');
    legend({'3D-GSPR','GESPAR'},'Location','northwest');

    figure(4); plot(ks,hist_all(:,2),'-o',ks,hist_all_1d(:,2),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Mean running time[sec]');
    legend({'3D-GSPR','GESPAR'},'Location','northwest');

    figure(5); plot(ks,hist_all(:,6),'-o',ks,hist_all_1d(:,6),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('NMSE');
    legend({'3D-GSPR','GESPAR'},'Location','northwest');
else
    figure(1); plot(ks,hist_all(:,1),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Recovery rate'); ylim([0 1]);
    legend({'3D-GSPR'},'Location','southwest');

    figure(2); plot(ks,hist_all(:,4),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Mean # iterations');
    legend({'3D-GSPR'},'Location','northwest');

    figure(3); plot(ks,hist_all(:,5),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Std # iteration');
    legend({'3D-GSPR'},'Location','northwest');

    figure(4); plot(ks,hist_all(:,2),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('Mean running time[sec]');
    legend({'3D-GSPR'},'Location','northwest');

    figure(5); plot(ks,hist_all(:,6),'-o','LineWidth',1.5,'MarkerSize',12);
    set(gca,'FontSize',15,'FontWeight','bold'); xlabel('Sparsity'); ylabel('NMSE');
    legend({'3D-GSPR'},'Location','northwest');
end

warning('on','all');
