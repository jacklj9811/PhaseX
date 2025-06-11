function [nmse, supp_match, runtime, iter] = run_gespar3d(x, dimlen, k, m, maxT, snr, verbose)
%RUN_GESPAR3D Perform one trial of 3D-GESPAR recovery
%c = |F x|^2 is measured and contaminated with AWGN of SNR 'snr'.
%The function returns NMSE, number of support matches, runtime and
%number of iterations used by GreedysparseRec_3d.

c = abs(fftn(reshape(x, dimlen, dimlen, dimlen))).^2;
cn = awgn(c(:), snr, 'measured');
measurementSet = 1:m;
fMin = inf;
xBest = zeros(size(x));
Tind = 0;
start = tic;
while Tind <= maxT
    [fVal, x_n, Tind] = GreedysparseRec_3d(cn, k, measurementSet, numel(x)/2, Tind, maxT, verbose);
    if fVal < fMin
        fMin = fVal;
        xBest = x_n;
        if fMin < 1e-4
            break;
        end
    end
end
runtime = toc(start);
matched = bestMatch3D(xBest, x);
matched = matched(:);
nmse = norm(x - matched)/norm(x);
supp_match = length(intersect(find(x), find(matched)));
if fMin >= 1e-4
    runtime = -runtime;
    Tind = -Tind;
end
iter = Tind;
end
