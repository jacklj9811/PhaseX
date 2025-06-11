function [nmse, supp_match, runtime, iter] = run_gespar1d(x, F, n, k, maxT, snr, verbose)
%RUN_GESPAR1D Perform one trial of GESPAR (1D version)

c = abs(fft(x)).^2;
cn = awgn(c, snr, 'measured');
ac = ifft(cn);
ac = ac(1:n/2);
iterations = 100;
start = tic;
itSoFar = 0;
trial = 0;
fValueMin = inf;
success = false;
while itSoFar < maxT
    noisy = 1;
    fThresh = 1e-3;
    randomWeights = 1;
    [x_n, fValue, its] = GESPAR_1DF(cn, n/2, k, iterations, verbose, F, ac, noisy, itSoFar, maxT, fThresh, randomWeights);
    trial = trial + 1;
    itSoFar = itSoFar + its;
    if fValue < 1e-4
        runtime = toc(start);
        success = true;
        fValueMin = fValue;
        x_best = x_n;
        break;
    end
    if fValue < fValueMin
        fValueMin = fValue;
        x_best = x_n;
    end
end
if ~success
    runtime = -toc(start);
    iter = -itSoFar;
else
    iter = itSoFar;
end
x_best = bestMatch(x_best, x);
nmse = norm(x - x_best)/norm(x);
supp_match = length(intersect(find(x), find(x_best)));
end
