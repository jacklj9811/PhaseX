% This does not work.
% This scripts performs a 2D Fourier GESPAR recovery example.
% It generates a random real matrix and tries to recover it from its (possibly noisy) 2D Fourier magnitude measurements. 
% 
% Important parameters:
% n= number of measurements.  The size of the matrix is sqrt(n) by sqrt(n)
% kVec = simulated sparsity level
% maxT = Number of replacements allowed per signal
% snr = noise added to measurements.  
close all;
clear;
% s = RandStream('mcg16807','Seed',0); % Seeding random number generation
% RandStream.setDefaultStream(s);

% Seeding random generator
s = RandStream('mcg16807','Seed',0);
% RandStream.setDefaultStream(s);
RandStream.setGlobalStream(s);
% n=32^2; % size of x
%\/
% dimlen=6;
dimlen=20;
n=dimlen^3;
%/\
m=n; % number of measurement vectors
snr=inf; 
% maxIt=1000;
% maxT=6400; % Total number of replacements
maxT=100; % Total number of replacements
display=1;
kind=0;
k=55; % Sparsity (number of nonzeros in x)
% k=4;
itsPer_kMax=1;
tic
x=zeros(n,1); 
locs=randperm(n);
x(locs(1:k))=(1+rand(k,1)).*(-1).^randi(2,[k,1]);% signal vector
% c=abs(fft2(reshape(x,sqrt(n),sqrt(n)))).^2;% ideal measurements
c=abs(fftn(reshape(x,dimlen,dimlen,dimlen))).^2;% ideal measurements
c=awgn(c(:),snr,'measured'); % noisy measurements
measurementSet=1:m; % Set of measurements (rows in M) used. 1:m uses all measurements 
fMin=inf;
xBest=zeros(size(x));
Tind=0;
verbose=1;
while Tind<=maxT      
    [fVal,x_n,Tind]=GreedysparseRec_3d(c,k,measurementSet,n,Tind,maxT,verbose);
    if fVal<fMin
        fMin=fVal;
        xBest=x_n;
        % xBestFlipped=bestMatch2D(xBest,x);
        if fMin<norm(c)/1e5 
                break;
        end
    end    
end
% err=1-abs(x'*xBestFlipped(:)/(norm(x)*norm(xBestFlipped(:))));
% fprintf('Error =%d\n',err);
% xBestFlipped=xBestFlipped(:);
t=toc;
% figure;subplot(2,2,[3:4]);stem(1:n,real(x),'--o');hold all;subplot(2,2,[3:4]);plot(1:n,real(xBestFlipped),'rx');xlabel('index');ylabel('value');xlim([0 n]);title('True/Recovered (column stacked)');legend('True signal','Recovered');
% subplot(2,2,1);imagesc(reshape(c/max(c),sqrt(n),sqrt(n)));colorbar;title('Measurements (Fourier magnitude)');colormap(bone);subplot(2,2,2);imagesc(reshape(x,sqrt(n),sqrt(n)));colormap(bone);title('True signal');colorbar;
