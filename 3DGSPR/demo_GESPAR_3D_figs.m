% This works.
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
warning('off','all');
% s = RandStream('mcg16807','Seed',0); % Seeding random number generation
% RandStream.setDefaultStream(s);

% Seeding random generator
s = RandStream('mcg16807','Seed',0);
% RandStream.setDefaultStream(s);
RandStream.setGlobalStream(s);
% n=32^2; % size of x
% snr=20; 
snr=inf; 
% maxT=6400; % Total number of replacements
% maxT=20;
% maxT=50;
maxT=1000;
display=1;
kind=0;
% k=20; % Sparsity (number of nonzeros in x)
% epoches=100;
epoches=30;
nmse_array=zeros(epoches,1);
n_supp_match=zeros(epoches,1);
t_array=zeros(epoches,1);
n_iter_array=zeros(epoches,1);
fprintf('epoches=%d\n',epoches);
dimlen=6;
n=dimlen^3;
% m=n; % number of measurement vectors
m=n; % number of measurement vectors
fprintf('dimlen=%d\n',dimlen);
% k=4;
% ks = 2:2:8;
% ks = [2:1:8];
% ks = [2:2:20];
ks = [2:2:26];
hist_all = zeros(length(ks),6);
nmse_array_1d=zeros(epoches,1);
n_supp_match_1d=zeros(epoches,1);
t_array_1d=zeros(epoches,1);
n_iter_array_1d=zeros(epoches,1);
hist_all_1d = hist_all;
count=1;
verbose=0;
% gespar_1d=1;
gespar_1d=0;
%
F=fft(eye(n));
%
for k = ks
	% if k > 30
	% 	maxT=100;
	% end
	for epoch = 1:epoches
		itsPer_kMax=1;
		x=zeros(n,1); 
		% locs=randperm(n);
		locs=randperm(n/2); % oversampling
		x(locs(1:k))=(1+rand(k,1)).*(-1).^randi(2,[k,1]);% signal vector
		% c=abs(fft2(reshape(x,sqrt(n),sqrt(n)))).^2;% ideal measurements

		% 3d-GESPAR
		c=abs(fftn(reshape(x,dimlen,dimlen,dimlen))).^2;% ideal measurements
		c=c(:);
		% cn=c(:);
		cn=awgn(c(:),snr,'measured'); % noisy measurements
		measurementSet=1:m; % Set of measurements (rows in M) used. 1:m uses all measurements 
		fMin=inf;
		xBest=zeros(size(x));
		Tind=0;
		t3d=tic;
		while Tind<=maxT
		    [fVal,x_n,Tind]=GreedysparseRec_3d(cn,k,measurementSet,n/2,Tind,maxT,verbose);
            if fVal<fMin
		        fMin=fVal;
		        xBest=x_n;
		        % if fMin<norm(c)/1e5
		        % if (fMin<2*norm(c-cn)^2 || fMin<1e-4)
		        if fMin<1e-4
		                break;
		        end
		    end    
        end
		t=toc(t3d);

		xBestFlipped=bestMatch3D(xBest,x);
		xBestFlipped=xBestFlipped(:);
		nmse_array(epoch)=norm(x-xBestFlipped)/norm(x);
		n_supp_match(epoch)=length(intersect(find(x),find(xBestFlipped)));
		% if fMin<norm(c)/1e5
		% if (fMin<2*norm(c-cn)^2 || fMin<1e-4) 
		if fMin<1e-4
			t_array(epoch) = t;
			n_iter_array(epoch) = Tind;
		else
			% t_array(epoch) = Inf;
			% n_iter_array(epoch) = Inf;			
			t_array(epoch) = -t;
			n_iter_array(epoch) = -Tind;
		end

		if gespar_1d == 1
			% GESPAR_1D
			% snr=1001; % snr>1000 is treated as noiseless
	    	c=abs(fft(x)).^2;% ideal measurements
	        % F=fft(eye(n));
	        success=0;
	        iterations=100;
	        % cn=c;
	        cn=awgn(c,snr,'measured'); % noised c
	        t1d=tic;%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        % cn=0.5*(cn+[cn(1);cn(end:-1:2)]);
	        ac=(ifft(cn));
	        ac=ac(1:n/2);
	        % ac=ac(1:n);
	        %%
	        itSoFar=0;
	        trial=0;
	        totIterations=maxT;
	        fValueMin=inf;
	        while (itSoFar<totIterations)
	            if (verbose)
	                fprintf('total replacements so far = %d \n',itSoFar);
	            end
	            %% using GESPAR to recover x from cn (noisy measurements)
	           noisy=1;
	           % noisy=0;
	           fThresh=1e-3;
	           randomWeights=1;
	           % randomWeights=0;
	           [x_n,fValue,its]=GESPAR_1DF(cn,n/2,k,iterations,verbose,F,ac,noisy,itSoFar,totIterations,fThresh,randomWeights);
	           % [x_n,fValue,its]=GESPAR_1DF(cn,n,k,iterations,verbose,F,ac,noisy,itSoFar,totIterations,fThresh,randomWeights);
	           trial=trial+1;
	           itSoFar=itSoFar+its;
	           if fValue<1e-4
	           % if (fValue<2*norm(c-cn)^2 || fValue<1e-4) % Breaking condition for a successful recovery
	               % t=round(toc*100)/100;
	               t=toc(t1d);
	               success=1;
	               fValueMin=fValue;
	               x_n_best=x_n;
				   t_array_1d(epoch) = t;
				   n_iter_array_1d(epoch) = itSoFar;
	               % fprintf('%d. succeeded! k = %d   total evaluations %d  in %d initial points took %2.2f secs\n',it,k,itSoFar,trial,t);
	             break;
	           end
	            if (fValue<fValueMin)
	               fValueMin=fValue;
	               x_n_best=x_n;
	           end

	        end
	        if (~success)
	            % t=round(toc*100)/100;
	            t=toc(t1d);
	   			% t_array_1d(epoch) = Inf;
				% n_iter_array_1d(epoch) = Inf;
				t_array_1d(epoch) = -t;
				n_iter_array_1d(epoch) = -itSoFar;
	        end
	        x_n_best=bestMatch(x_n_best,x);
	        nmse_array_1d(epoch)=norm(x-x_n_best)/norm(x);
			n_supp_match_1d(epoch)=length(intersect(find(x),find(x_n_best)));
		end
	end
	% recovery_rate = length(find(t_array<Inf)) / length(t_array);
	% mean_t_of_success = mean(t_array(find(t_array<Inf)));
	% std_t_of_success = std(t_array(find(t_array<Inf)));
	% mean_niter_of_success = mean(n_iter_array(find(n_iter_array<Inf)));
	% std_niter_of_success = std(n_iter_array(find(n_iter_array<Inf)));
	recovery_rate = length(find(t_array>0)) / length(t_array);
	mean_t_of_success = mean(abs(t_array));
	std_t_of_success = std(abs(t_array));
	mean_niter_of_success = mean(abs(n_iter_array));
	std_niter_of_success = std(abs(n_iter_array));
	nmse=mean(nmse_array);
	mean_supp_match=mean(n_supp_match)/length(find(x));
	is_supp_match=mean(n_supp_match==length(find(x)));
	hist_all(count,:)=[recovery_rate,mean_t_of_success,std_t_of_success,mean_niter_of_success,std_niter_of_success,nmse];
	fprintf('method: 3d-GESPAR\n');
	fprintf('k=%d, recovery rate=%.2f\nWhen succeed, mean t=%.4f secs, std t=%.4f, mean iteration=%.2f, std iteration=%.2f\n',k,recovery_rate,mean_t_of_success,std_t_of_success,mean_niter_of_success,std_niter_of_success);

	if gespar_1d == 1
		% recovery_rate = length(find(t_array_1d<Inf)) / length(t_array_1d);
		% mean_t_of_success = mean(t_array_1d(find(t_array_1d<Inf)));
		% std_t_of_success = std(t_array_1d(find(t_array_1d<Inf)));
		% mean_niter_of_success = mean(n_iter_array_1d(find(n_iter_array_1d<Inf)));
		% std_niter_of_success = std(n_iter_array_1d(find(n_iter_array_1d<Inf)));
		recovery_rate = length(find(t_array_1d>0)) / length(t_array_1d);
		mean_t_of_success = mean(abs(t_array_1d));
		std_t_of_success = std(abs(t_array_1d));
		mean_niter_of_success = mean(abs(n_iter_array_1d));
		std_niter_of_success = std(abs(n_iter_array_1d));
		nmse=mean(nmse_array_1d);
		hist_all_1d(count,:)=[recovery_rate,mean_t_of_success,std_t_of_success,mean_niter_of_success,std_niter_of_success,nmse];
		fprintf('method: GESPAR-1D\n');
		fprintf('k=%d, recovery rate=%.2f\nWhen succeed, mean t=%.4f secs, std t=%.4f, mean iteration=%.2f, std iteration=%.2f\n',k,recovery_rate,mean_t_of_success,std_t_of_success,mean_niter_of_success,std_niter_of_success);
	else
		fprintf('method: GESPAR-1D, too weak to be displayed.\n');
	end
	% if recovery_rate < 0.03
	% 	gespar_1d=0;
	% 	epoches=10;
	% end

	count=count+1;
end
if gespar_1d == 1
	figure(1);
	plot(ks,hist_all(:,1),'-o',ks,hist_all_1d(:,1),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Recovery rate','FontSize',15,'FontWeight','bold');
	ylim([0 1]);
	legend({'3D-GSPR','GESPAR'},'Location','southwest','FontSize',15,'FontWeight','bold');

	figure(2);
	plot(ks,hist_all(:,4),'-o',ks,hist_all_1d(:,4),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Mean # iterations','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR','GESPAR'},'Location','northwest','FontSize',15,'FontWeight','bold');

	figure(3);
	plot(ks,hist_all(:,5),'-o',ks,hist_all_1d(:,5),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Std # iteration','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR','GESPAR'},'Location','northwest','FontSize',15,'FontWeight','bold');

	figure(4);
	plot(ks,hist_all(:,2),'-o',ks,hist_all_1d(:,2),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Mean running time[sec]','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR','GESPAR'},'Location','northwest','FontSize',15,'FontWeight','bold');

	figure(5);
	plot(ks,hist_all(:,6),'-o',ks,hist_all_1d(:,6),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('NMSE','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR','GESPAR'},'Location','northwest','FontSize',15,'FontWeight','bold');
else
	figure(1);
	plot(ks,hist_all(:,1),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Recovery rate','FontSize',15,'FontWeight','bold');
	ylim([0 1]);
	legend({'3D-GSPR'},'Location','southwest','FontSize',15,'FontWeight','bold');

	figure(2);
	plot(ks,hist_all(:,4),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Mean # iterations','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR'},'Location','northwest','FontSize',15,'FontWeight','bold');

	figure(3);
	plot(ks,hist_all(:,5),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Std # iteration','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR'},'Location','northwest','FontSize',15,'FontWeight','bold');

	figure(4);
	plot(ks,hist_all(:,2),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('Mean running time[sec]','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR'},'Location','northwest','FontSize',15,'FontWeight','bold');

	figure(5);
	plot(ks,hist_all(:,6),'-o','LineWidth',1.5,'MarkerSize',12);
	set(gca,'FontSize',15,'FontWeight','bold');
	xlabel('Sparsity','FontSize',15,'FontWeight','bold');
	ylabel('NMSE','FontSize',15,'FontWeight','bold');
	legend({'3D-GSPR'},'Location','northwest','FontSize',15,'FontWeight','bold');
end
warning('on','all');




