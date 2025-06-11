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
% s = RandStream('mcg16807','Seed',0); % Seeding random number generation
% RandStream.setDefaultStream(s);

% Seeding random generator
s = RandStream('mcg16807','Seed',0);
% RandStream.setDefaultStream(s);
RandStream.setGlobalStream(s);
% n=32^2; % size of x
%\/
% dimlen=6;
dimlen=40;
n=dimlen^3;
%/\
m=n; % number of measurement vectors
snr=inf; 
% maxIt=1000;
maxT=6400; % Total number of replacements
% maxT=100; % Total number of replacements
display=1;
kind=0;
k=92; % Sparsity (number of nonzeros in x)
%
gitar = [5 4 4 1.0;
3 3 3 1.0;
4 4 7 1.0;
 3 4 12 1.0;
 4 4 9 1.0;
 5 4 0 1.0;
 2 4 0 1.0;
 6 3 2 1.0;
 3 4 2 1.0;
 2 3 0 1.0;
 2 4 5 1.0;
 4 3 5 1.0;
 4 4 14 1.0;
 4 3 1 1.0;
 1 4 2 1.0;
 3 4 4 1.0;
 5 4 2 1.0;
 4 4 11 1.0;
 5 4 6 1.0;
 3 4 0 1.0;
 5 3 3 1.0;
 4 4 5 1.0;
 6 4 1 1.0;
 2 4 3 1.0;
 3 3 6 1.0;
 4 4 8 1.0;
 4 3 0 1.0;
 3 3 4 1.0;
 2 3 1 1.0;
 6 3 5 1.0;
 6 4 3 1.0;
 2 4 1 1.0;
 5 3 1 1.0;
 3 3 2 1.0;
 4 3 6 1.0;
 4 4 3 1.0;
 4 4 13 1.0;
 5 4 5 1.0;
 3 4 6 1.0;
 5 3 4 1.0;
 6 4 2 1.0;
 5 3 2 1.0;
 4 3 4 1.0;
 3 3 0 1.0;
 4 4 12 1.0;
 4 3 2 1.0;
 2 3 2 1.0;
3 4 7 1.0;
 1 4 1 1.0;
 3 4 14 1.0;
3 3 5 1.0;
 4 4 0 1.0;
 4 4 1 1.0;
 3 4 3 1.0;
 4 4 2 1.0;
 5 3 0 1.0;
 2 3 4 1.0;
 5 4 1 1.0;
 6 3 3 1.0;
 3 4 5 1.0;
 4 4 6 1.0;
 4 4 10 1.0;
 5 4 3 1.0;
 2 3 3 1.0;
 2 4 4 1.0;
 6 3 1 1.0;
 5 3 5 1.0;
 3 4 1 1.0;
 2 3 6 1.0;
 6 4 0 1.0;
 4 3 3 1.0;
 2 3 5 1.0;
 2 4 2 1.0;
 4 4 4 1.0;
3 3 1 1.0;
 2 4 6 1.0;
 3 4 8 1.0;
6 3 0 1.0;
 5 3 6 1.0;
 1 3 1 1.0;
 3 4 10 1.0;
 3 4 13 1.0;
 6 4 5 1.0;
 6 4 6 1.0;
 6 4 4 1.0;
 1 3 2 1.0;
 4 4 15 1.0;
 6 3 4 1.0;
 3 4 9 1.0;
 1 4 0 1.0;
3 4 11 1.0;
 1 3 0 1.0];
gitar_35_idxs=gitar(:,1)*dimlen*dimlen+gitar(:,2)*dimlen+gitar(:,3)+1;
one=[5 5 5  ;
25  25  25  ;
10	15	14	;
11	14	14	;
11	15	14	;
12	13	14	;
12	14	14	;
12	15	14	;
13	12	14	;
13	13	14	;
13	14	14	;
13	15	14	;
14	14	14	;
15	14	14	;
16	14	14	;
17	14	14	;
14	15	14	;
15	15	14	;
16	15	14	;
17	15	14	;
18	12	14	;
18	13	14	;
18	14	14	;
18	15	14	;
18	16	14	;
18	17	14	;
19	12	14	;
19	13	14	;
19	14	14	;
19	15	14	;
19	16	14	;
19	17	14	;
10	15	15	;
11	14	15	;
11	15	15	;
12	13	15	;
12	14	15	;
12	15	15	;
13	12	15	;
13	13	15	;
13	14	15	;
13	15	15	;
14	14	15	;
15	14	15	;
16	14	15	;
17	14	15	;
14	15	15	;
15	15	15	;
16	15	15	;
17	15	15	;
18	12	15	;
18	13	15	;
18	14	15	;
18	15	15	;
18	16	15	;
18	17	15	;
19	12	15	;
19	13	15	;
19	14	15	;
19	15	15	;
19	16	15	;
19	17	15	;
10	15	16	;
11	14	16	;
11	15	16	;
12	13	16	;
12	14	16	;
12	15	16	;
13	12	16	;
13	13	16	;
13	14	16	;
13	15	16	;
14	14	16	;
15	14	16	;
16	14	16	;
17	14	16	;
14	15	16	;
15	15	16	;
16	15	16	;
17	15	16	;
18	12	16	;
18	13	16	;
18	14	16	;
18	15	16	;
18	16	16	;
18	17	16	;
19	12	16	;
19	13	16	;
19	14	16	;
19	15	16	;
19	16	16	;
19	17	16	;
];
one_35_idxs=one(:,1)*dimlen*dimlen+one(:,2)*dimlen+one(:,3)+1;
e_fig=[5 5 5  ;
25  25  25  ;
12	14	14	;
12	15	14	;
12	16	14	;
13	13	14	;
13	14	14	;
13	16	14	;
13	17	14	;
14	12	14	;
14	13	14	;
14	17	14	;
14	18	14	;
15	12	14	;
15	13	14	;
15	14	14	;
15	15	14	;
15	16	14	;
15	17	14	;
15	18	14	;
16	12	14	;
16	13	14	;
17	12	14	;
17	13	14	;
17	17	14	;
17	18	14	;
18	13	14	;
18	14	14	;
18	16	14	;
18	17	14	;
19	14	14	;
19	15	14	;
19	16	14	;
12	14	15	;
12	15	15	;
12	16	15	;
13	13	15	;
13	14	15	;
13	16	15	;
13	17	15	;
14	12	15	;
14	13	15	;
14	17	15	;
14	18	15	;
15	12	15	;
15	13	15	;
15	14	15	;
15	15	15	;
15	16	15	;
15	17	15	;
15	18	15	;
16	12	15	;
16	13	15	;
17	12	15	;
17	13	15	;
17	17	15	;
17	18	15	;
18	13	15	;
18	14	15	;
18	16	15	;
18	17	15	;
19	14	15	;
19	15	15	;
19	16	15	;
12	14	16	;
12	15	16	;
12	16	16	;
13	13	16	;
13	14	16	;
13	16	16	;
13	17	16	;
14	12	16	;
14	13	16	;
14	17	16	;
14	18	16	;
15	12	16	;
15	13	16	;
15	14	16	;
15	15	16	;
15	16	16	;
15	17	16	;
15	18	16	;
16	12	16	;
16	13	16	;
17	12	16	;
17	13	16	;
17	17	16	;
17	18	16	;
18	13	16	;
18	14	16	;
18	16	16	;
18	17	16	;
19	14	16	;
19	15	16	;
19	16	16	;
];
e_3d=zeros(dimlen,dimlen,dimlen);
for idx = 1:size(e_fig,1)
	loc=e_fig(idx,:);
	e_3d(loc(1),loc(2),loc(3))=1;
end
e_idxs=find(e_3d);
%
glass=[12	10	19	;
12	11	19	;
12	12	19	;
12	18	19	;
12	19	19	;
12	20	19	;
13	9	19	;
13	10	19	;
13	11	19	;
13	12	19	;
13	13	19	;
13	17	19	;
13	18	19	;
13	19	19	;
13	20	19	;
13	21	19	;
14	8	19	;
14	9	19	;
14	10	19	;
14	11	19	;
14	12	19	;
14	13	19	;
14	14	19	;
14	16	19	;
14	17	19	;
14	18	19	;
14	19	19	;
14	20	19	;
14	21	19	;
14	22	19	;
15	8	19	;
15	9	19	;
15	10	19	;
15	11	19	;
15	12	19	;
15	13	19	;
15	14	19	;
15	16	19	;
15	17	19	;
15	18	19	;
15	19	19	;
15	20	19	;
15	21	19	;
15	22	19	;
16	8	19	;
16	9	19	;
16	10	19	;
16	11	19	;
16	12	19	;
16	13	19	;
16	17	19	;
16	18	19	;
16	19	19	;
16	20	19	;
16	21	19	;
16	22	19	;
17	9	19	;
17	10	19	;
17	11	19	;
17	12	19	;
17	13	19	;
17	17	19	;
17	18	19	;
17	19	19	;
17	20	19	;
17	21	19	;
];
glass_mid=[16	14	19	;
17	15	19	;
16	16	19	;
];
glass_stick=[16	7	13	;
16	23	13	;
16	7	14	;
16	23	14	;
16	7	15	;
16	23	15	;
15	7	16	;
15	23	16	;
15	7	17	;
15	23	17	;
15	7	18	;
15	23	18	;
15	7	19	;
15	23	19	;
];
glass_feet=[15	7	12	;
15	23	12	;
14	7	11	;
14	23	11	;
];
glass_nose=[16	14	18	;
15	14	18	;
16	16	18	;
15	16	18	;
];
glass_3d=zeros(dimlen,dimlen,dimlen);
for idx = 1:size(glass,1)
	loc=glass(idx,:);
	glass_3d(loc(1),loc(2),loc(3))=45;
end
glass_3d = permute(glass_3d,[2 3 1]);
glass_3d = flip(glass_3d,2);
glass_3d(1,1,1)=1;
glass_3d(30,30,30)=90;

glass_mid_3d=zeros(dimlen,dimlen,dimlen);
for idx = 1:size(glass_mid,1)
	loc=glass_mid(idx,:);
	glass_mid_3d(loc(1),loc(2),loc(3))=15;
end
glass_mid_3d = permute(glass_mid_3d,[2 3 1]);
glass_mid_3d = flip(glass_mid_3d,2);

glass_stick_3d=zeros(dimlen,dimlen,dimlen);
for idx = 1:size(glass_stick,1)
	loc=glass_stick(idx,:);
	glass_stick_3d(loc(1),loc(2),loc(3))=15;
end
glass_stick_3d = permute(glass_stick_3d,[2 3 1]);
glass_stick_3d = flip(glass_stick_3d,2);

glass_feet_3d=zeros(dimlen,dimlen,dimlen);
for idx = 1:size(glass_feet,1)
	loc=glass_feet(idx,:);
	glass_feet_3d(loc(1),loc(2),loc(3))=45;
end
glass_feet_3d = permute(glass_feet_3d,[2 3 1]);
glass_feet_3d = flip(glass_feet_3d,2);

glass_nose_3d=zeros(dimlen,dimlen,dimlen);
for idx = 1:size(glass_nose,1)
	loc=glass_nose(idx,:);
	glass_nose_3d(loc(1),loc(2),loc(3))=15;
end
glass_nose_3d = permute(glass_nose_3d,[2 3 1]);
glass_nose_3d = flip(glass_nose_3d,2);

glass_idxs=find(glass_3d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
hold on
vals=nonzeros(glass_3d);
[xB,yB,zB]=ind2sub(size(glass_3d),find(glass_3d));
h1 = scatter3(xB,yB,zB,400*ones(length(xB),1),vals,'o','filled','MarkerEdgeColor','w');
vals=nonzeros(glass_mid_3d);
[xB,yB,zB]=ind2sub(size(glass_mid_3d),find(glass_mid_3d));
h2 = scatter3(xB,yB,zB,130*ones(length(xB),1),vals,'d','filled','MarkerEdgeColor','y');
vals=nonzeros(glass_stick_3d);
[xB,yB,zB]=ind2sub(size(glass_stick_3d),find(glass_stick_3d));
h3 = scatter3(xB,yB,zB,230*ones(length(xB),1),vals,'*','filled','MarkerEdgeColor','y');
vals=nonzeros(glass_feet_3d);
[xB,yB,zB]=ind2sub(size(glass_feet_3d),find(glass_feet_3d));
h4 = scatter3(xB,yB,zB,500*ones(length(xB),1),vals,'<','filled','MarkerEdgeColor','y');
vals=nonzeros(glass_nose_3d);
[xB,yB,zB]=ind2sub(size(glass_nose_3d),find(glass_nose_3d));
h5 = scatter3(xB,yB,zB,230*ones(length(xB),1),vals,'s','filled','MarkerEdgeColor','y');

alpha1 = 0.3;
alpha2 = 0.8;
set(h1, 'MarkerEdgeAlpha', alpha1, 'MarkerFaceAlpha', alpha1)
set(h4, 'MarkerEdgeAlpha', alpha2, 'MarkerFaceAlpha', alpha2)
xlabel('x');
ylabel('y');
zlabel('z');
colormap hsv
colorbar
grid on
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % k=4;
% itsPer_kMax=1;
% tic
% x=zeros(n,1); 
% % locs=randperm(n);
% % locs=gitar_35_idxs;
% locs=one_35_idxs;
% % x(locs(1:k))=(1+rand(k,1)).*(-1).^randi(2,[k,1]);% signal vector
% % x(locs(1:k))=1;% signal vector
% x(locs)=1:size(locs,1);% signal vector
% x_3d=reshape(x,dimlen,dimlen,dimlen);
% % c=abs(fft2(reshape(x,sqrt(n),sqrt(n)))).^2;% ideal measurements
% c=abs(fftn(reshape(x,dimlen,dimlen,dimlen))).^2;% ideal measurements
% c=awgn(c(:),snr,'measured'); % noisy measurements
% measurementSet=1:m; % Set of measurements (rows in M) used. 1:m uses all measurements 
% fMin=inf;
% xBest=zeros(size(x));
% Tind=0;
% verbose=1;
% while Tind<=maxT      
%     [fVal,x_n,Tind]=GreedysparseRec_3d(c,k,measurementSet,n,Tind,maxT,verbose);
%     if fVal<fMin
%         fMin=fVal;
%         xBest=x_n;
%         % xBestFlipped=bestMatch2D(xBest,x);
%         if fMin<norm(c)/1e5 
%                 break;
%         end
%     end    
% end
% xBestFlipped=bestMatch3D(x_n,x);
% xBestFlipped=xBestFlipped(:);
% % err=1-abs(x'*xBestFlipped(:)/(norm(x)*norm(xBestFlipped(:))));
% % fprintf('Error =%d\n',err);
% % xBestFlipped=xBestFlipped(:);
% t=toc;

% locsBest=find(xBestFlipped);
% vals=nonzeros(xBestFlipped);
% xB3d=zeros(dimlen,dimlen,dimlen);
% xB3d(locsBest)=vals;
% [xB,yB,zB]=ind2sub(size(xB3d),find(xB3d));

% hold on
% % h = scatter3(one(:,1),one(:,2),one(:,3),500*ones(90+2,1),92:-1:1,'s','filled','MarkerEdgeColor','k');
% % alpha = 0.5;
% % set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
% % colorbar
% h = scatter3(xB,yB,zB,500*ones(90+2,1),vals,'s','filled','MarkerEdgeColor','k');
% alpha = 0.5;
% set(h, 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha)
% colorbar
% hold off

% % figure;subplot(2,2,[3:4]);stem(1:n,real(x),'--o');hold all;subplot(2,2,[3:4]);plot(1:n,real(xBestFlipped),'rx');xlabel('index');ylabel('value');xlim([0 n]);title('True/Recovered (column stacked)');legend('True signal','Recovered');
% % subplot(2,2,1);imagesc(reshape(c/max(c),sqrt(n),sqrt(n)));colorbar;title('Measurements (Fourier magnitude)');colormap(bone);subplot(2,2,2);imagesc(reshape(x,sqrt(n),sqrt(n)));colormap(bone);title('True signal');colorbar;


% % [5 4 4 1.0;
% % 3 3 3 1.0;
% % 4 4 7 1.0;
% %  3 4 12 1.0;
% %  4 4 9 1.0;
% %  5 4 0 1.0;
% %  2 4 0 1.0;
% %  6 3 2 1.0;
% %  3 4 2 1.0;
% %  2 3 0 1.0;
% %  2 4 5 1.0;
% %  4 3 5 1.0;
% %  4 4 14 1.0;
% %  4 3 1 1.0;
% %  1 4 2 1.0;
% %  3 4 4 1.0;
% %  5 4 2 1.0;
% %  4 4 11 1.0;
% %  5 4 6 1.0;
% %  3 4 0 1.0;
% %  5 3 3 1.0;
% %  4 4 5 1.0;
% %  6 4 1 1.0;
% %  2 4 3 1.0;
% %  3 3 6 1.0;
% %  4 4 8 1.0;
% %  4 3 0 1.0;
% %  3 3 4 1.0;
% %  2 3 1 1.0;
% %  6 3 5 1.0;
% %  6 4 3 1.0;
% %  2 4 1 1.0;
% %  5 3 1 1.0;
% %  3 3 2 1.0;
% %  4 3 6 1.0;
% %  4 4 3 1.0;
% %  4 4 13 1.0;
% %  5 4 5 1.0;
% %  3 4 6 1.0;
% %  5 3 4 1.0;
% %  6 4 2 1.0;
% %  5 3 2 1.0;
% %  4 3 4 1.0;
% %  3 3 0 1.0;
% %  4 4 12 1.0;
% %  4 3 2 1.0;
% %  2 3 2 1.0;
% % 3 4 7 1.0;
% %  1 4 1 1.0;
% %  3 4 14 1.0;
% % 3 3 5 1.0;
% %  4 4 0 1.0;
% %  4 4 1 1.0;
% %  3 4 3 1.0;
% %  4 4 2 1.0;
% %  5 3 0 1.0;
% %  2 3 4 1.0;
% %  5 4 1 1.0;
% %  6 3 3 1.0;
% %  3 4 5 1.0;
% %  4 4 6 1.0;
% %  4 4 10 1.0;
% %  5 4 3 1.0;
% %  2 3 3 1.0;
% %  2 4 4 1.0;
% %  6 3 1 1.0;
% %  5 3 5 1.0;
% %  3 4 1 1.0;
% %  2 3 6 1.0;
% %  6 4 0 1.0;
% %  4 3 3 1.0;
% %  2 3 5 1.0;
% %  2 4 2 1.0;
% %  4 4 4 1.0;
% % 3 3 1 1.0;
% %  2 4 6 1.0;
% %  3 4 8 1.0;
% % 6 3 0 1.0;
% %  5 3 6 1.0;
% %  1 3 1 1.0;
% %  3 4 10 1.0;
% %  3 4 13 1.0;
% %  6 4 5 1.0;
% %  6 4 6 1.0;
% %  6 4 4 1.0;
% %  1 3 2 1.0;
% %  4 4 15 1.0;
% %  6 3 4 1.0;
% %  3 4 9 1.0;
% %  1 4 0 1.0;
% % 3 4 11 1.0;
% %  1 3 0 1.0];