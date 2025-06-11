function xBest=bestMatch2D(x1,x2)
% bestMatch finds the best match between x1 that matches x2 the best in
% terms of the following 2D DFT ambiguities: circular shift, global phase, flipping
[n,m]=size(x1);
x1=reshape(x1,sqrt(n),sqrt(n));
x2=reshape(x2,sqrt(n),sqrt(n));
minErr=inf;
phaseVec=linspace(0,pi,2);
for kk=1:sqrt(n)
    for jj=1:sqrt(n)
        for phaseInd=phaseVec
            for flip=0:1
                if (flip)
                    x1shift=fliplr(flipud(circshift(x1,[kk jj]).*exp(1i*phaseInd)));
                else
                    x1shift=circshift(x1,[kk jj]).*exp(1i*phaseInd);
                end
                err=norm(x2-x1shift,'fro');
                if err<minErr
                    xBest=x1shift;
                    minErr=err;
                end
            end
        end
    end
end

