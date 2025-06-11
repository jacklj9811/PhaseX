function xBest=bestMatch3D(x1,x2)
% bestMatch finds the best match between x1 that matches x2 the best in
% terms of the following 3D DFT ambiguities: circular shift, global phase, flipping
[n,m]=size(x1);
dimlen=round(n^(1/3));
x1=reshape(x1,dimlen,dimlen,dimlen);
x2=reshape(x2,dimlen,dimlen,dimlen);
minErr=inf;
phaseVec=linspace(0,pi,2);
for kk=1:dimlen
    for jj=1:dimlen
        for qq=1:dimlen
            for phaseInd=phaseVec
                for isflip=0:1
                    if (isflip)
                        x1shift=flip(flip(flip(circshift(x1,[kk jj qq]).*exp(1i*phaseInd),1),2),3);
                    else
                        x1shift=circshift(x1,[kk jj qq]).*exp(1i*phaseInd);
                    end
                    dis = x2-x1shift;
                    err=norm(dis(:));
                    if err<minErr
                        xBest=x1shift;
                        minErr=err;
                    end
                end
            end
        end
    end
end

