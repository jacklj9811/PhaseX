function [fMin,x_n,Tind] = GreedysparseRec_3d(c,k,measurementSet,n,Tind,maxT,verbose)
%GreedysparseRec Finds a k sparse solution x to c=|Fx|^2, for 2D signals
% verbose=1;
p=randperm(n);
supp=p(1:k);
iterations=1000; % GN iterations
c=c(measurementSet);
% w=(1+(rand(length(measurementSet),1)<0.5))'; % Use random weights
w=(1+(rand(length(measurementSet),1)<2))'; % Use random weights
% [x_k,iterationEr]=GN_3d(supp,c,n,randn(k,1)+1i*rand(k,1),iterations,w); % Initial guess - Gauss Newton
[x_k,iterationEr]=GN_3d(supp,c,2*n,randn(k,1)+1i*rand(k,1),iterations,w); % Initial guess - Gauss Newton
fMin=WG_cost_3d(c,x_k,w);
it=0;
while(1)
    it=it+1;
    %% Main iteration
    [junk,idx]=sort(abs(x_k(supp)));
    supp=supp(idx); % Sorting supp from min(abs(x_k)) to max
    fGrad=WGf_QU_G_Gradient_3d(c,x_k);
    offSupp=setdiff(1:n,supp);
    [junk,idx]=sort(-abs(fGrad(offSupp)));
    offSupp=offSupp(idx);
    pSupp=1:length(supp);
    pOffSupp=1:length(offSupp);
    improved=0;
        for iInd=1:length(supp)
            i=supp(pSupp(iInd)); %Index to remove
            for jInd=1:min(1,length(pOffSupp))
                j=offSupp(pOffSupp(jInd)); % Index to insert
                %% Check replacement
                suppTemp=supp;
                suppTemp((suppTemp==i))=j;
                Tind=Tind+1;
                %% Solve GN with given support
                % xTemp=GN_3d(suppTemp,c,n,x_k(suppTemp),iterations,w);
                xTemp=GN_3d(suppTemp,c,n*2,x_k(suppTemp),iterations,w);
                fTemp=WG_cost_3d(c,xTemp,w);
                if fTemp<fMin
                    if (verbose)
                        fprintf('it: %d, T: %d, Replaced %d with %d   f= %3.3f\n',it, Tind,i,j,fTemp);
                    end  
                    x_k=xTemp;
                    x_n=x_k;
                    supp=suppTemp;
                    improved=1;
                    % localopttest = (fMin-fTemp)/fMin; 
                    fMin=fTemp;
                    if fTemp<1e-3
                        if (verbose) fprintf('******************************************Success!, iteration=%d\n',it);end
                        return;
                    % else
                    %     % new
                    %     if 1
                    %         if localopttest < 1e-3
                    %             improved=0;
                    %             fprintf('But the improvement is too small. May stuck into local optimal. - trying new initial guess\n');
                    %         end
                    %         x_n=x_k;
                    %         return  
                    %     end
                    %     % new
                    end

                    break;
                end
            end
            if (improved) 
                break;
            end
        end
        if (~improved || Tind>maxT)
            if (verbose)
                fprintf('no possible improvement - trying new initial guess\n');
            end
            x_n=x_k;
            return                
        end
end
x_n=x_k;
end

