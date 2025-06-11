function [x,errVec]=GN_3d(S,c,n,x0,iterations,w)
    %\/
    dimlen=round(n^(1/3));
    %/\
    x=zeros(n,1);
    x(S)=x0;
    s=0.5; 
    W=dftmtx(dimlen);
    ind=0;
    %%% Building MM matrix on the fly
    for indS=S
        alpha=ceil(indS/dimlen^2);
        indS=indS-(alpha-1)*dimlen^2;
        beta=ceil(indS/dimlen);
        gamma=mod(indS,dimlen);
        if (gamma==0)
            gamma=dimlen;
        end
        ind=ind+1;
        MM(:,ind)=kron(W(:,alpha),kron(W(:,beta),W(:,gamma)));
    end
    %%% GN
    for i=1:iterations
        s=min([2*s,1]);
        z=fftn(reshape(x,dimlen,dimlen,dimlen));
        z=z(:);
        XM=w.*(z');
        B=bsxfun(@times,real(z).*sqrt(w'),real(MM))+bsxfun(@times,imag(z).*sqrt(w'),imag(MM));
        b=sqrt(w').*(c+abs(z).^2);
        xold=x;
        fold=WG_cost_3d(c,xold,w);
        x=zeros(n,1);
        x(S)=2*B\b;
        xnew=x;
        while((WG_cost_3d(c,xold+s*(xnew-xold),w)>fold))% && (s>1e-5))
            s=0.5*s;
        end
        x=xold+s*(xnew-xold);
        errVec(i)=fold;
        if (norm(x-xold)<1e-4)
            return
        end
    end
end