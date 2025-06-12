function x=DGN(w,S,c,n,x0,iterations,F)
%%% Damped Gauss Newton
x=zeros(2*n,1);
x(S)=x0; % S=support set, x0=z0 random vector/initial guess in article
s=0.5; % s=t0 in article
for i=1:iterations % iteration=L in article
    s=min([2*s,1]); % now s is u in article
    y=fft(x); % x->fft->y=hat(x)=Fx on page 5 in article, NOT y in article! NOT!
    % (IN ARTICLE NOT HERE)we have x'A_ix=|hat(x)_i|^2=z'B_iz where
    % (IN ARTICLE NOT HERE)B_i=U_S'A_iU_S, x=U_Sz
    % (FiUs)*(fft(x))=sim=(a+bi)*(c+di)=(a-bi)(c+di)=(ac+bd)+(ad-bc)i[ignore imaginary part for that we don't need imaginary part of z, while imaginary part of z is completely generated from the imaginary part of J]
    b=sqrt(w).*(abs(y).^2+c); % b=b in article, w=w in article, c=y in article
    B=bsxfun(@times,real(y).*sqrt(w),real(F(:,S)))... %<=>(real(y).*sqrt(w)).*real(F(:,S))
        +bsxfun(@times,imag(y).*sqrt(w),imag(F(:,S))); %<=>similarly...
    xold=x;
    fold=objectiveFun(w,c,xold);
    x=zeros(2*n,1);
    x(S)=2*B\b; % solve min Bx-b, output x then get x/2
    if rank(B)<length(S)
        pause; %  temporarily stops MATLAB? execution and waits for the user to press any key
    end
    xnew=x;
    while((objectiveFun(w,c,xold+s*(xnew-xold))>fold))% && (s>1e-5))
        s=0.5*s; % now s is (0.5)^m u in article
    end
    x=xold+s*(xnew-xold);
    if (norm(x-xold)<1e-4)
        return
    end
end



