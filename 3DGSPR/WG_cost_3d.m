function out=WG_cost_3d(c,x,w)
% Returns the squared 2 norm of the consistency error 
[m,n]=size(x);
dimlen=round(m^(1/3));
x=reshape(x,dimlen,dimlen,dimlen);
y=abs(fftn(x)).^2;
% out=sqrt(sum(w'.*((y(:)-c).^2)));
out=(sum(w'.*((y(:)-c).^2)))^(1/3);
end


