function out=WGf_QU_G_Gradient_3d(c,x)
% Returns gradient of cost function around current estimate x
[m,n]=size(x);
dimlen=round(m^(1/3));
z=fftn(reshape(x,dimlen,dimlen,dimlen));
c=reshape(c,dimlen,dimlen,dimlen);
out=ifftn((abs(z).^2-c).*z);
out=out(:);
end

