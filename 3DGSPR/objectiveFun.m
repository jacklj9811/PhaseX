function out=objectiveFun(w,c,x)
% Returns the squared 2 norm of the consistency error 
out=sum(w.*((abs(fft(x)).^2-c).^2)); % f(x) on page "5" in article
% c is y
end

