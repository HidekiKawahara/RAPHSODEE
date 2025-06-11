function y = whiteDetNoise(n,fs)
%%
yF = zeros(n,1);
fx = (0:n-1)'/n*fs;
tmp = fx(fx<=fs/2)*0+1;
tmp = tmp.*exp(1i*2*pi*rand(length(tmp),1));
tmp(1) = 1;
yF(1:length(tmp)) = tmp;
yF(end:-1:end-length(tmp)+2) = conj(tmp(2:end));
if rem(n,2) == 0
yF(n/2+1) = 1;
end
y = real(ifft(yF));
end