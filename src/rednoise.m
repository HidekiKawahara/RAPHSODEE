function y = rednoise(n,fs)
%%
yF = zeros(n,1);
fx = (0:n-1)'/n*fs;
tmp = 1.0 ./ max(40,fx(fx<=fs/2));
tmp = tmp.*exp(1i*2*pi*rand(length(tmp),1));
yF(1:length(tmp)) = tmp;
yF(end:-1:end-length(tmp)+2) = conj(tmp(2:end));
y = real(ifft(yF));
end