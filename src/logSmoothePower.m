function smootheSpec = logSmoothePower(fx,flx,spec)
cumPw = cumsum(abs(spec).^2*fx(2));
flxH = flx*2^(1/12);
flxL = flx*2^(-1/12);
bwList = (flxH - flxL);
smootheSpec = (interp1(fx,cumPw,flxH,"linear","extrap") ...
    - interp1(fx,cumPw,flxL,"linear","extrap"))./bwList(:);
end