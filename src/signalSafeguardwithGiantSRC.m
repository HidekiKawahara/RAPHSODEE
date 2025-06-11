function [y, output] = signalSafeguardwithGiantSRC(xOriginal,fsIn,fsOut,thLevel,varargin)
startTic = tic;
if fsIn == fsOut
    fs = fsIn;
    fsIn = 96000;
    fsOut = 44100;
    xIn = xOriginal;
else
    xIn = samplingRateConvByDFTwin(xOriginal,fsIn,fsOut);
    fs = fsOut;
end
baseLength = fsIn*fsOut/gcd(fsIn,fsOut)^2;
nData = baseLength*ceil(length(xIn)/baseLength);
[~, nChannel] = size(xIn);
x = [xIn;zeros(nData-length(xIn),nChannel)];
%% 全体を離散フーリエ変換します
xF = fft(x);
fx = (0:nData-1)'/nData*fs;
fxBi = fx;
fxBi(fxBi >fs/2) = fxBi(fxBi >fs/2)-fs;
%% 周波数平滑化用のガウス関数を設計します（対数波数軸）
% 1/6オクターブ離れた位置でレベルが1/sqrt(2)になるように設定する
flxBi = log2(abs(fxBi));
flxBi(1) = flxBi(2);
fc = 20*2.0 .^(0:1/24:log2(fs/2/20))';
flc = log2(fc);
nChannel = length(flc);
fBank = zeros(nData,nChannel);
avPw = zeros(nChannel,1);
for ii = 1:nChannel
    flb = (1/12)/sqrt(log(sqrt(2)));
    fBank(:,ii) = exp(-(flc(ii)-flxBi).^2/flb^2);
    avPw(ii) = sum(fBank(:,ii).*abs(xF).^2)/sum(fBank(:,ii));
end
absFxBi = abs(fxBi);
avPwi = interp1([0;fc],[-50;10*log10(avPw)],absFxBi,"linear","extrap");
%% 保護のレベルを設定します
cumPw = 10*log10(cumsum(abs(xF).^2)/sum(abs(xF).^2));
fLow = max(fx(cumPw<cumPw(end)-23));
%fLow = 70;
if fsIn < fsOut
    fHigh = fsIn/2-3000;
else
    fHigh = fs/2-3000;
end
if nargin == 5
    fHigh = varargin{1};
end
%thLevel = -10;
lowFloor = max(avPwi(fx <= fLow)+thLevel);
highFloor = max(avPwi(fxBi >= fHigh)+thLevel);
guardShaper = flxBi*0;
guardShaper(absFxBi<=fLow) = lowFloor;
guardShaper(absFxBi>=fHigh) = highFloor;
flxBiSeg = flxBi(absFxBi>fLow & absFxBi <fHigh);
guardShaper(absFxBi>fLow & absFxBi <fHigh) =  lowFloor ...
    + (highFloor-lowFloor)/(max(flxBiSeg)-flxBiSeg(1))*(flxBiSeg-flxBiSeg(1));
guardShaper = max(guardShaper,avPwi+thLevel);
%%
guardThreshold = 10.0 .^(guardShaper/20);
xFfix = xF;
xFfix(abs(xF)<guardThreshold) = guardThreshold(abs(xF)<guardThreshold) ...
    .* xFfix(abs(xF)<guardThreshold) ./abs(xFfix(abs(xF)<guardThreshold));
xSg = real(ifft(xFfix));

%% output setting
if nargout == 2
    figure;
    set(gcf,"Position",[680   632   611   246])
    semilogx(fx, 20*log10(abs(xF)),"LineWidth",2);grid on;
    set(gca,"xlim",[20 fs/2],"LineWidth",2,"fontsize",14);
    hold on
    semilogx(fx, 20*log10(abs(fft(xSg))),"LineWidth",2);grid on;
    semilogx(fc,10*log10(avPw),"LineWidth",2);
    axis([20 fs/2 -40 60])
    legend("original","safeguard","smoothed av.","Orientation","horizontal","Location","southwest");
    xlabel("frequency (Hz)")
    ylabel("level (dB)")
    y = xSg;
    output.nargout = nargout;
    output.nargin = nargin;
    output.fLow = fLow;
    output.fHign = fHigh;
    output.x = x;
    output.xF = xF;
    output.fx = fx;
    output.avPw = avPw;
    output.fc = fc;
    output.avPwi = 10.^(avPwi/10);
    output.elapsedTime = toc(startTic);
else
    y = xSg;
end
end