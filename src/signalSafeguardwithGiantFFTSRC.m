function [y, output] ...
    = signalSafeguardwithGiantFFTSRC(xOriginal,fsIn,fsOut,thLevel,varargin)
% signal safeguarding with giant FFTs sampling-rate conversion
%  This sampling rate conversion does not introduce aliasing effects.
%
% Pattern 1
%  y = signalSafeguardwithGiantFFTSRC(xOriginal,fsIn,fsOut,thLevel)
%  Arguments
%      xOriginal : input data vector with sampling frequency fsIn
%      fsIn      : sampling frequency of input signal (Hz)
%      fsOut     : sampling frequency of output signal (Hz)
%      thLevel   : threshould level from average spectrum (dB)
%  Output
%      y    : converted output with sampling frquency fsOut
%
% Pattern 2
%  [y, output] = signalSafeguardwithGiantFFTSRC(xOriginal,fsIn,fsOut,thLevel)
% Output
%      output  : structure with detailed debug data and shows spectrum figure
%
% Pattern 3
%  [y, output] = signalSafeguardwithGiantFFTSRC(xOriginal,fsIn,fsOut,thLevel,fHigh)
% Arguments (additional)
%      fHigh  : high frequency end of safeguarding (Hz)
%               default Nyquist frequency - 3000
% Output
%      output  : structure with detailed debug data and shows spectrum figure
%
% Pattern 4
%  [y, output] = signalSafeguardwithGiantFFTSRC(xOriginal,fsIn,fsOut,thLevel,fHigh,dispOn)
% Arguments (additional)
%      dispOn  : Figure display switch, 0:off (default), 1:on
% Output
%      output  : structure with detailed debug data and shows spectrum figure
%
% Pattern 5
%  [y, output] = signalSafeguardwithGiantFFTSRC(xOriginal,fsIn,fsOut,thLevel,fHigh,dispOn,Optmz)
% Arguments (additional)
%      Optmz   : Buffer length optimization, 0:off (default), 1:on
%                This optimization changes the safeguarded signal length
%                to optimize performance.
% Output
%      output  : structure with detailed debug data and shows spectrum figure


%Copyright 2025 Hideki Kawahara
%
%Licensed under the Apache License, Version 2.0 (the "License");
%you may not use this file except in compliance with the License.
%You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS,
%WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%See the License for the specific language governing permissions and
%limitations under the License.

startTic = tic;
Optmz = 0;
if nargin == 7
    Optmz = varargin{3};
end
if fsIn == fsOut % buffer length tuning
    xIn = xOriginal;
    lData = length(xIn);
    if Optmz && rem(lData,11025)>0 && rem(lData,8000)>0
        nData = min([2^ceil(log2(lData)),2^ceil(log2(lData/3))*3, ...
            2^ceil(log2(lData/5))*5,2^ceil(log2(lData/7))*7]);
    else
        nData = lData;
    end
else
    xIn = samplingRateConvByDFTwin(xOriginal,fsIn,fsOut);
    nData = length(xIn);
end
fs = fsOut;
[~, nChannel] = size(xIn);
x = [xIn;zeros(nData-length(xIn),nChannel)];
%% giant FFT: Fourier transform of whole signal
xF = fft(x);
fx = (0:nData-1)'/nData*fs;
fxBi = fx;
fxBi(fxBi >fs/2) = fxBi(fxBi >fs/2)-fs;
xF(abs(xF)==0) = eps;
%% Smoothed spectrum using Gaussian smoother
% Attenuation at 1/6 octave from the center is 1/sqrt(2)
flxBi = log2(abs(fxBi));
flxBi(1) = flxBi(2);
fc = 20*2.0 .^(0:1/24:log2(fs/2/20))';
flc = log2(fc);
%nChannel = length(flc);
%fBank = zeros(nData,nChannel);
%avPw = zeros(nChannel,1);
%for ii = 1:nChannel
%    flb = (1/12)/sqrt(log(sqrt(2)));
%    fBank(:,ii) = exp(-(flc(ii)-flxBi).^2/flb^2);
%    avPw(ii) = sum(fBank(:,ii).*abs(xF).^2)/sum(fBank(:,ii));
%end
flb = (1/12)/sqrt(log(sqrt(2)));
fBank = exp(-(flc' - flxBi).^2 / flb^2);
avPw = sum(fBank .* abs(xF).^2, 1) ./ sum(fBank, 1);
avPw = avPw';
absFxBi = abs(fxBi);
avPwi = interp1([0;fc],10*log10(avPw([1 1:end])),absFxBi,"linear","extrap");
%% Set signal safeguarding shaper
cumPw = 10*log10(cumsum(abs(xF).^2)/sum(abs(xF).^2));
fLow = max(fx(cumPw<cumPw(end)-23));
if isempty(fLow) || fLow < 10
    fLow = 10;
end
%fLow = 70;
displayOn = 0;
if fsIn < fsOut
    fHigh = fsIn/2-3000;
else
    fHigh = fs/2-3000;
end
if nargin == 5
    fHigh = varargin{1};
elseif nargin == 6
    fHigh = varargin{1};
    displayOn = varargin{2};
elseif nargin == 7
    fHigh = varargin{1};
    displayOn = varargin{2};
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
    .* xF(abs(xF)<guardThreshold) ./abs(xF(abs(xF)<guardThreshold));
xSg = real(ifft(xFfix));

%% output setting
if displayOn
    xStdDb = 20*log10(std(xF));
    constantGuard = 10.0 .^((xStdDb+thLevel)/20);
    xFfixC = xF;
    xFfixC(abs(xF)<constantGuard) = constantGuard * ...
        xFfixC(abs(xF)<constantGuard) ./abs(xFfixC(abs(xF)<constantGuard));
    xSgConst = real(ifft(xFfixC));
    maxLevel = max(20*log10(abs(xF)));
    figure;
    varSGFighandle = gcf;
    semilogx(fx, 20*log10(abs(xF)),"LineWidth",2);grid on;
    set(gca,"LineWidth",2,"fontsize",14);
    hold on
    semilogx(fx, 20*log10(abs(fft(xSg))),"LineWidth",2);grid on;
    semilogx(fc,10*log10(avPw),"LineWidth",2);
    axis([20 fs/2 maxLevel+[-95 5]])
    legend("original","safeguard","smoothed av.","Orientation","horizontal","Location","southwest");
    xlabel("frequency (Hz)")
    ylabel("level (dB)")
    varSGhandle = gca;

    %
    figure;
    constSGFighandle = gcf;
    semilogx(fx, 20*log10(abs(xF)),"LineWidth",2);grid on;
    set(gca,"LineWidth",2,"fontsize",14);
    hold on
    semilogx(fx, 20*log10(abs(fft(xSgConst))),"LineWidth",2);grid on;
    semilogx(fc,xStdDb+0*fc,"LineWidth",2);
    %yline(xStdDb,"LineWidth",2);
    axis([20 fs/2 maxLevel+[-95 5]])
    legend("original","safeguard","average level","Orientation","horizontal","Location","southwest");
    xlabel("frequency (Hz)")
    ylabel("level (dB)")
    constSGhandle = gca;
    if nargout == 2
        output.varSGhandle = varSGhandle;
        output.constSGhandle = constSGhandle;
        output.varSGFighandle = varSGFighandle;
        output.constSGFighandle = constSGFighandle;
        output.xSgConst = xSgConst;
    end
end

if nargout == 2
    y = xSg;
    output.nargout = nargout;
    output.nargin = nargin;
    output.fLow = fLow;
    output.fHigh = fHigh;
    output.x = x;
    output.xF = xF;
    output.fx = fx;
    output.avPw = avPw;
    output.fc = fc;
    output.avPwi = 10.^(avPwi/10);
    output.whiteDetNoise = whiteDetNoise(nData,fs);
    output.elapsedTime = toc(startTic);
else
    y = xSg;
end
end