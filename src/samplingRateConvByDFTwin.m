function [y,output] = samplingRateConvByDFTwin(x,fsIn,fsOut,varargin)
% sampling rate converter based on one shot DFT
%  This sampling rate conversion does not introduce aliasing effects.
%
% Pattern 1
%  y = samplingRateConvByDFTwin(x,fsIn,fsOut)
%  Argument
%      x    : input data vector or matrix, with sampling frequency fsIn
%      fsIn : sampling frequency of input signal (Hz)
%      fsOut: sampling frequency of output signal (Hz)
%  Output
%      y    : converted output or matrix, with sampling frquency fsOut
%
% Pattern 2
%  [y,output] = samplingRateConvByDFTwin(x,fsIn,fsOut)
% Output
%      output  : structure with detailed debug data
%
% Pattern 3
%  [y,output] = samplingRateConvByDFTwin(x,fsIn,fsOut,transW)
% Artument (additional)
%      transW  : transition width at high-end (Hz), default is 4000
% Output
%      output  : structure with detailed debug data
%
% Although this function was developed before reading the following
% article, this is an implementation of the same idea. Highly recommended
% to read
%
% Välimäki, V., & Bilbao, S. (2023). Giant FFTs for Sample-Rate Conversion.
% Journal of the Audio Engineering Society, 71(3), 88–99.
% DOI: https://doi.org/10.17743/jaes.2022.0061

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
[nData, nChannel] = size(x);
transW = 4000;
if nargin == 4
    transW = varargin{1};
end
baseLength = fsIn*fsOut/gcd(fsIn,fsOut)^2;
l_buffer = baseLength*ceil(length(x)/baseLength);
xExt = [x;zeros(l_buffer-nData,nChannel)];
xF = fft(xExt)/(l_buffer);
fxx = (0:l_buffer-1)'/l_buffer*fsIn;
l_bufferOut = l_buffer*fsOut/fsIn;
tmpFs = min(fsIn,fsOut);
fxxTrim = fxx(fxx<= tmpFs/2);
shaper = ones(length(fxxTrim),1);
tmp = (fxxTrim(fxxTrim>tmpFs/2-transW)-(tmpFs/2-transW))/transW;
shaper(fxxTrim>tmpFs/2-transW) = 0.5+0.5*cos(tmp*pi);
xOutFHalf = xF(fxx<= tmpFs/2,:);
if fsIn > fsOut
    xOutFHalf = xOutFHalf.*shaper;
    if rem(l_bufferOut,2) == 0
        xOutFHalf(end,:) = real(xOutFHalf(end,:));
        xOutF = [xOutFHalf;conj(xOutFHalf(end-1:-1:2,:))];
    else
        xOutF = [xOutFHalf;conj(xOutFHalf(end:-1:2,:))];
    end
else
    xOutF = zeros(l_bufferOut,nChannel);
    l_half = length(xOutFHalf);
    xOutFHalf = xOutFHalf.*shaper;
    xOutF(1:l_half,:) = xOutFHalf;
    xOutF(end:-1:end-l_half+2,:) = conj(xOutFHalf(2:end,:));
end
yExt = ifft(xOutF)*(l_bufferOut);
nDataOut = round(nData*fsOut/fsIn);
y = yExt(1:nDataOut,:);
if nargout > 1
    output.baseLength = baseLength;
    output.l_buffer = l_buffer;
    output.originalLength = nData;
    output.xExt = xExt;
    output.fxx = fxx;
    output.fxxHalf = fxx(fxx<=fsOut/2);
    output.xOutF = xOutF;
    output.l_bufferOut = l_bufferOut;
    output.y = y;
    output.elapsedTime = toc(startTic);
    output.shaper = shaper;
    output.transitionWidth = transW;
    output.nargout = nargout;
    output.nargin = nargin;
end
end