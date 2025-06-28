function [y, output] = samplingRateConvGiantFFT(x, fsIn, fsOut, options)
%samplingRateConvGiantFFT Resample a signal using a 'Giant FFT' method.
%
%   This function implements aliasing-free sample rate conversion (SRC) by
%   performing a single FFT over a zero-padded version of the entire input
%   signal. The spectrum is then truncated or extended to the target sample
%   rate, and an IFFT reconstructs the time-domain signal.
%
%   This method is functionally equivalent to the one described in:
%   Välimäki, V., & Bilbao, S. (2023). "Giant FFTs for Sample-Rate
%   Conversion." Journal of the Audio Engineering Society, 71(3), 88–99.
%   DOI: https://doi.org/10.17743/jaes.2022.0061
%
%  Refatored by Gemini Advanced and edited by Hideki Kawahara
%  28 June, 2025
%
%   SYNTAX:
%   y = samplingRateConvGiantFFT(x, fsIn, fsOut)
%   y = samplingRateConvGiantFFT(x, fsIn, fsOut, transitionWidthHz=2000)
%   [y, output] = samplingRateConvGiantFFT(...)
%
%   INPUTS:
%   x                (double matrix) : Input signal matrix, where each
%                                      column is a channel.
%   fsIn             (double scalar) : Sampling frequency of the input
%                                      signal in Hz.
%   fsOut            (double scalar) : Desired sampling frequency of the
%                                      output signal in Hz.
%
%   OPTIONS:
%   transitionWidthHz (double scalar): Transition band width in Hz for the
%                                      low-pass filter applied before
%                                      downsampling or upsampling.
%                                      Default is 2000 Hz.
%
%   OUTPUTS:
%   y                (double matrix) : Resampled output signal.
%   output           (struct)        : Optional struct with detailed data
%                                      for debugging and analysis.
%
%   Copyright 2025 Hideki Kawahara
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

arguments
    x double
    fsIn (1,1) double {mustBeNumeric, mustBePositive}
    fsOut (1,1) double {mustBeNumeric, mustBePositive}
    options.transitionWidthHz (1,1) double {mustBeNumeric, mustBePositive} = 2000
end

startTic = tic;
[numSamplesIn, numChannels] = size(x);

%% 1. Calculate Optimal FFT Length for a Common Frequency Grid
% To ensure both fsIn and fsOut fall on the same discrete frequency grid,
% we calculate a common base length and pad the signal to a multiple of it.
% This prevents spectral leakage from the resampling process itself.
baseFftLength = fsIn * fsOut / (gcd(fsIn, fsOut)^2);
fftLengthIn = baseFftLength * ceil(numSamplesIn / baseFftLength);

% Zero-pad the input signal to the calculated FFT length.
paddedInput = [x; zeros(fftLengthIn - numSamplesIn, numChannels)];

%% 2. Perform Forward FFT
% The FFT is normalized by the input length. The IFFT will be scaled by
% the output length to maintain the signal's amplitude.
inputSpectrum = fft(paddedInput) / fftLengthIn;
inputFrequencies = (0:fftLengthIn-1)' / fftLengthIn * fsIn;

%% 3. Design Low-Pass Anti-Aliasing/Imaging Filter
% A low-pass filter is applied in the frequency domain to prevent aliasing
% (in downsampling) or imaging (in upsampling). A half-cosine window
% creates a smooth transition band.
nyquistLimit = min(fsIn, fsOut);
targetNyquist = nyquistLimit / 2;

% Select frequencies up to the target Nyquist
% manually fixed start
[minv,bestID]=min(abs(inputFrequencies-targetNyquist));
if minv>10^(-10) % fragile: fails when signal is longer then 28 hous
    positiveFreqsIdx = inputFrequencies <= targetNyquist;
else
    positiveFreqsIdx = 1:bestID;
end
% manually fixed end
positiveFreqs = inputFrequencies(positiveFreqsIdx);

% Design the shaper (low-pass filter)
shaper = ones(length(positiveFreqs), 1);
transitionStartFreq = targetNyquist - options.transitionWidthHz;

% Check if the transition band is valid
if transitionStartFreq > 0
    is_in_transition = positiveFreqs > transitionStartFreq;

    % Calculate normalized frequencies (0 to 1) within the transition band
    normalizedTransitionFreqs = (positiveFreqs(is_in_transition) - transitionStartFreq) / options.transitionWidthHz;

    % Apply a half-cosine window shape to the transition band
    shaper(is_in_transition) = 0.5 + 0.5 * cos(normalizedTransitionFreqs * pi);
end

% Apply the filter to the positive frequencies of the input spectrum
outputSpectrumHalf = inputSpectrum(positiveFreqsIdx, :);
outputSpectrumHalf = outputSpectrumHalf .* shaper;

%% 4. Construct the Output Spectrum (Upsampling vs. Downsampling)
fftLengthOut = round(fftLengthIn * fsOut / fsIn); % Use round for precision
outputSpectrum = zeros(fftLengthOut, numChannels, 'like', 1j); % Pre-allocate

if fsIn > fsOut % Downsampling
    % The new spectrum is a truncated version of the filtered input spectrum.
    % We must reconstruct the full spectrum ensuring Hermitian symmetry for a
    % real-valued output signal.

    numPositiveBins = size(outputSpectrumHalf, 1);
    outputSpectrum(1:numPositiveBins, :) = outputSpectrumHalf;

    if rem(fftLengthOut, 2) == 0 % Even FFT length
        % The Nyquist bin must be real.
        outputSpectrum(numPositiveBins, :) = real(outputSpectrum(numPositiveBins, :));
        % Create conjugate symmetric part
        outputSpectrum(numPositiveBins+1:end, :) = conj(outputSpectrum(numPositiveBins-1:-1:2, :));
    else % Odd FFT length
        % Create conjugate symmetric part
        outputSpectrum(numPositiveBins+1:end, :) = conj(outputSpectrum(numPositiveBins:-1:2, :));
    end
else % Upsampling
    % The new spectrum is a zero-padded version of the input spectrum.
    numPositiveBins = size(outputSpectrumHalf, 1);

    % Place the positive frequencies
    outputSpectrum(1:numPositiveBins, :) = outputSpectrumHalf;

    % Place the conjugate symmetric negative frequencies
    if rem(fftLengthOut, 2) == 0 % Even FFT length
        outputSpectrum(end:-1:end-numPositiveBins+2, :) = conj(outputSpectrumHalf(2:numPositiveBins, :));
    else % Odd FFT length
        outputSpectrum(end:-1:end-numPositiveBins+2, :) = conj(outputSpectrumHalf(2:end,:));
    end
end

%% 5. Perform Inverse FFT and Truncate to Final Length
% Perform IFFT and scale by the output length to get the final amplitude.
paddedOutput = ifft(outputSpectrum) * fftLengthOut;

% Truncate the zero-padded signal to the correct output length
numSamplesOut = round(numSamplesIn * fsOut / fsIn);
y = paddedOutput(1:numSamplesOut, :);

% Ensure output is real, removing any minor imaginary part from numerical error
if isreal(x)
    y = real(y);
end

%% 6. Populate Output Structure for Debugging (if requested)
if nargout > 1
    output.baseFftLength = baseFftLength;
    output.fftLengthIn = fftLengthIn;
    output.fftLengthOut = fftLengthOut;
    output.originalLength = numSamplesIn;
    output.outputLength = numSamplesOut;
    output.paddedInput = paddedInput;
    output.inputFrequencies = inputFrequencies;
    output.outputSpectrum = outputSpectrum;
    output.filterShaper = shaper;
    output.positiveFrequencies = positiveFreqs;
    output.transitionWidthHz = options.transitionWidthHz;
    output.elapsedTime = toc(startTic);
end

end
