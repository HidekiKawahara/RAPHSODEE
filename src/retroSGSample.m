%% Pink noise simulation using Queen Mary University RIR data
% The original test signal is periodic pink noise made from sum of sinusoid
%
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

roomList = {'greathall','octagon','classroom'};
roomID = 2;
roomName = roomList{roomID};
rowList = {'F','B'};% Front and Back
rowID = 2;
rowName = rowList{rowID};

% Please replace this based on you download
baseDir =  ['/Volumes/HD-CD-A/RIRdatabase/QMLirDB/' roomName 'Omni/Omni/'];
switch rowName
    case rowList{1}
        switch roomName
            case roomList{1}
                fileRIR = 'x06y00.wav';
            case roomList{2}
                fileRIR = 'x06y00.wav';
            case roomList{3}
                fileRIR = '30x00y.wav';
        end
    case rowList{2}
        switch roomName
            case roomList{1}
                fileRIR = 'x06y12.wav';
            case roomList{2}
                fileRIR = 'x06y12.wav';
            case roomList{3}
                fileRIR = '30x45y.wav';
        end
end
[hIn,fsRIR] = audioread([baseDir fileRIR]);
%% set simulation fs to 44100
fs = 44100;
h = samplingRateConvByDFTwin(hIn,fsRIR,fs);
period = length(h);
nRepeat = 20;
nData = nRepeat*period;
xOrg = pinknoise(nData);
shaper = zeros(nData,1);
shaper(period+(1:(nRepeat-3)*period)) = 1;
shaper = fftfilt(hanning(period/2),shaper);
shaper = shaper/max(shaper);
xOrg = xOrg.*shaper;
%%
yOrg = fftfilt(h,xOrg);
thLevel = -5;
fHigh = 9000;
[xSg,prmSG] = signalSafeguardwithGiantFFTSRC(xOrg,fs,fs,thLevel,fHigh);
ySg = fftfilt(h,xSg); % system output with the safeguarded signal input
% Following adds trailing zeros to make the input signal length
% equar to the safeguarded signal
xOrgExt = [xOrg;zeros(10*period,1)];
xOrgExt = xOrgExt(1:length(xSg));
yOrgExt = fftfilt(h,xOrgExt); % system output with the original signal input

hSg = ifft(fft(ySg)./fft(xSg));
hOrg = ifft(fft(yOrgExt)./fft(xOrgExt));
hOrgSg = ifft(fft(yOrgExt)./fft(xSg));
%% measurement simulation generata "addedNoise" to simulate
% the observation noise
% Then, evaluate estimation error

baseNoise = pinknoise(length(xSg)); % you can try other noises
baseNoise = baseNoise/std(baseNoise)*std(xSg);
snr = 10; % SNR for observation noise.
mag = 10^(-snr/20);
addedNoise = baseNoise*mag;

hSgNoise = ifft(fft(ySg+addedNoise)./fft(xSg));
hOrgNoise = ifft(fft(yOrgExt+addedNoise)./fft(xOrgExt));
hOrgSgNoise = ifft(fft(yOrgExt+addedNoise)./fft(xSg));

errorOriginal = 20*log10(std(hOrgNoise(1:fs)-h(1:fs))/std(h(1:fs)));
errorRetroSG = 20*log10(std(hOrgSgNoise(1:fs)-h(1:fs))/std(h(1:fs)));
errorSG = 20*log10(std(hSgNoise(1:fs)-h(1:fs))/std(h(1:fs)));

%% Display the results
tx = (1:length(xSg))/fs;
figure;plot(tx,20*log10(abs(hOrgNoise)));grid on
hold on;
plot(tx,20*log10(abs(hOrgSgNoise)));
plot(tx,20*log10(abs(hSgNoise)));
plot(tx(1:length(h)),20*log10(abs(h)));
axis([0 4 -100 0])
xlabel('time (s)')
ylabel('level (dB)');
title("decay of impulse response "+ " SNR:" + num2str(10) + " (dB)")
legend("Original","retro SG","SG input","Truth", ...
    "Location","northeast");
text(1,-13, "original error: "+num2str(errorOriginal) + " dB")
text(1,-23, "retro SG error: "+num2str(errorRetroSG) + " dB")
text(1,-33, "SG error: "+num2str(errorSG) + " dB")
%
fx = (0:fs-1)';
figure;semilogx(fx,20*log10(abs(fft(hOrgNoise(1:fs)))));grid on;
hold on;
semilogx(fx,20*log10(abs(fft(hOrgSgNoise(1:fs)))));
semilogx(fx,20*log10(abs(fft(hSgNoise(1:fs)))));
semilogx(fx,20*log10(abs(fft(h(1:fs)))));
axis([10 fs/2 -70 20])
legend("Original","SG input","retro SG","Truth", ...
    "Location","southwest");
text(100,-43, "original error: "+num2str(errorOriginal) + " dB")
text(100,-53, "retro SG error: "+num2str(errorRetroSG) + " dB")
text(100,-63, "SG error: "+num2str(errorSG) + " dB")
