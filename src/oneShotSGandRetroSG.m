%% script for ASJ AA
[fileNr, pathNr] = uigetfile('*.wav'); % select reference structured signal file
[xOrg, fs] = audioread(string(pathNr)+string(fileNr));
playerRecorder = audioPlayerRecorder("SampleRate",fs,"BitDepth","24-bit integer");
devices = getAudioDevices(playerRecorder);
disp(devices);
t = datetime(datetime,'format','yyyyMMddHHmmss');
dirName = "simulIOSGndRetroSG" + string(t);
mkdir(dirName);
thLevel = input("Threshold level (dB):");
fHighIn = input("High freq. limit (Hz):");
[y, output] = signalSafeguardwithGiantFFTSRC(xOrg,fs,fs,thLevel,fHighIn,1);
disp("SG completed");
%%
devId = input("select device ID:");
outputSG = basicPlayRecLoopMono(fs,y,devices{devId});
disp("SG i/o completed");
outputRaw = basicPlayRecLoopMono(fs,xOrg,devices{devId});
disp("Raw i/o completed");
%%
xInSG = outputSG.outputSignal(1:length(xOrg));
xOutSG = outputSG.acquiredSignal(1:length(xOrg));
xOutRaw = outputRaw.acquiredSignal(1:length(xOrg));
iRespSG = ifft(fft(xOutSG)./fft(xInSG));
iRespRawRSG = ifft(fft(xOutRaw)./fft(xInSG));
%%
nData = length(iRespSG);
figure;
set(gcf,"Position",[1 1 1074 300])
plot((1:fs)/fs,20*log10(abs(iRespSG(1:fs))),"LineWidth",2);grid on;
set(gca,"FontSize",13,"LineWidth",2)
axis([0 1 -100 0])
xlabel("time(s)")
ylabel("level (dB)")
eval(['print -dpng -r200 ' char(dirName) '/ipBySG.png']);
%%
figure;
set(gcf,"Position",[1 1 1074 300])
plot((1:fs)/fs,20*log10(abs(iRespRawRSG(1:fs))),"LineWidth",2);grid on;
hold on;plot((1:fs)/fs,20*log10(abs(iRespSG(1:fs))),"LineWidth",2);grid on;
set(gca,"FontSize",13,"LineWidth",2)
axis([0 1 -100 0])
xlabel("time(s)")
ylabel("level (dB)")
legend("Retro-SG","SG -5dB")
eval(['print -dpng -r200 ' char(dirName) '/sgAndRetroSGip.png']);
%%
[~,peakIdx] = max(abs(iRespSG));
headLTI = peakIdx - round(0.05*fs);
lIresp = round(0.2*fs);
fxSG = (0:length(xOrg)-1)'/length(xOrg)*fs;
ltiSG = 20*log10(abs(fft(iRespSG(headLTI+(1:lIresp)),length(xOrg))));
rtvSG = 20*log10(abs(fft(iRespSG(round(nData/2)+(1:lIresp)),length(xOrg))));
ltiRawSG = 20*log10(abs(fft(iRespRawRSG(headLTI+(1:lIresp)),length(xOrg))));
rtvRawSG = 20*log10(abs(fft(iRespRawRSG(round(nData/2)+(1:lIresp)),length(xOrg))));
%%
lowLimit = output.fLow;
highLimit = output.fHigh;
figure;
semilogx(fxSG,ltiSG,"LineWidth",2);grid on
hold on;semilogx(fxSG,rtvSG,"LineWidth",2);grid on;
xline(lowLimit,"LineWidth",2,"Color",'g');
xline(highLimit,"LineWidth",2,"Color",'g');
set(gca,"FontSize",13,"LineWidth",2)
axis([20 fs/2 -60 20])
xlabel("Frequency (Hz)")
ylabel("level (dB)")
legend("LTI","RTV","Location","northwest")
title(string(fileNr) + " safeguard, Thlevel:" + num2str(thLevel) + " dB")
eval(['print -dpng -r200 ' char(dirName) '/fRespSG.png']);
%%
lowLimit = output.fLow;
highLimit = output.fHigh;
figure;
semilogx(fxSG,ltiRawSG,"LineWidth",2);grid on
hold on;semilogx(fxSG,rtvRawSG,"LineWidth",2);grid on;
xline(lowLimit,"LineWidth",2,"Color",'g');
xline(highLimit,"LineWidth",2,"Color",'g');
set(gca,"FontSize",13,"LineWidth",2)
axis([20 fs/2 -60 20])
xlabel("Frequency (Hz)")
ylabel("level (dB)")
legend("LTI","RTV","Location","northwest")
title(string(fileNr) + " retro-safeguard, ThLevel: " + num2str(thLevel) + " dB")
eval(['print -dpng -r200 ' char(dirName) '/fRespRetroSG.png']);