%% script for ASJ AA
fs = 44100;
playerRecorder = audioPlayerRecorder("SampleRate",fs,"BitDepth","24-bit integer");
devices = getAudioDevices(playerRecorder);
disp(devices);
t = datetime(datetime,'format','yyyyMMddHHmmss');
dirName = "recordingRetroSG" + string(t);
mkdir(dirName);
thLevel = input("Threshold level (dB):");
fHighIn = input("High freq. limit (Hz):");
testDuration = input("Test duration:");
fileNr = input("Describe input: ","s");
xOrg = zeros(testDuration*fs,1);
%%
devId = input("select device ID:");
pause(3)
sound(0.5*pinknoise(1000));
outputRaw = basicPlayRecLoop(fs,xOrg,devices{devId});
disp("Raw i/o completed");
sound(0.5*pinknoise(1000))
%%
xOutRaw = outputRaw.acquiredSignal(1:length(xOrg),:);
%% 
[y, output] = signalSafeguardwithGiantFFTSRC(xOutRaw(:,2),fs,fs,thLevel,fHighIn,1);
disp("SG completed");
%% 
xInSG = y(:,1);
iRespRawRSG = ifft(fft(xOutRaw)./fft(xInSG));
nData = length(iRespRawRSG);
%%
figure;
%set(gcf,"Position",[1 1 1074 300])
plot((1:fs)/fs,20*log10(abs(iRespRawRSG(1:fs,:))),"LineWidth",2);grid on;
set(gca,"FontSize",13,"LineWidth",2)
axis([0 1 -100 0])
xlabel("time(s)")
ylabel("level (dB)")
eval(['print -dpng -r200 ' char(dirName) '/ipByRetroSG.png']);
%%
[~,peakIdx] = max(abs(iRespRawRSG(:,2)));
headLTI = 0;
lIresp = round(0.2*fs);
fxSG = (0:length(xInSG)-1)'/length(xInSG)*fs;
ltiRawSG = 20*log10(abs(fft(iRespRawRSG(headLTI+(1:lIresp),:),length(xOrg))));
rtvRawSG = 20*log10(abs(fft(iRespRawRSG(round(nData/2)+(1:lIresp),:),length(xOrg))));
ltiRawSGRF = fft(iRespRawRSG(headLTI+(1:lIresp),:),length(xOrg));
rtvRawSGRF = fft(iRespRawRSG(round(nData/2)+(1:lIresp),:),length(xOrg));
%%
lowLimit = output.fLow;
highLimit = output.fHigh;
maxGain = max(ltiRawSG(fxSG>lowLimit & fxSG < highLimit,1));
figure;
semilogx(fxSG,[ltiRawSG(:,1)-maxGain ltiRawSG(:,2)],"LineWidth",2);grid on
hold on;semilogx(fxSG,[rtvRawSG(:,1)-maxGain rtvRawSG(:,2)],"LineWidth",2);grid on;
xline(lowLimit,"LineWidth",2,"Color",'g');
xline(highLimit,"LineWidth",2,"Color",'g');
set(gca,"FontSize",13,"LineWidth",2)
axis([20 fs/2 -60 20])
xlabel("Frequency (Hz)")
ylabel("level (dB)")
legend("Ch1-Ch2 LTI","Raw-SG LIT","Ch1-Ch2 RTV","Raw-SG RTV","Location","northwest")
title(string(fileNr) + " rawRecording, ThLevel: " + num2str(thLevel) + " dB")
eval(['print -dpng -r200 ' char(dirName) '/fRespRetroSG.png']);

%%
figure;
semilogx(fxSG,ltiRawSG(:,1)-ltiRawSG(:,2),"LineWidth",2);grid on;
hold on;semilogx(fxSG,ltiSG(:,1)-ltiSG(:,2),"LineWidth",2);grid on;
xline(lowLimit,"LineWidth",2,"Color",'g');
xline(highLimit,"LineWidth",2,"Color",'g');
set(gca,"FontSize",13,"LineWidth",2)
axis([20 fs/2 -70 10])
xlabel("Frequency (Hz)")
ylabel("level (dB)")
legend("LTI Retro SG","LIT using SG","Location","northwest")
title(string(fileNr) + " rawRecording, ThLevel: " + num2str(thLevel) + " dB")
eval(['print -dpng -r200 ' char(dirName) '/gainRetroSG.png']);