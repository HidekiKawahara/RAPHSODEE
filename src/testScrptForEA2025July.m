%% SG test script for EAjuly Sapporo
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

[fileN, pathN] = uigetfile('*.wav');
disp(fileN);
[xIn, fsIn] = audioread(string(pathN)+string(fileN));
infoStr = audioinfo(string(pathN)+string(fileN));
fsOut = 44100;
thLevel = input("Set thershold:");
fHigh = input("Set high-end frequency:");
[xSg,outPrm] = signalSafeguardwithGiantFFTSRC(xIn(:,1),fsIn,fsOut,thLevel,fHigh);
disp(outPrm);
figure(outPrm.varSGFighandle);
title(outPrm.varSGhandle,string(fileN)+" Th:"+num2str(thLevel)+" (dB)");
drawnow
pause(1)
filesgSpecVar = [fileN(1:end-4) num2str(thLevel) 'Var.png'];
print('-dpng','-r200',filesgSpecVar);
pause(1)

figure(outPrm.constSGFighandle);
title(outPrm.constSGhandle,string(fileN)+" Th:"+num2str(thLevel)+" (dB)");
drawnow
pause(1)
filesgSpecConst = [fileN(1:end-4) num2str(thLevel) 'Const.png'];
print('-dpng','-r200',filesgSpecConst);
%%
fs = fsOut;
roomList = {'greathall','octagon','classroom'};
roomID = 2;
roomName = roomList{roomID};
rowList = {'F','B'};% Front and Back
rowID = 1;
rowName = rowList{rowID};
simNoiseList = {'white','pink','red'};
noiseId = 3;
simNoise = simNoiseList{noiseId};

l_data = length(xSg);
xOrg = outPrm.x;
xSgConst = outPrm.xSgConst;
xWhite = outPrm.whiteDetNoise;

%baseDir =  ['/Volumes/HD-CD-A/RIRdatabase/QMLirDB/' roomName 'Omni/Omni/'];
baseDir =  ['~/Downloads/octagonOmni/Omni/'];
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
h = samplingRateConvByDFTwin(hIn,fsRIR,fsOut);
hExt = [h;zeros(l_data-length(h),1)];
hExtF = fft(hExt);
xOrgIn = [xOrg;xOrg;xOrg;xOrg];
xSgIn = [xSg;xSg;xSg;xSg];
xSgConstIn = [xSgConst;xSgConst;xSgConst;xSgConst];
xWhiteIn = [xWhite;xWhite;xWhite;xWhite];
xWhiteIn = xWhiteIn/std(xWhiteIn)*std(xOrgIn);
%
snrList = -6;%0:10:300;
%
switch simNoise
    case simNoiseList{1}
        %xn = randn(length(xOrgIn),1);
        xn = whiteDetNoise(length(xOrgIn),fs);
    case simNoiseList{2}
        %xn = pinknoise(length(xOrgIn));
        xn = pinkDetNoise(length(xOrgIn),fs);
    case simNoiseList{3}
        xn = rednoise(length(xOrgIn),fs);
end
xn = xn/std(xn);
%
fxTest = (0:l_data-1)/l_data*fs;
wavErrList = zeros(l_data,length(snrList),4);
respErrList = zeros(l_data,length(snrList),4);
hEstList = zeros(l_data,length(snrList),4);
%specErrList = zeros(l_data,length(snrList));
telRangeSel = fxTest > 300 & fxTest < 3400;
for ii = 1:length(snrList)
    tic
    yRaw = fftfilt(hExt,xOrgIn);
    addedNoise = xn*10^(-snrList(ii)/20)*std(yRaw);
    yOrg = fftfilt(hExt,xOrgIn)+addedNoise;
    ySg = fftfilt(hExt,xSgIn)+addedNoise;
    ySgConst = fftfilt(hExt,xSgConstIn)+addedNoise;
    yWhite = fftfilt(hExt,xWhiteIn)+addedNoise;
    selIdx = l_data+(1:l_data);
    hEstOrgF = (fft(yOrg(selIdx))./fft(xOrgIn(selIdx)));
    hEstSgF = (fft(ySg(selIdx))./fft(xSgIn(selIdx)));
    hEstConstF = (fft(ySgConst(selIdx))./fft(xSgConstIn(selIdx)));
    hEstWhiteF = (fft(yWhite(selIdx))./fft(xWhiteIn(selIdx)));
    hEstOrg = ifft(hEstOrgF);
    hEstSg = ifft(hEstSgF);
    hEstConst = ifft(hEstConstF);
    hEstWhite = ifft(hEstWhiteF);
    %
    hEstOrgFnxt = (fft(yOrg(selIdx+l_data))./fft(xOrgIn(selIdx)));
    hEstSgFnxt = (fft(ySg(selIdx+l_data))./fft(xSgIn(selIdx)));
    hEstConstFnxt = (fft(ySgConst(selIdx+l_data))./fft(xSgConstIn(selIdx)));
    hEstWhiteFnxt = (fft(yWhite(selIdx+l_data))./fft(xWhiteIn(selIdx)));
    hEstOrgnxt = ifft(hEstOrgFnxt);
    hEstSgnxt = ifft(hEstSgFnxt);
    hEstConstnxt = ifft(hEstConstFnxt);
    hEstWhitenxt = ifft(hEstWhiteFnxt);
    %
    wavErrList(:,ii,:) = [hEstOrg-hExt, hEstSg-hExt, hEstConst-hExt, ...
        hEstWhite-hExt];
    respErrList(:,ii,:) = ...
        [hEstOrgnxt-hEstOrg, hEstSgnxt-hEstSg, hEstConstnxt-hEstConst, ...
        hEstWhitenxt-hEstWhite];
    hEstList(:,ii,:) = ...
        [hEstOrgnxt+hEstOrg, hEstSgnxt+hEstSg, hEstConstnxt+hEstConst, ...
        hEstWhitenxt+hEstWhite]/2;
    toc
end
%%
figure;
set(gcf,"Position",[879   499   423   423])
plot(snrList,20*log10(std(wavErrList(:,:,1)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
hold on;
plot(snrList,20*log10(std(wavErrList(:,:,2)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
plot(snrList,20*log10(std(wavErrList(:,:,3)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
plot(snrList,20*log10(std(wavErrList(:,:,4)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
axis([0 snrList(end) -350 120]);
set(gca,"LineWidth",2,"fontsize",14);
legend("original","SG-spec","SG-const","white");
xlabel("SNR (dB)")
ylabel("error St.Dev. (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle(string(fileN)+" Th:"+num2str(thLevel)+" (dB)");
%figure;plot(snrList,specErrList);grid on;
fileRepNameSD = [fileN(1:end-4) num2str(thLevel) roomName rowName ...
    fileRIR(1:end-4) simNoise 'SD.png'];
pause(1);
print('-dpng','-r200',fileRepNameSD);

figure;
set(gcf,"Position",[879   499   423   423])
plot(snrList,20*log10(max(wavErrList(:,:,1)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
hold on;
plot(snrList,20*log10(max(wavErrList(:,:,2)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
plot(snrList,20*log10(max(wavErrList(:,:,3)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
plot(snrList,20*log10(max(wavErrList(:,:,4)))-20*log10(max(abs(hExt))), ...
    "o-","LineWidth",2);grid on;
axis([0 snrList(end) -350 120]);
set(gca,"LineWidth",2,"fontsize",14);
legend("original","SG-spec","SG-const","white");
xlabel("SNR (dB)")
ylabel("error Max.Dev. (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle(string(fileN)+" Th:"+num2str(thLevel)+" (dB)");
fileRepNameMaxD = [fileN(1:end-4) num2str(thLevel) roomName rowName ...
    fileRIR(1:end-4) simNoise 'MaxD.png'];
print('-dpng','-r200',fileRepNameMaxD);
%%
figure;
set(gcf,"Position",[680   602   560   327]);
semilogx(fxTest,20*log10(abs(fft(xOrg))),"LineWidth",2);grid on;
hold on;
semilogx(fxTest,20*log10(abs(fft(xSgConst))),"LineWidth",2);
semilogx(fxTest,20*log10(abs(fft(xSg))),"LineWidth",2);

xnpw = whiteDetNoise(l_data,fs);
xnpw = xnpw/std(xnpw);
xnpp = pinkDetNoise(l_data,fs);
xnpp = xnpp/std(xnpp);
xnp = rednoise(l_data,fs);
xnp = xnp/std(xnp);
snr = 10;
semilogx(fxTest,20*log10(abs(fft(xnp*10^(-snr/20)*std(xOrgIn)))), ...
    "--","LineWidth",2);
semilogx(fxTest,20*log10(abs(fft(xnpp*10^(-snr/20)*std(xOrgIn)))), ...
    "--","LineWidth",2);
semilogx(fxTest,20*log10(abs(fft(xnpw*10^(-snr/20)*std(xOrgIn)))), ...
    "--","LineWidth",2);
axis([40 fs/2 [-30 60]+50])
set(gca,"LineWidth",2,"fontsize",14);
xlabel("frequency (Hz)")
ylabel("level (dB)")
legend("original","const SG","var SG","red","pink","white","Location","northeast");
title(string(fileN)+" Th:"+num2str(thLevel)+" (dB)  SNR:" ...
    + num2str(snr)+" (dB)","Interpreter","none")
pause(1)
specSGandNoise = [fileN(1:end-4) num2str(thLevel) 'dB_SNR' ...
    num2str(snr) 'dB.png'];
print('-dpng','-r200',specSGandNoise);
%%
snrId = 3;
figure;
semilogx(fxTest,20*log10(abs(fft(hExt))),"LineWidth",2);grid on;
hold on;
semilogx(fxTest,20*log10(abs(fft(respErrList(:,snrId,2))))-3,"LineWidth",2);
semilogx(fxTest,20*log10(abs(fft(respErrList(:,snrId,3))))-3,"LineWidth",2);
semilogx(fxTest,20*log10(abs(fft(respErrList(:,snrId,3))))-4,"LineWidth",2);
axis([40 fs/2 -40 25])
set(gca,"LineWidth",2,"fontsize",14);
xlabel("frequency (Hz)")
ylabel("level (dB)")
legend("RIR spec","var SG Err","con SG Err","white Err");
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle(string(fileN)+" Th:"+num2str(thLevel)+" (dB)  SNR:" + ...
    num2str(snrList(snrId))+" (dB)")
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise num2str(snrList(snrId)) 'dB.png'];
pause(1)
print('-dpng','-r200',fileDiagName);
%%
% indpect detail of the estimated impulse response
% 1: original, 2: var SG, 3: const SG, 4, white sin sum
hrespSg20dBSNR = hEstList(:,2,2);
hrespSgConst20dBSNR = hEstList(:,2,3);
hrespRSW20dBSNR = hEstList(:,2,4); % random sin sum white
hrespSg200dBSNR = hEstList(:,20,2);
txH = (1:length(h))/fs;
tx = (1:length(hrespSg20dBSNR))/fs;
idxH = (1:length(h));
figure;
set(gcf,"Position",[680   602   560   327]);
plot(tx,20*log10(abs(hrespSg20dBSNR)),"LineWidth",2);grid on;hold on;
plot(tx,20*log10(abs(hrespSg200dBSNR)),"LineWidth",2);
axis([0 4 -80 3])
set(gca,"LineWidth",2,"fontsize",14);
xlabel("time (s)")
ylabel("level (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle(string(fileN)+" Th:"+num2str(thLevel));
legend("SNR 20dB","SNR 200dB");
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'pwrChk.png'];
pause(1)
print('-dpng','-r200',fileDiagName);
%%
hCorrSg = fftfilt(hrespSg20dBSNR(fs/2:-1:1),hrespSg20dBSNR);
hCorrSg200 = fftfilt(hrespSg200dBSNR(fs/2:-1:1),hrespSg200dBSNR);
figure;
set(gcf,"Position",[680   602   560   327]);
plot(tx-0.5,20*log10(abs(hCorrSg)/max(hCorrSg)), ...
    tx-0.5,20*log10(abs(hCorrSg200)/max(hCorrSg200)),"LineWidth",2);grid on;
set(gca,"LineWidth",2,"fontsize",14);
xlabel("lag (s)")
ylabel("abs. corr. (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle("Corr. with initial 0.5s " + string(fileN)+" Th:"+num2str(thLevel));
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'CrossCorr.png'];
legend("SNR 20dB","SNR 200dB",Location="southwest");
pause(1)
print('-dpng','-r200',fileDiagName);
%%
figure;
semilogx(20*log10(abs(fft(hrespSg20dBSNR((1:fs))))),"LineWidth",2);grid on;
hold on;
semilogx(20*log10(abs(fft(hrespSg20dBSNR(2.5*fs+(1:fs))))),"LineWidth",2);
axis([40 fs/2 -40 25]);
set(gca,"LineWidth",2,"fontsize",14);
xlabel("frequency (Hz)")
ylabel("level (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle("LTI and err. Sg.var " + string(fileN)+" Th:"+num2str(thLevel));
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'SgVar.png'];
legend("0s-0.5s","3s-3.5s");
pause(1)
print('-dpng','-r200',fileDiagName);
%%
figure;
semilogx(20*log10(abs(fft(hrespSgConst20dBSNR((1:fs))))),"LineWidth",2);grid on;
hold on;
semilogx(20*log10(abs(fft(hrespSgConst20dBSNR(2.5*fs+(1:fs))))),"LineWidth",2);
axis([40 fs/2 -40 25]);
set(gca,"LineWidth",2,"fontsize",14);
xlabel("frequency (Hz)")
ylabel("level (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle("LTI and err. Sg.con " + string(fileN)+" Th:"+num2str(thLevel));
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'SgCon.png'];
legend("0s-0.5s","3s-3.5s");
pause(1)
print('-dpng','-r200',fileDiagName);
%%
figure;
semilogx(20*log10(abs(fft(hrespRSW20dBSNR((1:fs))))),"LineWidth",2);grid on;
hold on;
semilogx(20*log10(abs(fft(hrespRSW20dBSNR(2.5*fs+(1:fs))))),"LineWidth",2);
axis([40 fs/2 -40 25]);
set(gca,"LineWidth",2,"fontsize",14);
xlabel("frequency (Hz)")
ylabel("level (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle("LTI and err. white ssum " + string(fileN)+" Th:"+num2str(thLevel));
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'White.png'];
legend("0s-0.5s","3s-3.5s");
pause(1)
print('-dpng','-r200',fileDiagName);
%%
% これは，うまくいかなかった。エレガントにみえたけれど。。。
tx = (1:length(hrespRSW20dBSNR))'/fs;
figure;
semilogx(20*log10(abs(fft(hrespSg20dBSNR).*exp(-1*tx))),"LineWidth",2);grid on;
hold on;
semilogx(20*log10(abs(fft(hrespSg20dBSNR).*exp(-1*tx(end:-1:1)))),"LineWidth",2);
%%
%soundsc(outPrm.x,fs)
yO = fftfilt(h,outPrm.x);
hF = ifft(fft(yO)./fft(xSg));
hRaw = ifft(fft(yO)./fft(outPrm.x));
tx = (1:length(hF))/fs;
figure;plot(tx,20*log10(abs(hRaw)));grid on;
hold on;plot(tx,20*log10(abs(hF)));
%%
figure;
lhf = length(hF);
fxh = (0:length(hF)-1)'/length(hF)*fs;
semilogx(fxh,20*log10(abs(fft(hF(1:fs/2),lhf))));grid on;
hold on
semilogx(fxh,20*log10(abs(fft(hF(3*fs+(1:fs/2)),lhf))));grid on;
%%
% indpect detail of the estimated impulse response
% with very long signal
% 1: original, 2: var SG, 3: const SG, 4, white sin sum
snrId = 1;
hrespSg20dBSNR = hEstList(:,1,2);
hrespSgConst20dBSNR = hEstList(:,1,3);
hrespRSW20dBSNR = hEstList(:,1,4); % random sin sum white
%hrespSg200dBSNR = hEstList(:,1,2);
txH = (1:length(h))/fs;
tx = (1:length(hrespSg20dBSNR))/fs;
idxH = (1:length(h));
figure;
set(gcf,"Position",[680   602   560   327]);
plot(tx,20*log10(abs(hrespSg20dBSNR)),"LineWidth",2);grid on;hold on;
plot((1:length(hExt))/fs,20*log10(abs(hExt)),"LineWidth",2);
axis([0 4 -80 3])
set(gca,"LineWidth",2,"fontsize",14);
xlabel("time (s)")
ylabel("level (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG SNR:" + num2str(snrList(snrId)) + " (dB)");
subtitle(string(fileN)+" Th:"+num2str(thLevel) + " (dB)","interpreter","none");
legend("estimated","true","fontsize",16);
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'pwrLongSg.png'];
pause(1)
print('-dpng','-r200',fileDiagName);
%%
figure;
semilogx(20*log10(abs(fft(hrespSg20dBSNR((1:fs))))),"LineWidth",2);grid on;
hold on;
semilogx(20*log10(abs(fft(hrespSg20dBSNR(end-3*fs+(1:fs))))),"LineWidth",2);
axis([40 fs/2 -40 25]);
set(gca,"LineWidth",2,"fontsize",14);
xlabel("frequency (Hz)")
ylabel("level (dB)")
title(string(roomName) + string(rowName)+ " " + string(fileRIR)+" " ...
    +string(simNoise)+"-noise-BG");
subtitle("LTI and err. Sg.var " + string(fileN)+" Th:"+num2str(thLevel), ...
    "Interpreter","none");
fileDiagName = [fileN(1:end-4) 'Th' num2str(thLevel) 'dB_' roomName ...
    rowName fileRIR(1:end-4) simNoise 'SgVar.png'];
legend("0s-0.5s","last 1s");
pause(1)
print('-dpng','-r200',fileDiagName);
%% Retrospective sound safefeguarding
%  This is the game changer.
%
%
% hEstOrgF = (fft(yOrg(selIdx))./fft(xOrgIn(selIdx))); % using original
hEstOrgRetroF = (fft(yOrg(selIdx))./fft(xSgIn(selIdx))); % replace with SG
% Above is retrospective use of signal safeguarding
hEstOrgRetro = ifft(hEstOrgRetroF);
%
figure;
tx3s = (1:3*fs)/fs;
txEx = (1:length(hExt))/fs;
plot(tx3s,20*log10(abs(hEstOrg(1:3*fs))));grid on
hold on
plot(tx3s,20*log10(abs(hEstOrgRetro(1:3*fs))));grid on
plot(txEx,20*log10(abs(hExt)));grid on
legend("Original","Retro","True");