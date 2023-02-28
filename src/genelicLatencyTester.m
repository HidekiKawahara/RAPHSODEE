function output = genelicLatencyTester(iResp,ioBlockL,driverName,fs)

%%

%       Copyright 2023 Hideki Kawahara
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

output.impulseResponse = iResp;
output.ioBlockLength = ioBlockL;
output.driverName = driverName;
output.samplingFrequency = fs;
%
fftl = 2^floor(log2(length(iResp)+2*ioBlockL));
x = iResp*0;
x(1) = 0.05;
hReverbe = fft(x,fftl);
playrec = audioPlayerRecorder(fs,"Device", driverName, ...
    "BitDepth","24-bit integer" ...
    ,"SupportVariableSize",true,"BufferSize",ioBlockL);
testDuration = 5; % seconds
nIteration = round(testDuration/(ioBlockL/fs));
currentIteration = 0;
buffVec = zeros(ioBlockL,1);
%------- generate test input signal LNN a kind of -----
simInput = zeros(ioBlockL*nIteration+fftl,1);
unitSig = zeros(fs,1);
tt = (1:fs)'/fs;
for ii = 1:fs/2
    unitSig = unitSig + (1/sqrt(ii))*sin(2*pi*ii*tt+randn(1,1)*2*pi);
end
for ii = 1:testDuration
    simInput((ii-1)*fs+(1:fs)) = unitSig;
end
simInput = simInput/max(abs(simInput))*0.5 + randn(length(simInput),1)/100000;
%------------
soundBuff = zeros(fftl,1);
timeBuffer = zeros(nIteration,1);
overRunRecord = zeros(nIteration,1);
underRunRecord = zeros(nIteration,1);
caputuredSignal = zeros(ioBlockL*nIteration+fftl,1);
pause(1)
while currentIteration < nIteration
    currentIteration = currentIteration + 1;
    tmpOut = buffVec+simInput((currentIteration-1)*ioBlockL+(1:ioBlockL));
    [readData,nUnderruns,nOverruns] = ...
        playrec([tmpOut tmpOut]);
    caputuredSignal((currentIteration-1)*ioBlockL+(1:ioBlockL)) = readData;
    underRunRecord(currentIteration) = nUnderruns;
    overRunRecord(currentIteration) = nOverruns;
    currentTime1 = datetime('now');
    soundBuff(1:fftl-ioBlockL) = soundBuff(ioBlockL+1:fftl);
    soundBuff(fftl-ioBlockL+1:fftl) = readData;
    resp = real(ifft(fft(soundBuff).*hReverbe));
    buffVec = resp(fftl-ioBlockL+1:fftl);
    currentTime2 = datetime('now');
    timeBuffer(currentIteration) = seconds(currentTime2-currentTime1);
end
%
respF = fft(caputuredSignal(fs+(1:fs)))./fft(simInput(fs+(1:fs)));
respFix = respF;
respFix(1) = 0;
respT = ifft(respFix);
tt = (1:fs)'/fs;
[~,maxidx] = max(abs(respT));
[~,secIdx] = max(abs(respT).*(tt > tt(maxidx)+0.01));
latency = tt(secIdx)-tt(maxidx);
release(playrec);
output.overRunRecord = overRunRecord;
output.underRunRecord = underRunRecord;
output.simInput = simInput;
output.caputuredSignal = caputuredSignal;
output.timeBuffer = timeBuffer;
output.totalOverRun = sum(output.overRunRecord);
output.totalUnderRun = sum(output.underRunRecord);
output.medianLatencyMs = median(timeBuffer)*1000;
output.latencyMs = latency*1000;
return;