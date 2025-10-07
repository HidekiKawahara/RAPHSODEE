function output = basicPlayRecLoopMono(fs,testSignal,audioDevice)
startTic = tic;
bufferSize = fs/20;
try
playerRecorder = audioPlayerRecorder("SampleRate",fs,"BitDepth","32-bit float", ...
    "SupportVariableSize",true, "BufferSize", bufferSize, "Device",audioDevice);
catch
    disp("Thid device does not support specified playRec functionality.");
    output.elapsedTime = toc(startTic);
    return;
end
nData = length(testSignal);
nextCount = 0;
baseIndex = 1:bufferSize;
totalUnderrun = 0;
totalOverrun = 0;
y = zeros(nData,1);
while nextCount < nData - bufferSize
    audioToDevide = testSignal(nextCount+baseIndex);
    [audioFromDevice,numUnderrun,numOverrun] = playerRecorder(audioToDevide);
    totalUnderrun = totalUnderrun + numUnderrun;
    totalOverrun = totalOverrun + numOverrun;
    y(nextCount+baseIndex,:) = audioFromDevice;
    nextCount = nextCount+bufferSize;
end
release(playerRecorder);
output.elapsedTime = toc(startTic);
output.playerRecorder = playerRecorder;
output.acquiredSignal = y;
output.outputSignal = testSignal;
output.totalUnderrun = totalUnderrun;
output.totalOverrun = totalOverrun;
end