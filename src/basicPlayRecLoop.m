function output = basicPlayRecLoop(fs,testSignal,audioDevice)
startTic = tic;
bufferSize = fs/20;
[~,nCol] = size(testSignal);
if nCol~=2
   testSignal = [testSignal(:,1) testSignal(:,1)]; 
end
try
playerRecorder = audioPlayerRecorder("SampleRate",fs,"BitDepth","32-bit float", ...
    "SupportVariableSize",true, "BufferSize", bufferSize, "RecorderChannelMapping", ...
    [1 2],"Device",audioDevice,"PlayerChannelMapping",[1 2]);
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
y = zeros(nData,2);
while nextCount < nData - bufferSize
    audioToDevide = testSignal(nextCount+baseIndex,:);
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