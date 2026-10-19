classdef acouticCheckerByGFFTRefactored_v2_exportedR < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        LeftPanel                     matlab.ui.container.Panel
        snapshotGUIButton             matlab.ui.control.Button
        OutputchannelButtonGroup      matlab.ui.container.ButtonGroup
        LAndRButton                   matlab.ui.control.RadioButton
        RButton                       matlab.ui.control.RadioButton
        LButton                       matlab.ui.control.RadioButton
        SpectrumdisplayButtonGroup    matlab.ui.container.ButtonGroup
        Octave12thRadioButton         matlab.ui.control.RadioButton
        Octave6thRadioButton          matlab.ui.control.RadioButton
        RawDFTButton                  matlab.ui.control.RadioButton
        LoadCyclicData                matlab.ui.control.Button
        STOP                          matlab.ui.control.Button
        START                         matlab.ui.control.Button
        InputChannelDropDown          matlab.ui.control.DropDown
        SamplingfrequencymenuLabel_2  matlab.ui.control.Label
        SamplingFrequencyDropDown     matlab.ui.control.DropDown
        SamplingfrequencymenuLabel    matlab.ui.control.Label
        SystemStatusLabel             matlab.ui.control.Label
        StatusLabel                   matlab.ui.control.Label
        IOdeviceselectionDropDown     matlab.ui.control.DropDown
        IOdeviceselectionLabel        matlab.ui.control.Label
        CenterPanel                   matlab.ui.container.Panel
        RecordingLamp                 matlab.ui.control.Lamp
        RecordingLampLabel            matlab.ui.control.Label
        CalibrationCh2Panel           matlab.ui.container.Panel
        CALValuedBEditField_2         matlab.ui.control.NumericEditField
        CALValuedBEditField_2Label    matlab.ui.control.Label
        MeasuredSPLDropDown_2         matlab.ui.control.DropDown
        MeasuredSPLDropDown_2Label    matlab.ui.control.Label
        CalibrateButton_2             matlab.ui.control.Button
        CalibrationCh1Panel           matlab.ui.container.Panel
        CALValuedBEditField           matlab.ui.control.NumericEditField
        CALValuedBEditFieldLabel      matlab.ui.control.Label
        MeasuredSPLDropDown           matlab.ui.control.DropDown
        MeasuredSPLDropDownLabel      matlab.ui.control.Label
        CalibrateButton               matlab.ui.control.Button
        Voicememo5sButton             matlab.ui.control.Button
        SetworkingfolderButton        matlab.ui.control.Button
        RightPanel                    matlab.ui.container.Panel
        BPFfcSlider                   matlab.ui.control.Slider
        BPFfcSliderLabel              matlab.ui.control.Label
        SnapshotLevelButton           matlab.ui.control.Button
        bpfPowerAxes                  matlab.ui.control.UIAxes
        BandSPLAxes                   matlab.ui.control.UIAxes
        SpectrumUIAxes                matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end

    events
        RequestDisplayUpdate
    end

    properties (Access = private)
        StreamTimer            % matlab timer for polling
        SampleRate = 48000     % Sampling frequency (Hz)
        FrameLength = 1024     % default value (must be integer in samples)
        nChannel = 1
        CyclicData             % Cyclic Audio Data, structure
        SupportedSamplingRate
        FileName
        FileFolder
        FilePath
        nChannelOut = 2
        maximumDuration = 300  % seconds
        DisplayUpdateListener  % listener handle
        BandSPLDisplay         % Band spectrum display
        raw1kHzBandLevel = 0
        captured1kHzBandLevel = 0
        calibrationData
        InitialCondition
        WorkingDirectory
        VoiceMemo
    end

    properties (Access = public)
        AS                     % audiostreamer object
        StreamingEnable
        streamingData
        AnalysisData
        SpectrumDisplay
        powerRespDisplay
    end

    methods (Access = private)

        function populateAudioDevices(app)
            app.CyclicData = struct;
            app.SamplingFrequencyDropDown.Enable = 'off';
            app.InputChannelDropDown.Enable = "off";
            app.LoadCyclicData.Enable = 'off';
            app.START.Enable = 'off';
            app.STOP.Enable = 'off';
            app.StreamingEnable = false;
            try
                apr = audioPlayerRecorder;
            catch ME
                uialert(app.UIFigure, "Cannot create audioPlayerRecorder: " + ME.message, "Audio Error");
                return
            end

            % Get full-duplex devices
            devices = getAudioDevices(apr);
            release(apr);

            % If no devices found, show message and populate with empty entry
            if isempty(devices)
                app.IOdeviceselectionDropDown.Items = {'<No full-duplex devices found>'};
                app.IOdeviceselectionDropDown.Enable = 'off';
                app.StatusLabel.Text = 'No full-duplex audio devices available';
                return
            end

            devices{end+1} = ' ---- ';

            % Populate DropDown and enable it
            app.IOdeviceselectionDropDown.Items = devices(:);
            app.IOdeviceselectionDropDown.Value = app.IOdeviceselectionDropDown.Items{end};
            app.IOdeviceselectionDropDown.Enable = 'on';
            app.IOdeviceselectionLabel.Enable = 'on';
            app.StatusLabel.Text = "Select IO device";
        end

        function results = SelectFileButtonPushed(app)
            [fname, fpath] = uigetfile({'*.wav';'*.aif';'*.WAV';'*.AIF'}, ...
                'Select file');
            if isequal(fname, 0)
                results = 0;
                return; % user canceled
            end
            app.FileName   = fname;
            app.FileFolder = fpath;
            app.FilePath   = fullfile(fpath, fname);
            results = 1;
        end

        function results = initializeGraphics(app)
            % Spectrum plot
            fx = (0:app.SampleRate-1)';
            app.AnalysisData.bpfFc = 1000;
            ltiGain = zeros(app.SampleRate, 1);
            rtvGain = zeros(app.SampleRate, 1) - 30;
            sdtiGain = zeros(app.SampleRate, 1) - 20;
            app.SpectrumDisplay.ltiHandle = ...
                semilogx(app.SpectrumUIAxes, fx, ltiGain, "LineWidth", 2);
            hold(app.SpectrumUIAxes, "on");
            app.SpectrumDisplay.rtvHandle = ...
                semilogx(app.SpectrumUIAxes, fx, rtvGain, "LineWidth", 2);
            app.SpectrumDisplay.sdtiHandle = ...
                semilogx(app.SpectrumUIAxes, fx, sdtiGain, "LineWidth", 2);
            axis(app.SpectrumUIAxes, [20 app.SampleRate/2 -80 20]);
            app.SpectrumDisplay.fc1kHandle = ...
                xline(app.SpectrumUIAxes, 1000, "g", "1000", "LineWidth", 3);
            grid(app.SpectrumUIAxes, "on");
            set(app.SpectrumUIAxes, "LineWidth", 2);
            ylabel(app.SpectrumUIAxes, "Level (rel. dB)");
            xlabel(app.SpectrumUIAxes, "frequency (Hz)");
            title(app.SpectrumUIAxes, "  ");
            drawnow;

            % BandSPL plot
            signalLevel = zeros(app.SampleRate, 1);
            disturbanceLevel = zeros(app.SampleRate, 1) - 100;
            snapshotLevel = zeros(app.SampleRate, 1) - 100;
            app.BandSPLDisplay.signalLevelHandle = ...
                semilogx(app.BandSPLAxes, fx, signalLevel, "LineWidth", 2);
            hold(app.BandSPLAxes, "on");
            app.BandSPLDisplay.distrbanceLevelHandle = ...
                semilogx(app.BandSPLAxes, fx, disturbanceLevel, "LineWidth", 2);
            app.BandSPLDisplay.indirectLevelHandle = ...
                semilogx(app.BandSPLAxes, fx, disturbanceLevel, "LineWidth", 2);
            app.BandSPLDisplay.snapshotLevelHandle = ...
                semilogx(app.BandSPLAxes, fx, snapshotLevel, "g", "LineWidth", 2);
            app.BandSPLDisplay.fc1kHandle = ...
                xline(app.BandSPLAxes, 1000, "g", "1000", "LineWidth", 3);
            axis(app.BandSPLAxes, [20 app.SampleRate/2 -110 0]);
            grid(app.BandSPLAxes, "on");
            set(app.BandSPLAxes, "LineWidth", 2);
            ylabel(app.BandSPLAxes, "SPL CAL.NO (dB)");
            xlabel(app.BandSPLAxes, "frequency (Hz)");
            title(app.BandSPLAxes, "  ");

            % Power decay plot
            tresp = (1:app.SampleRate)' / app.SampleRate;
            app.powerRespDisplay.pResp1kHandle = ...
                plot(app.bpfPowerAxes, tresp, 0*tresp, "LineWidth", 2);
            hold(app.bpfPowerAxes, "on");
            app.powerRespDisplay.pRespFcHandle = ...
                plot(app.bpfPowerAxes, tresp, 0*tresp, "g", "LineWidth", 2);
            grid(app.bpfPowerAxes, "on");
            axis(app.bpfPowerAxes, [0 1 -70 0]);
            set(app.bpfPowerAxes, "LineWidth", 2);
            xlabel(app.bpfPowerAxes, "time (s)");
            ylabel(app.bpfPowerAxes, "level (dB)");
            title(app.bpfPowerAxes, "Power response");
            drawnow;
            results = 1;
        end

        function results = updateGraphics(app)
            % Spectrum display
            fs = app.SampleRate;
            fftl = app.AnalysisData.impulseLength;
            fx = (0:fftl-1)' / fftl * fs;
            app.AnalysisData.frequencyAxis = fx;
            app.AnalysisData.frequencyAxix = fx; % Kept for backward compatibility
            ltiGain = zeros(fftl, 1);
            rtvGain = zeros(fftl, 1) - 30;
            sdtiGain = zeros(fftl, 1) - 120;
            set(app.SpectrumDisplay.ltiHandle, "XData", fx, "YData", ltiGain);
            set(app.SpectrumDisplay.rtvHandle, "XData", fx, "YData", rtvGain);
            set(app.SpectrumDisplay.sdtiHandle, "XData", fx, "YData", sdtiGain);

            % Slider for spectrum display (guarded for +/-11 filter-bank margin)
            app.BPFfcSlider.Enable = "on";
            nBands = app.AnalysisData.gfStr12.nChannel;
            fcSet = app.AnalysisData.gfStr12.fcSet;
            [~, idx1kHz] = min(abs(fcSet - 1000));
            minIdx = min(12, nBands);
            maxIdx = max(minIdx, nBands - 11);
            app.BPFfcSlider.Limits = [minIdx, maxIdx];
            app.BPFfcSlider.Value = min(max(idx1kHz, minIdx), maxIdx);

            % BandSPL display
            signalLevel = zeros(fftl, 1);
            indirectLevel = zeros(fftl, 1) - 200;
            disturbanceLevel = zeros(fftl, 1) - 200;
            snapshotLevel = zeros(fftl, 1) - 200;
            set(app.BandSPLDisplay.signalLevelHandle, "XData", fx, "YData", signalLevel);
            set(app.BandSPLDisplay.distrbanceLevelHandle, "XData", fx, "YData", disturbanceLevel);
            set(app.BandSPLDisplay.indirectLevelHandle, "XData", fx, "YData", indirectLevel);
            set(app.BandSPLDisplay.snapshotLevelHandle, "XData", fx, "YData", snapshotLevel);
            app.SpectrumDisplay.frequencyAxis = (0:app.AnalysisData.dataLength-1)' ...
                / app.AnalysisData.dataLength * fs;
            drawnow limitrate;
            results = 1;
        end

        function analysisAndUpdateGraphicsHandler(app, ~, ~)
            if isempty(app.AS)
                return;
            elseif ~isPlaying(app.AS)
                if ~isempty(app.StreamTimer) && isvalid(app.StreamTimer)
                    stop(app.StreamTimer);
                else
                    pause(0.1);
                    return;
                end
                pause(0.1);
                stopHousekeeping(app);
                return;
            end

            fs = app.SampleRate;
            fxx = app.SpectrumDisplay.frequencyAxis;
            fxImpulse = app.AnalysisData.frequencyAxis;
            fcSet = app.AnalysisData.gfStr6.fcSet;
            fcSetL = fcSet * 2^(-1/2);
            fcSetH = fcSet * 2^(1/2);
            capturedData = read(app.AS);
            nInData = size(capturedData, 1);
            cycleLength = app.streamingData.tmpBufferLength;
            sourceBuffer = app.streamingData.wholeSource;

            if app.AnalysisData.headIdx + nInData > size(sourceBuffer, 1) || ...
               app.AnalysisData.headIdx + cycleLength > size(sourceBuffer, 1)
                if ~isempty(app.StreamTimer) && isvalid(app.StreamTimer)
                    stop(app.StreamTimer);
                end
                pause(0.2);
                stop(app.AS);
                stopHousekeeping(app);
                return;
            end

            sourceId = 1;
            if app.RButton.Value == 1 && app.nChannelOut == 2
                sourceId = 2;
            end

            if nInData > 0
                app.streamingData.capturedBuffer(app.AnalysisData.headIdx + (1:nInData), :) = ...
                    capturedData;
                app.AnalysisData.headIdx = app.AnalysisData.headIdx + nInData;

                if app.AnalysisData.headIdx >= cycleLength
                    tmpHeadIdx = app.AnalysisData.headIdx;
                    rangeIdx = (tmpHeadIdx - cycleLength + 1):tmpHeadIdx;
                    if app.AnalysisData.firstTime
                        app.AnalysisData.captured = sourceBuffer(rangeIdx, sourceId);
                        app.InitialCondition.captured = ...
                            app.streamingData.capturedBuffer(rangeIdx, 1);
                        app.InitialCondition.source = ...
                            sourceBuffer(rangeIdx, sourceId);
                    else
                        app.AnalysisData.captured = ...
                            app.streamingData.capturedBuffer(rangeIdx, 1);
                    end
                    app.AnalysisData.source = sourceBuffer(rangeIdx, sourceId);
                    epsSpec = 1e-12;
                    capturedSpectrum = fft(app.AnalysisData.captured, app.AnalysisData.Nfft);
                    capturedPower = abs(capturedSpectrum).^2 / app.AnalysisData.Nfft^2;
                    octaveBandLevel = interp1(fxx, cumsum(capturedPower), fcSetH, "linear", "extrap") ...
                        - interp1(fxx, cumsum(capturedPower), fcSetL, "linear", "extrap");
                    octaveBandLevel = max(octaveBandLevel, epsSpec);

                    app.AnalysisData.transferFunction = ...
                        capturedSpectrum ...
                        ./ fft(app.AnalysisData.source + epsSpec, app.AnalysisData.Nfft);
                    app.AnalysisData.impulseResponse = ...
                        ifft(app.AnalysisData.transferFunction, "symmetric");

                    fftl = app.AnalysisData.impulseLength;
                    impulseResp = app.AnalysisData.impulseResponse;
                    Nimp = numel(impulseResp);
                    [~, idxPeak] = max(abs(impulseResp));
                    biasL = app.AnalysisData.halfDataLength;

                    % Modular index wrapping to prevent out-of-bounds access
                    indIndices = mod((0:fftl-1) + idxPeak + round(fs*0.015) - 1, Nimp) + 1;
                    rtvIndices = mod((0:fftl-1) + biasL - 1, Nimp) + 1;

                    ltiSpec = fft(impulseResp(1:fftl));
                    indirectSpec = fft(impulseResp(indIndices));
                    rtvSpec = fft(impulseResp(rtvIndices));

                    ltiPw = abs(ltiSpec).^2;
                    indirectPw = abs(indirectSpec).^2;
                    rtvPw = abs(rtvSpec).^2;

                    octaveLTIBandLvl = max(interp1(fxImpulse, cumsum(ltiPw), fcSetH, "linear", "extrap") ...
                        - interp1(fxImpulse, cumsum(ltiPw), fcSetL, "linear", "extrap"), epsSpec);
                    octaveRTVBandLvl = max(interp1(fxImpulse, cumsum(rtvPw), fcSetH, "linear", "extrap") ...
                        - interp1(fxImpulse, cumsum(rtvPw), fcSetL, "linear", "extrap"), epsSpec);
                    octaveIndBandLvl = max(interp1(fxImpulse, cumsum(indirectPw), fcSetH, "linear", "extrap") ...
                        - interp1(fxImpulse, cumsum(indirectPw), fcSetL, "linear", "extrap"), epsSpec);

                    relRTVBandLvl = octaveRTVBandLvl ./ octaveLTIBandLvl;
                    relRTVBandLvlDB = 10*log10(relRTVBandLvl);
                    relIndBandLvl = octaveIndBandLvl ./ octaveLTIBandLvl;
                    relIndBandLvlDB = 10*log10(relIndBandLvl);

                    if app.RawDFTButton.Value == 1
                        frequencyAxis = app.AnalysisData.frequencyAxis;
                        ltiGain = 20*log10(max(abs(ltiSpec), epsSpec));
                        rtvGain = 20*log10(max(abs(rtvSpec), epsSpec));
                    elseif app.Octave6thRadioButton.Value == 1
                        frequencyAxis = app.AnalysisData.gfStr6.fcSet;
                        ltiGain = ...
                            10*log10(max(app.AnalysisData.gfStr6.filterWeight * abs(ltiSpec).^2, epsSpec));
                        rtvGain = ...
                            10*log10(max(app.AnalysisData.gfStr6.filterWeight * abs(rtvSpec).^2, epsSpec));
                    else
                        frequencyAxis = app.AnalysisData.gfStr12.fcSet;
                        ltiGain = ...
                            10*log10(max(app.AnalysisData.gfStr12.filterWeight * abs(ltiSpec).^2, epsSpec));
                        rtvGain = ...
                            10*log10(max(app.AnalysisData.gfStr12.filterWeight * abs(rtvSpec).^2, epsSpec));
                    end

                    set(app.SpectrumDisplay.ltiHandle, "XData", frequencyAxis, "YData", ltiGain);
                    set(app.SpectrumDisplay.rtvHandle, "XData", frequencyAxis, "YData", rtvGain);

                    octaveBandLevelDB = 10*log10(octaveBandLevel) ...
                        + app.calibrationData.ch1CalibrationDb + 3;
                    set(app.BandSPLDisplay.signalLevelHandle, "XData", fcSet, "YData", octaveBandLevelDB);
                    set(app.BandSPLDisplay.distrbanceLevelHandle, "XData", fcSet, "YData", ...
                        octaveBandLevelDB + relRTVBandLvlDB);
                    set(app.BandSPLDisplay.indirectLevelHandle, "XData", fcSet, "YData", ...
                        octaveBandLevelDB + relIndBandLvlDB);

                    % Power decay calculation
                    resp = ifft(app.AnalysisData.fc1kHzFilterF .* ltiSpec);
                    powerResp1k = 10*log10(max(sum(abs(resp).^2, 2), epsSpec));
                    powerResp1k = powerResp1k - max(powerResp1k);
                    tx = (1:fftl) / app.SampleRate;
                    set(app.powerRespDisplay.pResp1kHandle, "XData", tx(1:10:end), "YData", powerResp1k(1:10:end));

                    resp = ifft(app.AnalysisData.fcSliderFilterF .* ltiSpec);
                    powerRespFc = 10*log10(max(sum(abs(resp).^2, 2), epsSpec));
                    powerRespFc = powerRespFc - max(powerRespFc);
                    set(app.powerRespDisplay.pRespFcHandle, "XData", tx(1:10:end), "YData", powerRespFc(1:10:end));

                    app.AnalysisData.firstTime = false;
                    drawnow limitrate;
                end
            end
        end

        function results = audioStreamingLoopSlim(app)
            cycleLength = app.streamingData.tmpBufferLength;
            sourceBuffer = app.streamingData.wholeSource;
            app.AnalysisData.Nfft = cycleLength;
            app.AnalysisData.firstTime = true;
            record(app.AS);
            play(app.AS, sourceBuffer, "non-blocking");
            app.AnalysisData.headIdx = 0; % 0-based accumulator so 1st write is 1:nInData
            app.StreamTimer = timer( ...
                ExecutionMode = "fixedRate", ...
                Period = 0.06, ...
                TimerFcn = @(~,~) timerProcessUpdate(app));
            start(app.StreamTimer);
            app.SnapshotLevelButton.Enable = "on";
            results = 1;
        end

        function timerProcessUpdate(app)
            notify(app, 'RequestDisplayUpdate');
        end

        function results = updateSamplingRate(app)
            fs = app.CyclicData.samplingFrequency;
            matchedID = find(app.SupportedSamplingRate == fs, 1, 'first');
            if isempty(matchedID)
                results = 0;
                disp(num2str(fs) + " Hz sampling is not supported.");
                return;
            else
                app.SamplingFrequencyDropDown.Value = ...
                    app.SamplingFrequencyDropDown.Items{matchedID};
                app.SampleRate = fs;
            end
            results = 1;
        end

        function results = stopHousekeeping(app)
            try
                if isprop(app, 'StreamTimer') && ~isempty(app.StreamTimer) && isvalid(app.StreamTimer)
                    stop(app.StreamTimer);
                    delete(app.StreamTimer);
                    app.StreamTimer = [];
                end
            catch ME
                disp("Error cleaning StreamTimer: " + ME.message);
            end
            app.StatusLabel.Text = 'Stopped: START or Select another file.';
            app.STOP.Enable = 'off';
            app.START.Enable = 'on';
            app.LoadCyclicData.Enable = 'on';
            app.Voicememo5sButton.Enable = "on";

            if app.nChannelOut > 1
                app.OutputchannelButtonGroup.Enable = "on";
                app.LAndRButton.Enable = "on";
                app.LButton.Enable = "on";
                app.RButton.Enable = "on";
            end
            results = 1;
        end

        function results = initializeWorkingDirectory(app)
            app.StatusLabel.Text = "Click 'Set working folder' to begin.";
            app.SetworkingfolderButton.Enable = "on";
            app.IOdeviceselectionDropDown.Enable = "off";
            app.IOdeviceselectionLabel.Enable = "off";
            drawnow;
            results = 1;
        end

        function results = createAudioStreamer(app)
            results = 0;
            dev = app.IOdeviceselectionDropDown.Value;
            if isempty(dev) || strcmp(dev, ' ---- ') || strcmp(dev, '<No full-duplex devices found>')
                uialert(app.UIFigure, "Please select a valid audio IO device.", "Audio Error");
                return;
            end
            app.StatusLabel.Text = ['Opening: ' dev];
            switch computer
                case {'MACI64','MACA64'}
                    driverStr = 'CoreAudio';
                case 'PCWIN64'
                    driverStr = 'ASIO';
                otherwise
                    uialert(app.UIFigure, "Unsupported operating system: " + string(computer), "Audio Error");
                    return
            end

            % Create audiostreamer for full-duplex
            try
                app.AS = audiostreamer("Mode", "full-duplex", "SampleRate", app.SampleRate, ...
                    "Player", dev, "Recorder", dev, "DeviceBitFormat", "int24", ...
                    "PlayerChannels", (1:app.nChannelOut), "RecorderChannels", (1:app.nChannel), ...
                    "DeviceBufferSize", app.FrameLength, "Driver", driverStr);
            catch ME
                uialert(app.UIFigure, "Cannot start audiostreamer: " + ME.message, "Audio Error");
                return
            end
            results = 1;
        end

        function results = enableCalibration(app)
            results = 0;
            switch app.nChannel
                case 1
                    app.calibrationData.ch1 = 0;
                    app.CalibrationCh1Panel.Enable = "on";
                    app.CalibrateButton.Enable = "on";
                    app.CALValuedBEditField.Enable = "off";
                    app.CALValuedBEditFieldLabel.Enable = "off";
                case 2
                    app.calibrationData.ch1 = 0;
                    app.CalibrationCh1Panel.Enable = "on";
                    app.CalibrateButton.Enable = "on";
                    app.CALValuedBEditField.Enable = "off";
                    app.CALValuedBEditFieldLabel.Enable = "off";
                    app.calibrationData.ch2 = 0;
                    app.CalibrationCh2Panel.Enable = "on";
                    app.CalibrateButton_2.Enable = "on";
                    app.CALValuedBEditField_2.Enable = "off";
                    app.CALValuedBEditField_2Label.Enable = "off";
                otherwise
                    app.calibrationData.ch1 = 0; % safeguard
                    return;
            end
            app.calibrationData.initialized = true;
            results = 1;
        end

        function results = generateCalibrationSignal(app)
            tic;
            fs = app.SampleRate;
            tt = (1:fs)' / fs;
            fc = 1000; % center frequency
            fcList = round(fc*2^(-1/2)):round(fc*2^(1/2));
            phases = rand(1, numel(fcList));
            x = sum(sin(2*pi * (tt * fcList + phases)), 2);
            unitSource = x / std(x) * 0.1;

            if app.nChannelOut == 2
                tmpBuffer = [unitSource, unitSource];
                if app.LButton.Value == 1
                    tmpBuffer(:, 2) = 0;
                elseif app.RButton.Value == 1
                    tmpBuffer(:, 1) = 0;
                end
            else
                tmpBuffer = unitSource;
            end

            nCycles = floor(app.maximumDuration * fs / size(unitSource, 1));
            x = repmat(tmpBuffer, nCycles, 1);
            toc;
            app.calibrationData.calibrationSignal = x;
            results = 1;
        end

        function results = designFcOctaveSet(app)
            fc = app.AnalysisData.bpfFc;
            fcSet = app.AnalysisData.gfStr12.fcSet;
            nBands = numel(fcSet);
            [~, idxFc] = min(abs(fc - fcSet));
            bandIndices = min(max(idxFc + (-11:11), 1), nBands);
            app.AnalysisData.tmpBPFset = ...
                app.AnalysisData.gfStr12.quadfilterWeight(:, bandIndices);
            results = 1;
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Hide the axes toolbar on the UIAxes
            app.SpectrumUIAxes.Toolbar.Visible = 'off';
            app.DisplayUpdateListener = addlistener(app, 'RequestDisplayUpdate', ...
                @(src,ev) analysisAndUpdateGraphicsHandler(app, src, ev));
            app.calibrationData.initialized = false;
            app.calibrationData.ch1CalibrationDb = 0;
            app.calibrationData.ch2CalibrationDb = 0;
            app.calibrationData.ch1Measured = 94;
            app.InitialCondition = struct('source', [], 'captured', []);
            app.WorkingDirectory = [];
            initializeGraphics(app);
            initializeWorkingDirectory(app);
        end

        % Value changed function: IOdeviceselectionDropDown
        function IOdeviceselectionDropDownValueChanged(app, event)
            selected = app.IOdeviceselectionDropDown.Value;
            if isempty(selected) || strcmp(selected, ' ---- ') || strcmp(selected, '<No full-duplex devices found>')
                return;
            end
            app.StatusLabel.Text = "Selected: " + selected;
            candidateRates = [8000 11025 16000 22050 32000 44100 48000 88200 96000 192000];
            supported = zeros(length(candidateRates), 1);
            ii = 0;
            tmpOutCh = 0;
            for r = candidateRates
                try
                    % Try to create an audiostreamer in full-duplex with this rate
                    as = audiostreamer("full-duplex", r, ...
                        Player = selected, ...
                        Recorder = selected);
                    tmpOutCh = as.MaxPlayerChannels;
                    release(as); % release if created
                    ii = ii + 1;
                    supported(ii) = r;
                catch
                    % Not supported at this rate
                end
            end
            supported = supported(1:ii);
            if isempty(supported)
                uialert(app.UIFigure, "No supported full-duplex sampling rates found for " + selected, "Device Error");
                return;
            end
            app.SupportedSamplingRate = supported;
            supportedStringArray = strtrim(string(num2str(supported)))';
            supportedStringArray{end+1} = ' ---- ';
            app.SamplingFrequencyDropDown.Items = supportedStringArray;
            app.SamplingFrequencyDropDown.Value = app.SamplingFrequencyDropDown.Items{end};
            app.SamplingFrequencyDropDown.Enable = 'on';
            app.SamplingfrequencymenuLabel.Enable = 'on';
            app.IOdeviceselectionDropDown.Enable = "off";
            if tmpOutCh > 1
                app.nChannelOut = 2; % up to stereo
                app.OutputchannelButtonGroup.Enable = "on";
                app.LAndRButton.Enable = "on";
                app.LButton.Enable = "on";
                app.RButton.Enable = "on";
            else
                app.nChannelOut = 1;
            end
            app.StatusLabel.Text = "Select Sampling frequency";
        end

        % Value changed function: SamplingFrequencyDropDown
        function SamplingFrequencyDropDownValueChanged(app, event)
            value = app.SamplingFrequencyDropDown.Value;
            if isempty(value) || strcmp(value, ' ---- ')
                return;
            end
            selectedDevice = app.IOdeviceselectionDropDown.Value;
            r = str2double(value);
            if isnan(r)
                return;
            end
            app.AS = audiostreamer("full-duplex", r, ...
                Player = selectedDevice, ...
                Recorder = selectedDevice);
            app.StatusLabel.Text = "Sampling frequency is " + value + " Hz";
            app.SampleRate = r;
            pause(0.2);
            app.SamplingFrequencyDropDown.Enable = "off";
            app.StatusLabel.Text = "Select Number of Input Channels";
            maxRecCh = app.AS.MaxRecorderChannels;
            release(app.AS);

            if maxRecCh > 1
                app.InputChannelDropDown.Enable = "on";
                app.SamplingfrequencymenuLabel_2.Enable = "on";
            else
                app.InputChannelDropDown.Value = app.InputChannelDropDown.Items{1};
                app.InputChannelDropDown.Enable = "off";
                InputChannelDropDownValueChanged(app, event);
            end
        end

        % Value changed function: InputChannelDropDown
        function InputChannelDropDownValueChanged(app, event)
            itemId = app.InputChannelDropDown.ValueIndex;
            app.nChannel = str2double(app.InputChannelDropDown.ItemsData{itemId});
            if app.nChannel > 0
                app.StatusLabel.Text = "Number of Input Channels: " + app.InputChannelDropDown.Value;
                pause(0.2);
                app.InputChannelDropDown.Enable = "off";
                app.StatusLabel.Text = "Ready to load cyclic data";
                app.LoadCyclicData.Enable = 'on';
                app.SpectrumdisplayButtonGroup.Enable = "on";
                app.Octave12thRadioButton.Enable = "on";
                app.Octave6thRadioButton.Enable = "on";
                app.RawDFTButton.Enable = "on";
                app.Voicememo5sButton.Enable = "on";
            end
        end

        % Button pushed function: START
        function STARTButtonPushed(app, event)
            app.Voicememo5sButton.Enable = "off";
            if createAudioStreamer(app) == 0
                app.Voicememo5sButton.Enable = "on";
                return;
            end

            % Set timer environment
            app.streamingData.headIndex = 1;
            app.streamingData.playIdx = 1;
            app.StatusLabel.Text = "Streaming ";
            app.StreamingEnable = true;
            fs = app.SampleRate;
            app.streamingData.capturedBuffer = ...
                zeros(fs * app.maximumDuration, app.nChannel);

            % Ensure mono column vector before channel routing
            unitSource = app.CyclicData.data(:, 1);
            sz = size(unitSource);
            if app.nChannelOut == 2
                tmpBuffer = [unitSource, unitSource];
                if app.LButton.Value == 1
                    tmpBuffer(:, 2) = 0;
                elseif app.RButton.Value == 1
                    tmpBuffer(:, 1) = 0;
                end
            else
                tmpBuffer = unitSource;
            end

            nCycles = floor(app.maximumDuration * fs / sz(1));
            app.streamingData.wholeSource = repmat(tmpBuffer, nCycles, 1);
            app.streamingData.tmpBuffer = tmpBuffer;
            app.streamingData.tmpBufferLength = size(tmpBuffer, 1);
            app.streamingData.total_numUnderrun = 0;
            app.streamingData.total_numOverrun = 0;

            % Update UI controls
            app.START.Enable = 'off';
            app.STOP.Enable = 'on';
            app.LoadCyclicData.Enable = 'off';
            app.OutputchannelButtonGroup.Enable = "off";
            app.LAndRButton.Enable = "off";
            app.LButton.Enable = "off";
            app.RButton.Enable = "off";

            % Start audio streaming and polling timer
            audioStreamingLoopSlim(app);
        end

        % Button pushed function: STOP
        function STOPButtonPushed(app, event)
            app.StreamingEnable = false;
            if ~isempty(app.StreamTimer) && isvalid(app.StreamTimer)
                stop(app.StreamTimer);
            end
            if ~isempty(app.AS)
                stop(app.AS);
            end
            pause(0.2);
            app.SnapshotLevelButton.Enable = "off";
            if ~app.calibrationData.initialized
                enableCalibration(app);
            end
            stopHousekeeping(app);

            if isempty(app.WorkingDirectory) || ~isfield(app.WorkingDirectory, 'fullpath')
                return;
            end

            tstr = string(datetime('now', 'Format', 'yyyyMMddHHmmss'));

            % Save initial cycle if available
            if isfield(app.InitialCondition, 'source') && ~isempty(app.InitialCondition.source) && ...
               isfield(app.InitialCondition, 'captured') && ~isempty(app.InitialCondition.captured)
                fNameSource = "initialSourceAT" + tstr + "r.wav";
                fNameCaptured = "initialCapturedAT" + tstr + "r.wav";
                sourceFullpath = fullfile(app.WorkingDirectory.fullpath, fNameSource);
                capturedFullpath = fullfile(app.WorkingDirectory.fullpath, fNameCaptured);
                audiowrite(sourceFullpath, app.InitialCondition.source, ...
                    app.SampleRate, "BitsPerSample", 24, 'Title', app.FileName);
                audiowrite(capturedFullpath, app.InitialCondition.captured, ...
                    app.SampleRate, "BitsPerSample", 24, 'Title', app.FileName);
            end

            % Save accumulated stream if data was captured
            if isfield(app.AnalysisData, 'headIdx') && app.AnalysisData.headIdx > 0
                headIdx = min(app.AnalysisData.headIdx, size(app.streamingData.wholeSource, 1));
                sourceData = app.streamingData.wholeSource(1:headIdx, :);
                capturedData = app.streamingData.capturedBuffer(1:headIdx, :);
                fNameSource = "sourceAT" + tstr + "r.wav";
                fNameCaptured = "capturedAT" + tstr + "r.wav";
                sourceFullpath = fullfile(app.WorkingDirectory.fullpath, fNameSource);
                capturedFullpath = fullfile(app.WorkingDirectory.fullpath, fNameCaptured);
                audiowrite(sourceFullpath, sourceData, ...
                    app.SampleRate, "BitsPerSample", 24, 'Title', app.FileName);
                audiowrite(capturedFullpath, capturedData, ...
                    app.SampleRate, "BitsPerSample", 24, 'Title', app.FileName);
            end
        end

        % Button pushed function: LoadCyclicData
        function LoadCyclicDataButtonPushed(app, event)
            results = SelectFileButtonPushed(app);
            if results == 0
                return;
            end
            [x, fs] = audioread(app.FilePath);
            if isempty(x)
                uialert(app.UIFigure, "Selected audio file is empty.", "File Error");
                return;
            end
            x = x(:, 1); % Enforce mono column vector for cyclic source
            app.CyclicData.name = app.FileName;
            app.CyclicData.data = x;
            app.CyclicData.samplingFrequency = fs;
            if app.SampleRate ~= fs
                results = updateSamplingRate(app);
                if results == 0
                    app.StatusLabel.Text = "This sampling rate is not supported. Select another file.";
                    return;
                end
            end
            title(app.SpectrumUIAxes, "File: " + app.FileName, "Interpreter", "none");
            app.AnalysisData.impulseLength = min(fs, floor(size(x, 1)/3));
            app.AnalysisData.dataLength = size(x, 1);
            app.AnalysisData.halfDataLength = floor(size(x, 1)/2);
            app.AnalysisData.frequencyAxis = ...
                (0:app.AnalysisData.impulseLength-1)' / app.AnalysisData.impulseLength * fs;
            app.AnalysisData.frequencyAxix = app.AnalysisData.frequencyAxis;
            frequencyAxisRow = app.AnalysisData.frequencyAxis';
            app.AnalysisData.gfStr6 = ...
                constantQgaussBPFweight(frequencyAxisRow, 1/3, 1/24);
            app.AnalysisData.gfStr12 = ...
                constantQgaussBPFweight(frequencyAxisRow, 1/12, 1/24);
            app.AnalysisData.bpfFc = 1000;
            designFcOctaveSet(app);
            app.AnalysisData.fc1kHzFilterF = app.AnalysisData.tmpBPFset;
            app.AnalysisData.fcSliderFilterF = app.AnalysisData.tmpBPFset;
            % Maintain backward-compatible field names if accessed externally
            app.AnalysisData.fc1kHzFulterF = app.AnalysisData.tmpBPFset;
            app.AnalysisData.fcSliderFulterF = app.AnalysisData.tmpBPFset;
            updateGraphics(app);
            app.StatusLabel.Text = "Ready to START";
            app.START.Enable = 'on';
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            % UIFIGURECLOSEREQUEST Safe cleanup when app window is closed.

            % Stop & delete app timer
            try
                if isprop(app, 'StreamTimer') && ~isempty(app.StreamTimer) && isvalid(app.StreamTimer)
                    stop(app.StreamTimer);
                    delete(app.StreamTimer);
                    app.StreamTimer = [];
                end
            catch ME
                disp("Error cleaning StreamTimer: " + ME.message);
            end

            % Stop and release audiostreamer if present
            try
                if isprop(app, 'AS') && ~isempty(app.AS)
                    try
                        stop(app.AS);
                    catch
                    end
                    try
                        if ismethod(app.AS, 'release')
                            release(app.AS);
                        end
                    catch
                    end
                    app.AS = [];
                end
            catch ME
                disp("Error cleaning audio stream object: " + ME.message);
            end

            % Delete capture buffers to free memory
            try
                if isprop(app, 'streamingData') && ~isempty(app.streamingData)
                    app.streamingData = [];
                end
            catch
            end

            % Finally delete the app
            try
                delete(app);
            catch ME
                try
                    if isprop(app, 'UIFigure') && isvalid(app.UIFigure)
                        delete(app.UIFigure);
                    end
                catch
                    disp("Error closing app window: " + ME.message);
                end
            end
        end

        % Button pushed function: SnapshotLevelButton
        function SnapshotLevelButtonPushed(app, event)
            xData = app.BandSPLDisplay.signalLevelHandle.XData;
            yData = app.BandSPLDisplay.signalLevelHandle.YData;
            set(app.BandSPLDisplay.snapshotLevelHandle, "XData", xData, "YData", yData);
        end

        % Button pushed function: SetworkingfolderButton
        function SetworkingfolderButtonPushed(app, event)
            selpath = uigetdir;
            if isequal(selpath, 0)
                return;
            end
            app.WorkingDirectory.base = selpath;
            tstr = string(datetime('now', 'Format', 'yyyyMMddHHmmss'));
            workingDirectoryName = "acTest" + tstr;
            app.WorkingDirectory.name = workingDirectoryName;
            fullDirPath = fullfile(selpath, workingDirectoryName);
            mkdir(fullDirPath);
            app.WorkingDirectory.fullpath = fullDirPath;
            app.SetworkingfolderButton.Enable = "off";
            app.snapshotGUIButton.Enable = "on";
            populateAudioDevices(app);
        end

        % Button pushed function: Voicememo5sButton
        function Voicememo5sButtonPushed(app, event)
            if createAudioStreamer(app) == 0
                return;
            end
            xSilence = zeros(app.SampleRate * 5, app.nChannelOut);
            record(app.AS);
            play(app.AS, xSilence, "non-blocking");
            backUpMessage = app.StatusLabel.Text;
            startTime = tic;
            while isvalid(app) && toc(startTime) <= 5.3
                pause(0.1);
                if isPlaying(app.AS)
                    app.RecordingLamp.Enable = "on";
                    app.RecordingLampLabel.Enable = "on";
                    timeLeft = max(0, round((5 - toc(startTime))*10)/10);
                    app.StatusLabel.Text = "Voice memo: " + num2str(timeLeft) + " (s) left";
                    drawnow limitrate;
                end
            end
            if ~isvalid(app)
                return;
            end
            app.StatusLabel.Text = backUpMessage;
            app.RecordingLamp.Enable = "off";
            app.RecordingLampLabel.Enable = "off";
            y = read(app.AS);
            app.VoiceMemo.data = y;
            app.VoiceMemo.SamplingRate = app.SampleRate;
            stop(app.AS);
            pause(0.2);
            release(app.AS);

            if ~isempty(y) && ~isempty(app.WorkingDirectory) && isfield(app.WorkingDirectory, 'fullpath')
                maxAmp = max(abs(y(:)));
                if maxAmp > 0
                    yNorm = (y / maxAmp) * 0.5;
                else
                    yNorm = y;
                end
                tstr = string(datetime('now', 'Format', 'yyyyMMddHHmmss'));
                fNameVoicememo = "voiceMemo" + tstr + "r.wav";
                voicememoFullpath = fullfile(app.WorkingDirectory.fullpath, fNameVoicememo);
                audiowrite(voicememoFullpath, yNorm, ...
                    app.SampleRate, "BitsPerSample", 24, 'Title', 'VoiceMemo');
            end
        end

        % Button pushed function: snapshotGUIButton
        function snapshotGUIButtonPushed(app, event)
            if isempty(app.WorkingDirectory) || ~isfield(app.WorkingDirectory, 'fullpath')
                return;
            end
            frameImage = getframe(app.UIFigure); % capture frame
            img = frame2im(frameImage);          % convert to image array
            tstr = string(datetime('now', 'Format', 'yyyyMMddHHmmss'));
            fName = "GUIsnapAT" + tstr + ".png";
            fullPathName = fullfile(app.WorkingDirectory.fullpath, fName);
            imwrite(img, fullPathName);
        end

        % Button pushed function: CalibrateButton
        function CalibrateButtonPushed(app, event)
            splLevels = [0, 95:-1:65];
            fs = app.SampleRate;
            supportedStringArray = strtrim(string(num2str(splLevels')))';
            app.MeasuredSPLDropDown.Items = supportedStringArray;
            app.MeasuredSPLDropDown.Value = app.MeasuredSPLDropDown.Items{1};

            app.MeasuredSPLDropDownLabel.Enable = "on";
            generateCalibrationSignal(app);
            if createAudioStreamer(app) == 0
                return;
            end
            record(app.AS);
            play(app.AS, app.calibrationData.calibrationSignal, "non-blocking");
            backUpMessage = app.StatusLabel.Text;
            startTime = tic;
            while isvalid(app) && isPlaying(app.AS)
                pause(0.1);
                app.RecordingLamp.Enable = "on";
                app.RecordingLampLabel.Enable = "on";
                if toc(startTime) > 3.1
                    app.MeasuredSPLDropDown.Enable = "on";
                end
                drawnow limitrate;
            end
            if ~isvalid(app)
                return;
            end
            y = read(app.AS);
            pause(0.2);
            release(app.AS);

            nSamples = size(y, 1);
            tailLength = min(nSamples, 3 * fs);
            if tailLength > 0
                tailSegment = y(end - tailLength + 1 : end, 1);
                observedSTDdB = 20*log10(max(std(tailSegment), 1e-12));
                app.CALValuedBEditField.Value = ...
                    app.calibrationData.ch1Measured - observedSTDdB;
                app.calibrationData.ch1CalibrationDb = ...
                    app.calibrationData.ch1Measured - observedSTDdB;
            end
            app.RecordingLamp.Enable = "off";
            app.RecordingLampLabel.Enable = "off";
            app.CalibrateButton.Enable = "off";
            app.MeasuredSPLDropDown.Enable = "off";
            app.MeasuredSPLDropDownLabel.Enable = "off";
            app.StatusLabel.Text = backUpMessage;
            drawnow;
        end

        % Value changed function: MeasuredSPLDropDown
        function MeasuredSPLDropDownValueChanged(app, event)
            value = app.MeasuredSPLDropDown.Value;
            app.calibrationData.ch1Measured = str2double(value);
            app.BandSPLAxes.YLim = [0 100];
            app.BandSPLAxes.YLabel.String = "CAL. SPL (dB)";
            if ~isempty(app.AS)
                stop(app.AS);
            end
        end

        % Value changed function: BPFfcSlider
        function BPFfcSliderValueChanged(app, event)
            fcSet = app.AnalysisData.gfStr12.fcSet;
            idx = min(max(round(app.BPFfcSlider.Value), 1), numel(fcSet));
            app.BPFfcSlider.Value = idx;
            fc = fcSet(idx);
            app.AnalysisData.bpfFc = fc;
            app.SpectrumDisplay.fc1kHandle.Value = fc;
            app.SpectrumDisplay.fc1kHandle.Label = num2str(fc, "%6.0f");
            app.BandSPLDisplay.fc1kHandle.Value = fc;
            app.BandSPLDisplay.fc1kHandle.Label = num2str(fc, "%6.0f");
            designFcOctaveSet(app);
            app.AnalysisData.fcSliderFilterF = app.AnalysisData.tmpBPFset;
            app.AnalysisData.fcSliderFulterF = app.AnalysisData.tmpBPFset;
        end

        % Value changing function: BPFfcSlider
        function BPFfcSliderValueChanging(app, event)
            fcSet = app.AnalysisData.gfStr12.fcSet;
            idx = min(max(round(event.Value), 1), numel(fcSet));
            fc = fcSet(idx);
            app.BPFfcSliderLabel.Text = num2str(fc, "%7.0f");
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {665, 665, 665};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {665, 665};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {289, '1x', 658};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1121 665];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {289, '1x', 658};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create IOdeviceselectionLabel
            app.IOdeviceselectionLabel = uilabel(app.LeftPanel);
            app.IOdeviceselectionLabel.Enable = 'off';
            app.IOdeviceselectionLabel.Position = [30 515 112 22];
            app.IOdeviceselectionLabel.Text = 'IO device selection';

            % Create IOdeviceselectionDropDown
            app.IOdeviceselectionDropDown = uidropdown(app.LeftPanel);
            app.IOdeviceselectionDropDown.ValueChangedFcn = createCallbackFcn(app, @IOdeviceselectionDropDownValueChanged, true);
            app.IOdeviceselectionDropDown.Enable = 'off';
            app.IOdeviceselectionDropDown.Position = [37 494 225 22];

            % Create StatusLabel
            app.StatusLabel = uilabel(app.LeftPanel);
            app.StatusLabel.Position = [32 544 229 68];
            app.StatusLabel.Text = 'StatusLabel';

            % Create SystemStatusLabel
            app.SystemStatusLabel = uilabel(app.LeftPanel);
            app.SystemStatusLabel.Position = [32 611 83 22];
            app.SystemStatusLabel.Text = 'System Status';

            % Create SamplingfrequencymenuLabel
            app.SamplingfrequencymenuLabel = uilabel(app.LeftPanel);
            app.SamplingfrequencymenuLabel.Enable = 'off';
            app.SamplingfrequencymenuLabel.Position = [30 460 112 22];
            app.SamplingfrequencymenuLabel.Text = 'Sampling frequency';

            % Create SamplingFrequencyDropDown
            app.SamplingFrequencyDropDown = uidropdown(app.LeftPanel);
            app.SamplingFrequencyDropDown.ValueChangedFcn = createCallbackFcn(app, @SamplingFrequencyDropDownValueChanged, true);
            app.SamplingFrequencyDropDown.Enable = 'off';
            app.SamplingFrequencyDropDown.Position = [37 439 225 22];

            % Create SamplingfrequencymenuLabel_2
            app.SamplingfrequencymenuLabel_2 = uilabel(app.LeftPanel);
            app.SamplingfrequencymenuLabel_2.Enable = 'off';
            app.SamplingfrequencymenuLabel_2.Position = [30 406 86 22];
            app.SamplingfrequencymenuLabel_2.Text = 'Input Channels';

            % Create InputChannelDropDown
            app.InputChannelDropDown = uidropdown(app.LeftPanel);
            app.InputChannelDropDown.Items = {'1', '2', '0'};
            app.InputChannelDropDown.ItemsData = {'1', '2', '0'};
            app.InputChannelDropDown.ValueChangedFcn = createCallbackFcn(app, @InputChannelDropDownValueChanged, true);
            app.InputChannelDropDown.Enable = 'off';
            app.InputChannelDropDown.Position = [37 385 225 22];
            app.InputChannelDropDown.Value = '0';

            % Create START
            app.START = uibutton(app.LeftPanel, 'push');
            app.START.ButtonPushedFcn = createCallbackFcn(app, @STARTButtonPushed, true);
            app.START.Enable = 'off';
            app.START.Position = [47 123 196 23];
            app.START.Text = 'START';

            % Create STOP
            app.STOP = uibutton(app.LeftPanel, 'push');
            app.STOP.ButtonPushedFcn = createCallbackFcn(app, @STOPButtonPushed, true);
            app.STOP.Enable = 'off';
            app.STOP.Position = [48 81 196 23];
            app.STOP.Text = 'STOP';

            % Create LoadCyclicData
            app.LoadCyclicData = uibutton(app.LeftPanel, 'push');
            app.LoadCyclicData.ButtonPushedFcn = createCallbackFcn(app, @LoadCyclicDataButtonPushed, true);
            app.LoadCyclicData.Enable = 'off';
            app.LoadCyclicData.Position = [49 340 196 23];
            app.LoadCyclicData.Text = 'Load cyclic data';

            % Create SpectrumdisplayButtonGroup
            app.SpectrumdisplayButtonGroup = uibuttongroup(app.LeftPanel);
            app.SpectrumdisplayButtonGroup.Enable = 'off';
            app.SpectrumdisplayButtonGroup.Title = 'Spectrum display';
            app.SpectrumdisplayButtonGroup.Position = [37 165 224 106];

            % Create RawDFTButton
            app.RawDFTButton = uiradiobutton(app.SpectrumdisplayButtonGroup);
            app.RawDFTButton.Enable = 'off';
            app.RawDFTButton.Text = 'Raw DFT ';
            app.RawDFTButton.Position = [11 60 195 22];
            app.RawDFTButton.Value = true;

            % Create Octave6thRadioButton
            app.Octave6thRadioButton = uiradiobutton(app.SpectrumdisplayButtonGroup);
            app.Octave6thRadioButton.Enable = 'off';
            app.Octave6thRadioButton.Text = '1/3 octave smoothed';
            app.Octave6thRadioButton.Position = [11 38 184 22];

            % Create Octave12thRadioButton
            app.Octave12thRadioButton = uiradiobutton(app.SpectrumdisplayButtonGroup);
            app.Octave12thRadioButton.Enable = 'off';
            app.Octave12thRadioButton.Text = '1/12 octave smoothed';
            app.Octave12thRadioButton.Position = [11 16 184 22];

            % Create OutputchannelButtonGroup
            app.OutputchannelButtonGroup = uibuttongroup(app.LeftPanel);
            app.OutputchannelButtonGroup.Enable = 'off';
            app.OutputchannelButtonGroup.Title = 'Output channel';
            app.OutputchannelButtonGroup.Position = [37 281 224 50];

            % Create LButton
            app.LButton = uiradiobutton(app.OutputchannelButtonGroup);
            app.LButton.Enable = 'off';
            app.LButton.Text = 'L';
            app.LButton.Position = [11 4 29 22];
            app.LButton.Value = true;

            % Create RButton
            app.RButton = uiradiobutton(app.OutputchannelButtonGroup);
            app.RButton.Enable = 'off';
            app.RButton.Text = 'R';
            app.RButton.Position = [75 4 30 22];

            % Create LAndRButton
            app.LAndRButton = uiradiobutton(app.OutputchannelButtonGroup);
            app.LAndRButton.Enable = 'off';
            app.LAndRButton.Text = 'L and R';
            app.LAndRButton.Position = [139 4 64 22];

            % Create snapshotGUIButton
            app.snapshotGUIButton = uibutton(app.LeftPanel, 'push');
            app.snapshotGUIButton.ButtonPushedFcn = createCallbackFcn(app, @snapshotGUIButtonPushed, true);
            app.snapshotGUIButton.Enable = 'off';
            app.snapshotGUIButton.Position = [47 15 100 23];
            app.snapshotGUIButton.Text = 'snapshot GUI';

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create SetworkingfolderButton
            app.SetworkingfolderButton = uibutton(app.CenterPanel, 'push');
            app.SetworkingfolderButton.ButtonPushedFcn = createCallbackFcn(app, @SetworkingfolderButtonPushed, true);
            app.SetworkingfolderButton.Position = [19 610 137 23];
            app.SetworkingfolderButton.Text = 'Set working folder';

            % Create Voicememo5sButton
            app.Voicememo5sButton = uibutton(app.CenterPanel, 'push');
            app.Voicememo5sButton.ButtonPushedFcn = createCallbackFcn(app, @Voicememo5sButtonPushed, true);
            app.Voicememo5sButton.Enable = 'off';
            app.Voicememo5sButton.Position = [19 543 137 23];
            app.Voicememo5sButton.Text = 'Voice memo 5s';

            % Create CalibrationCh1Panel
            app.CalibrationCh1Panel = uipanel(app.CenterPanel);
            app.CalibrationCh1Panel.Enable = 'off';
            app.CalibrationCh1Panel.Title = 'Calibration Ch-1';
            app.CalibrationCh1Panel.Position = [13 336 152 165];

            % Create CalibrateButton
            app.CalibrateButton = uibutton(app.CalibrationCh1Panel, 'push');
            app.CalibrateButton.ButtonPushedFcn = createCallbackFcn(app, @CalibrateButtonPushed, true);
            app.CalibrateButton.Enable = 'off';
            app.CalibrateButton.Position = [21 108 110 23];
            app.CalibrateButton.Text = 'Calibrate';

            % Create MeasuredSPLDropDownLabel
            app.MeasuredSPLDropDownLabel = uilabel(app.CalibrationCh1Panel);
            app.MeasuredSPLDropDownLabel.HorizontalAlignment = 'right';
            app.MeasuredSPLDropDownLabel.Enable = 'off';
            app.MeasuredSPLDropDownLabel.Position = [34 83 84 22];
            app.MeasuredSPLDropDownLabel.Text = 'Measured SPL';

            % Create MeasuredSPLDropDown
            app.MeasuredSPLDropDown = uidropdown(app.CalibrationCh1Panel);
            app.MeasuredSPLDropDown.ValueChangedFcn = createCallbackFcn(app, @MeasuredSPLDropDownValueChanged, true);
            app.MeasuredSPLDropDown.Enable = 'off';
            app.MeasuredSPLDropDown.Position = [21 58 110 22];

            % Create CALValuedBEditFieldLabel
            app.CALValuedBEditFieldLabel = uilabel(app.CalibrationCh1Panel);
            app.CALValuedBEditFieldLabel.HorizontalAlignment = 'right';
            app.CALValuedBEditFieldLabel.Enable = 'off';
            app.CALValuedBEditFieldLabel.Position = [33 34 86 22];
            app.CALValuedBEditFieldLabel.Text = 'CAL. Value: dB';

            % Create CALValuedBEditField
            app.CALValuedBEditField = uieditfield(app.CalibrationCh1Panel, 'numeric');
            app.CALValuedBEditField.ValueDisplayFormat = '%9.1f';
            app.CALValuedBEditField.Enable = 'off';
            app.CALValuedBEditField.Position = [35 10 81 22];

            % Create CalibrationCh2Panel
            app.CalibrationCh2Panel = uipanel(app.CenterPanel);
            app.CalibrationCh2Panel.Enable = 'off';
            app.CalibrationCh2Panel.Title = 'Calibration Ch-2';
            app.CalibrationCh2Panel.Position = [13 157 152 165];

            % Create CalibrateButton_2
            app.CalibrateButton_2 = uibutton(app.CalibrationCh2Panel, 'push');
            app.CalibrateButton_2.Enable = 'off';
            app.CalibrateButton_2.Position = [21 108 110 23];
            app.CalibrateButton_2.Text = 'Calibrate';

            % Create MeasuredSPLDropDown_2Label
            app.MeasuredSPLDropDown_2Label = uilabel(app.CalibrationCh2Panel);
            app.MeasuredSPLDropDown_2Label.HorizontalAlignment = 'right';
            app.MeasuredSPLDropDown_2Label.Enable = 'off';
            app.MeasuredSPLDropDown_2Label.Position = [34 83 84 22];
            app.MeasuredSPLDropDown_2Label.Text = 'Measured SPL';

            % Create MeasuredSPLDropDown_2
            app.MeasuredSPLDropDown_2 = uidropdown(app.CalibrationCh2Panel);
            app.MeasuredSPLDropDown_2.Enable = 'off';
            app.MeasuredSPLDropDown_2.Position = [21 58 110 22];

            % Create CALValuedBEditField_2Label
            app.CALValuedBEditField_2Label = uilabel(app.CalibrationCh2Panel);
            app.CALValuedBEditField_2Label.HorizontalAlignment = 'right';
            app.CALValuedBEditField_2Label.Enable = 'off';
            app.CALValuedBEditField_2Label.Position = [33 34 86 22];
            app.CALValuedBEditField_2Label.Text = 'CAL. Value: dB';

            % Create CALValuedBEditField_2
            app.CALValuedBEditField_2 = uieditfield(app.CalibrationCh2Panel, 'numeric');
            app.CALValuedBEditField_2.ValueDisplayFormat = '%9.1f';
            app.CALValuedBEditField_2.Enable = 'off';
            app.CALValuedBEditField_2.Position = [35 10 81 22];

            % Create RecordingLampLabel
            app.RecordingLampLabel = uilabel(app.CenterPanel);
            app.RecordingLampLabel.HorizontalAlignment = 'right';
            app.RecordingLampLabel.Enable = 'off';
            app.RecordingLampLabel.Position = [37 575 60 22];
            app.RecordingLampLabel.Text = 'Recording';

            % Create RecordingLamp
            app.RecordingLamp = uilamp(app.CenterPanel);
            app.RecordingLamp.Enable = 'off';
            app.RecordingLamp.Position = [112 575 20 20];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Create SpectrumUIAxes
            app.SpectrumUIAxes = uiaxes(app.RightPanel);
            title(app.SpectrumUIAxes, 'Title')
            xlabel(app.SpectrumUIAxes, 'X')
            ylabel(app.SpectrumUIAxes, 'Y')
            zlabel(app.SpectrumUIAxes, 'Z')
            app.SpectrumUIAxes.Box = 'on';
            app.SpectrumUIAxes.Position = [32 361 593 296];

            % Create BandSPLAxes
            app.BandSPLAxes = uiaxes(app.RightPanel);
            title(app.BandSPLAxes, 'Title')
            xlabel(app.BandSPLAxes, 'X')
            ylabel(app.BandSPLAxes, 'Y')
            zlabel(app.BandSPLAxes, 'Z')
            app.BandSPLAxes.Box = 'on';
            app.BandSPLAxes.Position = [32 42 316 265];

            % Create bpfPowerAxes
            app.bpfPowerAxes = uiaxes(app.RightPanel);
            title(app.bpfPowerAxes, 'Title')
            xlabel(app.bpfPowerAxes, 'X')
            ylabel(app.bpfPowerAxes, 'Y')
            zlabel(app.bpfPowerAxes, 'Z')
            app.bpfPowerAxes.Box = 'on';
            app.bpfPowerAxes.Position = [347 42 300 265];

            % Create SnapshotLevelButton
            app.SnapshotLevelButton = uibutton(app.RightPanel, 'push');
            app.SnapshotLevelButton.ButtonPushedFcn = createCallbackFcn(app, @SnapshotLevelButtonPushed, true);
            app.SnapshotLevelButton.Enable = 'off';
            app.SnapshotLevelButton.Position = [32 15 100 23];
            app.SnapshotLevelButton.Text = 'Snapshot Level';

            % Create BPFfcSliderLabel
            app.BPFfcSliderLabel = uilabel(app.RightPanel);
            app.BPFfcSliderLabel.HorizontalAlignment = 'right';
            app.BPFfcSliderLabel.Position = [6 322 48 22];
            app.BPFfcSliderLabel.Text = 'BPF fc';

            % Create BPFfcSlider
            app.BPFfcSlider = uislider(app.RightPanel);
            app.BPFfcSlider.MajorTicks = [];
            app.BPFfcSlider.MajorTickLabels = {''};
            app.BPFfcSlider.ValueChangedFcn = createCallbackFcn(app, @BPFfcSliderValueChanged, true);
            app.BPFfcSlider.ValueChangingFcn = createCallbackFcn(app, @BPFfcSliderValueChanging, true);
            app.BPFfcSlider.MinorTicks = 0:100;
            app.BPFfcSlider.Step = 1;
            app.BPFfcSlider.Enable = 'off';
            app.BPFfcSlider.Position = [70 326 547 3];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = acouticCheckerByGFFTRefactored_v2_exportedR

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end