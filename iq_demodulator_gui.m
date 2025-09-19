function iq_demodulator_gui()
    %%% System Parameters
    fs = 1e6;                % Sample rate
    fc = 100e3;              % Carrier frequency
    Rs = 100e3;              % Symbol rate
    M = 4;                   % QPSK
    k = log2(M);
    sps = fs / Rs;
    msg = 'HELLO MATLAB GUI';

    bits = de2bi(double(msg), 8, 'left-msb')';
    bits = bits(:);
    dataSym = bi2de(reshape(bits, k, []).', 'left-msb');
    modSig = qammod(dataSym, M, 'UnitAveragePower', true);
    txSig = upfirdn(modSig, rcosdesign(0.35, 4, sps), sps);
    t = (0:length(txSig)-1)' / fs;
    passbandSig = real(txSig .* exp(1j*2*pi*fc*t));

    % UI Initialization
    f = figure('Name', 'Real-Time I/Q Demodulator', 'NumberTitle', 'off', 'Position', [100 100 1200 700]);

    % Sliders for SNR, CFO, and Phase
    uicontrol('Style', 'text', 'String', 'SNR (dB)', 'Position', [50 650 100 20]);
    snrSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 30, 'Value', 20, 'Position', [50 630 120 20]);

    uicontrol('Style', 'text', 'String', 'CFO (Hz)', 'Position', [200 650 100 20]);
    cfoSlider = uicontrol('Style', 'slider', 'Min', -1e4, 'Max', 1e4, 'Value', 0, 'Position', [200 630 120 20]);

    uicontrol('Style', 'text', 'String', 'Phase Offset (°)', 'Position', [350 650 120 20]);
    phaseSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 180, 'Value', 0, 'Position', [350 630 120 20]);

    % Axes for Spectrum, Constellation, and Message
    ax1 = subplot(3,1,1); title('Spectrum Analyzer'); colormap jet
    ax2 = subplot(3,1,2); title('I/Q Constellation');
    ax3 = subplot(3,1,3); axis off;
    textBox = text(0.1, 0.5, '', 'FontSize', 14, 'FontWeight', 'bold');

    % Processing Loop
    while ishandle(f)
        % Slider Values
        snr = snrSlider.Value;
        cfo = cfoSlider.Value;
        phaseOffset = deg2rad(phaseSlider.Value);

        % Channel Impairments
        impairedSig = awgn(passbandSig, snr, 'measured');
        impairedSig = impairedSig .* exp(1j*(2*pi*cfo*t + phaseOffset));

        % I/Q Demodulation
        I = real(impairedSig) .* cos(2*pi*fc*t);
        Q = -real(impairedSig) .* sin(2*pi*fc*t);
        baseband = I + 1j*Q;

        lpf = designfilt('lowpassfir', 'PassbandFrequency', 0.45*Rs, ...
                         'StopbandFrequency', 0.55*Rs, 'SampleRate', fs);
        baseband = filter(lpf, baseband);

        rxSym = downsample(baseband, sps, round(sps/2));
        rxSym = rxSym .* exp(-1j*phaseOffset); % simple phase correction

        rxDataSym = qamdemod(rxSym, M, 'UnitAveragePower', true);
        rxBits = de2bi(rxDataSym, k, 'left-msb')';
        rxBits = rxBits(:);

        rxBytes = reshape(rxBits(1:floor(length(rxBits)/8)*8), 8, []).';
        try
            rxMsg = char(bi2de(rxBytes, 'left-msb'))';
        catch
            rxMsg = '[Decoding Error]';
        end

        % Update Plots
        subplot(ax1);
        spectrogram(impairedSig, 256, 200, 256, fs, 'yaxis');
        title(['Spectrum Analyzer | SNR = ', num2str(snr), ' dB, CFO = ', num2str(cfo), ' Hz']);

        subplot(ax2);
        plot(rxSym, 'o');
        title('I/Q Constellation');
        xlabel('In-Phase'); ylabel('Quadrature'); axis equal; grid on;

        subplot(ax3);
        set(textBox, 'String', ['Decoded Message: ', rxMsg]);

        drawnow;
    end
end
