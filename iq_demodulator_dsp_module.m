function iq_live_input_demodulator_nodsp()
    %% Load IQ Input
    try
        load('iq_input.mat', 'rxIQ', 'fs');  % Expected complex baseband and sample rate
    catch
        error('Expected input file "iq_input.mat" with variables rxIQ and fs');
    end

    if ~isvector(rxIQ) || ~isnumeric(fs)
        error('rxIQ must be complex vector, fs must be scalar.');
    end

    %% System Parameters
    M = 4;              % QPSK
    Rs = 100e3;         % Symbol rate
    sps = fs / Rs;
    k = log2(M);        % Bits per symbol

    %% Matched Filtering and Downsampling
    rrc = rcosdesign(0.35, 4, sps, 'normal');
    rxFiltered = filter(rrc, 1, rxIQ);

    offset = round(sps/2);
    rxSymbols = downsample(rxFiltered, sps, offset);

    %% Demodulation
    demodData = qamdemod(rxSymbols, M, 'UnitAveragePower', true);
    rxBits = de2bi(demodData, k, 'left-msb')';
    rxBits = rxBits(:);

    %% Message Recovery
    try
        rxBytes = reshape(rxBits(1:floor(length(rxBits)/8)*8), 8, []).';
        rxMsg = char(bi2de(rxBytes, 'left-msb'))';
    catch
        rxMsg = '[Decoding Error]';
    end

    %% Visualization
    figure('Name','IQ Demodulation Results','NumberTitle','off');

    % Spectrum Plot
    subplot(2,1,1);
    N = length(rxIQ);
    f = linspace(0, fs/2, N/2);
    P = abs(fft(rxIQ));
    plot(f/1e6, 20*log10(P(1:N/2)/max(P)));
    title('Input Spectrum'); xlabel('Frequency (MHz)'); ylabel('Power (dB)');
    grid on;

    % Constellation
    subplot(2,1,2);
    plot(real(rxSymbols), imag(rxSymbols), 'o');
    title('I/Q Constellation'); xlabel('In-Phase'); ylabel('Quadrature');
    axis equal; grid on;

    %% Message Output
    disp(['[Decoded Message]: ', rxMsg]);
end
