function [t, dsss_signal, dsss_data_upsampled] = generate_dsss_waveform(data, fs)
    % Parameters
    carrier_freq = 9.25e9;  % Carrier frequency in Hz
    data_rate = 75e6;       % Data rate in bits per second
    spreading_factor = 8;   % Spreading factor for DSSS
    
    % Generate PN sequence (length spreading_factor)
    pn_sequence = randi([0 1], 1, spreading_factor) * 2 - 1; % BPSK (-1, 1)
    
    % Oversampling factor based on fs and data rate
    samples_per_bit = round(fs / data_rate);
    samples_per_chip = round(samples_per_bit / spreading_factor);
    
    % Expand data to match spreading sequence
    spread_data = repelem(data * 2 - 1, spreading_factor); % BPSK mapping
    dsss_data = kron(spread_data, pn_sequence); % Spreading operation
    
    % Upsample to match the sampling frequency
    dsss_data_upsampled = repelem(dsss_data, max(1, samples_per_chip));
    
    % Time vector
    t = (0:length(dsss_data_upsampled)-1) / fs;
    
    % Carrier modulation
    carrier = cos(2 * pi * carrier_freq * t); % Carrier waveform
    dsss_signal = dsss_data_upsampled .* carrier; % Modulation
end