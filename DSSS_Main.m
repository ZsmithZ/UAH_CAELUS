% Main script to call function and analyze waveform
clc; clear; close all;

% Simulation parameters
fs = 10e9; % Sampling frequency
num_bits = 1000;
data = randi([0 1], 1, num_bits); % Random binary data

% Generate DSSS waveform
[t, dsss_signal, dsss_data_upsampled] = generate_dsss_waveform(data, fs);

% Add dithering (adjust SNR)
snr = 30; % Set desired SNR in dB
dsss_signal_noisy = awgn(dsss_signal, snr, 'measured');

% Plot input data
figure;
stem(data(1:100), 'filled');
xlabel('Bit Index');
ylabel('Bit Value');
title('Input Data');
grid on;

% Plot power spectral density (PSD)
figure;
hold on;
pwelch(dsss_signal, [], [], [], fs, 'centered');
pwelch(dsss_signal_noisy, [], [], [], fs, 'centered');
hold off;
title('Power Spectral Density of DSSS Signal');
legend('Original Signal', 'Noisy Signal');

% Plot I/Q constellation
iq_data = dsss_data_upsampled + 1j * hilbert(dsss_data_upsampled);
figure;
hold on;
scatter(real(iq_data(1:1000)), imag(iq_data(1:1000)), '.');
iq_data_noisy = dsss_signal_noisy(1:1000) + 1j * hilbert(dsss_signal_noisy(1:1000));
scatter(real(iq_data_noisy), imag(iq_data_noisy), 'r.');
hold off;
xlabel('In-phase');
ylabel('Quadrature');
title('I/Q Constellation');
legend('Original', 'Noisy');
grid on;
