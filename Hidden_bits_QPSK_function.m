% MATLAB Simulation of Hidden Bit Encoding in QPSK
% Simulates three methods: Phase Perturbation, Power Level Encoding, and Differential Phase Shift

clear; clc; close all;

%% Simulation Parameters
N = 1e5; % Number of symbols
EbN0_dB = 0:2:20; % SNR range in dB
M = 4; % QPSK

% Phase Perturbation Parameters
delta_theta = pi/36; % Small phase perturbation

% Power Level Encoding Parameters
delta_A = 0.05; % Small amplitude variation

% Differential Phase Shift Parameters
delta_phi = pi/36;

% Adjustable Parameters
operating_frequency = 2.4e9; % Operating frequency in Hz (e.g., 2.4 GHz)
channel_bandwidth = 20e6; % Channel bandwidth in Hz (e.g., 20 MHz)
FEC_coding_rate = 3/4; % Forward Error Correction rate

% Compute Spectral Efficiency
modulation_order = log2(M); % Bits per symbol
spectral_efficiency = (modulation_order * FEC_coding_rate) / (channel_bandwidth/1e6); % bps/Hz

BER_standard = zeros(size(EbN0_dB));
BER_hidden_phase = zeros(size(EbN0_dB));
BER_hidden_power = zeros(size(EbN0_dB));
BER_hidden_diff = zeros(size(EbN0_dB));

% Function for BER Calculation
function [BER_standard, BER_hidden_phase, BER_hidden_power, BER_hidden_diff] = compute_BER(EbN0_dB, N, M, delta_theta, delta_A, delta_phi)
    BER_standard = zeros(size(EbN0_dB));
    BER_hidden_phase = zeros(size(EbN0_dB));
    BER_hidden_power = zeros(size(EbN0_dB));
    BER_hidden_diff = zeros(size(EbN0_dB));
    
    for i = 1:length(EbN0_dB)
        EbN0 = 10^(EbN0_dB(i)/10);
        noise_var = 1/(2*EbN0);
        
        %% Standard QPSK Modulation
        bits = randi([0 1], N*2, 1);
        symbols = (1-2*bits(1:2:end)) + 1j*(1-2*bits(2:2:end));
        symbols = symbols / sqrt(2);
        
        %% Phase Perturbation Encoding
        hidden_bits_phase = randi([0 1], N, 1);
        perturbation = (hidden_bits_phase * 2 - 1) * delta_theta;
        symbols_phase = symbols .* exp(1j * perturbation);
        
        %% Power Level Encoding
        hidden_bits_power = randi([0 1], N, 1);
        amplitude_variation = 1 + (hidden_bits_power * 2 - 1) * delta_A;
        symbols_power = symbols .* amplitude_variation;
        
        %% Differential Phase Shift Encoding
        hidden_bits_diff = randi([0 1], N, 1);
        diff_shift = (hidden_bits_diff * 2 - 1) * delta_phi;
        symbols_diff = symbols .* exp(1j * diff_shift);
        
        %% Add Noise
        noise = sqrt(noise_var) * (randn(N, 1) + 1j * randn(N, 1));
        received_standard = symbols + noise;
        received_phase = symbols_phase + noise;
        received_power = symbols_power + noise;
        received_diff = symbols_diff + noise;
    end
end

[BER_standard, BER_hidden_phase, BER_hidden_power, BER_hidden_diff] = compute_BER(EbN0_dB, N, M, delta_theta, delta_A, delta_phi);
