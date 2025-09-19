%% Live Script for DSP Analysis of QPSK Signal
% This script utilizes the DSP toolbox to generate, process, and analyze a
% QPSK signal with real-time PSD visualization and input/output waveform comparison.

clear; clc; close all;

%% Parameters
N = 1e5; % Number of symbols
M = 4; % QPSK modulation order
EbN0_dB = 10; % SNR in dB
Fs = 1e6; % Sampling frequency
Ts = 1/Fs; % Sampling period

% Generate random bits and map to QPSK symbols
bits = randi([0 1], N*2, 1); 
symbols = (1-2*bits(1:2:end)) + 1j*(1-2*bits(2:2:end));
symbols = symbols / sqrt(2);

%% Add Noise
EbN0 = 10^(EbN0_dB/10);
noise_var = 1/(2*EbN0);
noise = sqrt(noise_var) * (randn(N, 1) + 1j * randn(N, 1));
received_signal = symbols + noise;

%% Create DSP System Objects
spectrumAnalyzer = dsp.SpectrumAnalyzer('SampleRate', Fs, 'ShowLegend', true, ...
    'Title', 'Power Spectral Density', 'YLimits', [-100 10]);

timeScope = dsp.TimeScope('SampleRate', Fs, 'TimeSpanSource', 'Property', ...
    'TimeSpan', N*Ts, 'YLimits', [-2 2], 'ShowLegend', true, 'Title', 'Input vs Output Waveform');

%% Run Live Analysis
for i = 1:1000:N
    idx = i:min(i+999, N);
    spectrumAnalyzer(received_signal(idx));
    timeScope([real(symbols(idx)), real(received_signal(idx))]);
    pause(0.1); % Small pause for visualization updates
end

%% Release System Objects
release(spectrumAnalyzer);
release(timeScope);
