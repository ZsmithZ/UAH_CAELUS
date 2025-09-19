% =========================================================================
% DigiSim Main Function
% =========================================================================
% Author:           Aaron Smith
% Organization:     Northrop Grumman PGS Antenna Design & Analysis
% Created:          12-31-2024
% Last Updated:     01-17-2025
% MATLAB Version:   24.1.0.2568132 (R2024a) Update 1
% 
% Description:
% This script performs [briefly describe the script's purpose and functionality].
% It includes data preprocessing, analysis, and result visualization steps.
%
% Inputs:
% - [Input Description, e.g., Input data files or parameters]
%
% Outputs:
% - [Output Description, e.g., Visualizations, processed data, etc.]
%
% Dependencies:
% - [List of required functions, toolboxes, or external files]
% =========================================================================

% process_signal_with_pfb.m
% Processes a 2.5 GHz signal through a polyphase filter bank and compares two responses.

clear; clc; close all;

%% Parameters
fs = 5e9;                  % Sampling frequency (5 GHz)
signal_freq = 2.5e9;       % Signal frequency (2.5 GHz)
num_channels = 16;         % Number of polyphase filter bank channels
filter_order = 128;        % Filter order for the prototype filter
decimation_factor = fs / num_channels;

% Generate time vector and signal
t = 0:1/fs:1e-6;           % Time vector (1 microsecond duration)
input_signal = cos(2*pi*signal_freq*t);  % Simulated 2.5 GHz cosine signal

% Ensure the input signal length is a multiple of num_channels
num_samples = floor(length(input_signal) / num_channels) * num_channels;
input_signal = input_signal(1:num_samples);

%% Polyphase Filter Bank Design
% Generate a prototype lowpass filter
prototype_filter = fir1(filter_order, 1/num_channels);

% Ensure the filter length is a multiple of the number of channels
filter_length = length(prototype_filter);
remainder = mod(filter_length, num_channels);
if remainder ~= 0
    padding = num_channels - remainder;
    prototype_filter = [prototype_filter, zeros(1, padding)];
end

% Create polyphase components
P = reshape(prototype_filter, num_channels, []);  % Polyphase matrix

%% Apply Polyphase Filter Bank
% Initialize channel outputs
filtered_signals = zeros(num_channels, num_samples / num_channels);

% Process signal through each polyphase branch
for k = 1:num_channels
    filtered_signal = filter(P(k, :), 1, input_signal);  % Filter input signal
    filtered_signals(k, :) = downsample(filtered_signal, num_channels);  % Downsample
end

%% Extract Two Responses
channel_1 = 3;  % Example: Extract channel 3
channel_2 = 7;  % Example: Extract channel 7
response_1 = filtered_signals(channel_1, :);
response_2 = filtered_signals(channel_2, :);

%% Plot Results
% Plot time-domain comparison
figure;
subplot(1, 2, 1);
plot(response_1, 'LineWidth', 1.5);
title(sprintf('Channel %d Response', channel_1));
xlabel('Samples');
ylabel('Amplitude');
grid on;

subplot(1, 2, 2);
plot(response_2, 'LineWidth', 1.5, 'Color', 'r');
title(sprintf('Channel %d Response', channel_2));
xlabel('Samples');
ylabel('Amplitude');
grid on;

% Plot frequency-domain comparison
figure;
subplot(1, 2, 1);
freq_response_1 = fftshift(abs(fft(response_1)));
plot(linspace(-fs/(2*num_channels), fs/(2*num_channels), length(freq_response_1)), freq_response_1, 'LineWidth', 1.5);
title(sprintf('Channel %d Frequency Response', channel_1));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

subplot(1, 2, 2);
freq_response_2
