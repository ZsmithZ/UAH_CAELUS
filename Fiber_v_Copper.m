%% Space-Based Payload Trade Study: 10G Ethernet vs. Optical Fiber
% Author: [Your Name]
% Date: [Current Date]
% Description: Performance analysis and visualization of key metrics

clear; clc; close all;

%% Define Parameters
distance = linspace(10, 1000, 100); % Distance in meters
bandwidth_eth = 10; % 10G Ethernet max bandwidth (Gbps)
bandwidth_fiber = 100; % Optical fiber WDM bandwidth (Gbps)

latency_eth_base = 1.5e-6; % Ethernet base latency (seconds per meter)
latency_fiber_base = 0.005e-6; % Optical fiber base latency (seconds per meter)

power_eth = 0.8 * distance; % Power consumption (W)
power_fiber = 0.2 * distance; % Optical interconnects consume less power

radiation_dose = linspace(1, 500, 100); % Total Ionizing Dose (krad)
ber_eth = 1e-9 * (1 + 0.001 * radiation_dose); % Ethernet BER increases with radiation
ber_fiber = 1e-12 * (1 + 0.0001 * radiation_dose); % Fiber is more resistant

%% Compute Metrics
throughput_eth = bandwidth_eth * (1 - ber_eth); % Adjusted throughput
throughput_fiber = bandwidth_fiber * (1 - ber_fiber);

latency_eth = latency_eth_base * distance; % Ethernet latency
latency_fiber = latency_fiber_base * distance; % Optical latency

jitter_eth = 0.1e-6 + 0.05e-6 * randn(size(distance)); % Ethernet jitter varies
jitter_fiber = 0.005e-6 + 0.002e-6 * randn(size(distance)); % Fiber has lower jitter

%% Plot Results
figure;
subplot(2,2,1);
plot(distance, throughput_eth, 'r', 'LineWidth', 2);
hold on;
plot(distance, throughput_fiber, 'b', 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Throughput (Gbps)');
title('Throughput vs. Distance');
legend('10G Ethernet', 'Optical Fiber');

subplot(2,2,2);
plot(distance, latency_eth * 1e6, 'r', 'LineWidth', 2);
hold on;
plot(distance, latency_fiber * 1e6, 'b', 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Latency (\mus)');
title('Latency vs. Distance');
legend('10G Ethernet', 'Optical Fiber');

subplot(2,2,3);
plot(radiation_dose, ber_eth, 'r', 'LineWidth', 2);
hold on;
plot(radiation_dose, ber_fiber, 'b', 'LineWidth', 2);
xlabel('Radiation Dose (krad)'); ylabel('Bit Error Rate (BER)');
title('BER vs. Radiation Exposure');
legend('10G Ethernet', 'Optical Fiber');

subplot(2,2,4);
plot(distance, power_eth, 'r', 'LineWidth', 2);
hold on;
plot(distance, power_fiber, 'b', 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Power Consumption (W)');
title('Power Consumption vs. Distance');
legend('10G Ethernet', 'Optical Fiber');

figure;
plot(distance, jitter_eth * 1e6, 'r', 'LineWidth', 2);
hold on;
plot(distance, jitter_fiber * 1e6, 'b', 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Jitter (\mus)');
title('Jitter Analysis');
legend('10G Ethernet', 'Optical Fiber');

%% Conclusion
fprintf('Trade Study Results:\n');
fprintf(' - Optical Fiber shows significant advantages in throughput, latency, BER, and power efficiency.\n');
fprintf(' - 10G Ethernet suffers from increased jitter, higher latency, and radiation susceptibility.\n');
fprintf(' - Optical Fiber is the preferred option for space-based high-speed payloads.\n');
