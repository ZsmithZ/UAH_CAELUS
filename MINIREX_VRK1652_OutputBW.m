
ibw_min = 300e6;
ibw_max = 4e9;
n_points = 100;
n_bits = 14;
oversample_factor = 2.5;
n_tiles = 4;
enable_IQ = true;

capture_time_rt = 1e-3;
dma_latency = 20e-6;
safety_margin = 1.2;

burst_time = 25e-6;
prf_duty_cycle = 0.1;

waveforms = struct( ...
    'BPSK',  struct('decimation', 4, 'bps', 1, 'symb_rate_factor', 0.4, 'burst_mode', false, 'designated_ibw', 10e6), ...
    'QPSK',  struct('decimation', 4, 'bps', 2, 'symb_rate_factor', 0.4, 'burst_mode', false, 'designated_ibw', 25e6), ...
    'OQPSK', struct('decimation', 4, 'bps', 2, 'symb_rate_factor', 0.4, 'burst_mode', false, 'designated_ibw', 20e6), ...
    'PSK8',  struct('decimation', 6, 'bps', 3, 'symb_rate_factor', 0.35, 'burst_mode', false, 'designated_ibw', 35e6), ...
    'QAM16', struct('decimation', 8, 'bps', 4, 'symb_rate_factor', 0.3, 'burst_mode', false, 'designated_ibw', 45e6), ...
    'LFM',   struct('decimation', 2, 'bps', NaN, 'symb_rate_factor', NaN, 'burst_mode', true, 'designated_ibw', 1e9) ...
);

% Calculation reference: Rs = sym_rate_factor * designated_ibw = XXX Msps
% Bit rate = XXX Msps * bit per symbol

ibw_range = linspace(ibw_min, ibw_max, n_points);
waveform_names = fieldnames(waveforms);
iq_factor = 2 ^ enable_IQ;

figure;
for i = 1:length(waveform_names)
    name = waveform_names{i};
    wf = waveforms.(name);
    is_burst = wf.burst_mode;
    fs = oversample_factor .* ibw_range;
    raw_bw = iq_factor .* fs .* n_bits .* n_tiles;

    if ~is_burst
        symb_rate = wf.symb_rate_factor .* ibw_range;
        eff_bw = symb_rate .* wf.bps .* n_tiles;
        mem_req = eff_bw .* (capture_time_rt + dma_latency) / 8 / 1e6;
        mem_req = mem_req * safety_margin;
        label_mode = 'Streaming (Symbol Rate)';
    else
        eff_bw = raw_bw ./ wf.decimation;
        mem_req = eff_bw .* burst_time / 8 / 1e6;
        mem_req = mem_req / prf_duty_cycle * safety_margin;
        label_mode = 'Burst Mode (Full BW)';
    end

    [~, idx_annot] = min(abs(ibw_range - wf.designated_ibw));
    annotated_bw = eff_bw(idx_annot) / 1e9;
    annotated_mem = mem_req(idx_annot);

    subplot(3, 2, i);
    yyaxis left;
    plot(ibw_range / 1e9, eff_bw / 1e9, 'b-', 'LineWidth', 2);
    ylabel('Output BW (Gbps)');
    ax = gca;
    ax.YColor = 'b';

    yyaxis right;
    plot(ibw_range / 1e9, mem_req, 'Color', [0 0.6 0], 'LineStyle', '--', 'LineWidth', 2);
    ylabel('Memory Req (MB)');
    ax = gca;
    ax.YColor = [0 0.6 0];
    ylim([0 max(mem_req) * 1.2]);

    title(sprintf('%s (%s)', name, label_mode), 'FontWeight', 'bold');
    xlabel('Instantaneous Bandwidth (GHz)');
    legend('Output BW (Gbps)', 'Memory Req (MB)', 'Location', 'northwest');
    grid on;
    xlim([ibw_min ibw_max] / 1e9);

    % Annotated vertical line at designated IBW
    xline(wf.designated_ibw / 1e9, 'r-.', ...
        sprintf('Designated IBW: %.1f MHz', wf.designated_ibw / 1e6), ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom', ...
        'Color', 'r', 'LineWidth', 1.5);

    % Single-line annotation textbox
    x_annot = wf.designated_ibw / 1e9 + 0.4;
    y_annot = 0.6 * max(mem_req) * 1.2;
    annotation_str = sprintf('%.1f Gbps / %.1f MB @ %.1f MHz IBW', ...
                             annotated_bw, annotated_mem, wf.designated_ibw/1e6);
    text(x_annot, y_annot, annotation_str, ...
         'EdgeColor', 'black', 'BackgroundColor', 'white', 'FontSize', 8);
end

sgtitle('VRK1652 Output Bandwidt / Memory Requirements');
