% Five-antenna single-carrier 16-QAM TR comparison on measured ARA CSI.
% No OFDM is used. The measured frequency-bin CSI is used only to estimate
% h(t); the communication waveform is single-carrier pulse-shaped 16-QAM.

clear; clc; close all;
rng(7, 'twister');

root_dir = fileparts(mfilename('fullpath'));
csi_root = fullfile(root_dir, 'sklk_csi_data_dec_24_to_may_25', 'csi_data');
if ~exist(csi_root, 'dir')
    csi_root = fullfile('C:\Users\mdsah\Desktop\Skylark Data', ...
        'sklk_csi_data_dec_24_to_may_25', 'csi_data');
end
out_dir = fullfile(root_dir, 'tr_5antenna_results');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

bandwidth_hz = 24e6;
sample_rate = bandwidth_hz;
samples_per_symbol = 128;
symbol_rate = sample_rate / samples_per_symbol;
rrc_rolloff = 0.25;
rrc_span_symbols = 12;

num_bits = 80000;
bits_per_symbol = 4;
snr_db = 0:2:18;
display_snr_db = 18;
equalizer_length = 48;
max_equalizer_design_snr_db = 20;
num_tx_antennas = 5;

assert(mod(num_bits, bits_per_symbol) == 0);
rrc = rootRaisedCosine(rrc_rolloff, rrc_span_symbols, samples_per_symbol);

%% Load one measured ARA impulse response
files = dir(fullfile(csi_root, '**', '*.h5'));
if isempty(files), error('No CSI .h5 files found under %s', csi_root); end

h_channels = [];
selected_file = "";
selected_meta = struct();
for n = 1:numel(files)
    csi_file = fullfile(files(n).folder, files(n).name);
    try
        [h_candidate, meta] = loadAraMultiAntennaChannels( ...
            csi_file, bandwidth_hz, num_tx_antennas);
        if size(h_candidate,1) >= 32 && size(h_candidate,2) == num_tx_antennas
            h_channels = h_candidate;
            selected_file = string(csi_file);
            selected_meta = meta;
            break;
        end
    catch
    end
end
if isempty(h_channels), error('No readable five-antenna CSI set found.'); end

% One common calibration preserves relative antenna gains and phases.
h_channels = h_channels / sqrt(mean(sum(abs(h_channels).^2,1)));

% Fixed-total-power equal-weight five-antenna transmission without TR.
h_no_tr = sum(h_channels,2) / sqrt(num_tx_antennas);

% Fixed-total-power five-antenna TR precoding and coherent receive sum.
h_tr_filters = conj(flipud(h_channels));
h_tr_filters = h_tr_filters / sqrt(sum(abs(h_tr_filters(:)).^2));
h_effective_tr = zeros(2*size(h_channels,1)-1,1);
for ant = 1:num_tx_antennas
    h_effective_tr = h_effective_tr + ...
        conv(h_tr_filters(:,ant),h_channels(:,ant));
end

fprintf('Selected CSI file:\n%s\n', selected_file);
fprintf('UE handle %g, UE channel %g, %g BS antennas, CSI bins %g\n', ...
    selected_meta.handle,selected_meta.ue_channel, ...
    selected_meta.num_antennas,selected_meta.num_bins);
fprintf('Selected BS antenna indices: '); fprintf('%g ',selected_meta.bs_ant_indices); fprintf('\n');
fprintf('Sample rate %.2f MHz, symbol rate %.2f Msym/s, no OFDM.\n', ...
    sample_rate/1e6, symbol_rate/1e6);

%% Single-carrier 16-QAM transmitter
tx_bits = randi([0 1], num_bits, 1);
tx_symbols = qam16Mod(tx_bits);
tx_up = zeros((numel(tx_symbols)-1)*samples_per_symbol + 1, 1);
tx_up(1:samples_per_symbol:end) = tx_symbols;
tx_waveform = conv(tx_up, rrc);

%% Simulate no-TR and TR cases
[ber_no_tr, evm_no_tr_db, rx_plot_no_tr] = simulateCase( ...
    tx_waveform, tx_symbols, tx_bits, rrc, h_no_tr, samples_per_symbol, ...
    equalizer_length, max_equalizer_design_snr_db, ...
    snr_db, display_snr_db);

[ber_tr, evm_tr_db, rx_plot_tr] = simulateCase( ...
    tx_waveform, tx_symbols, tx_bits, rrc, h_effective_tr, samples_per_symbol, ...
    equalizer_length, max_equalizer_design_snr_db, ...
    snr_db, display_snr_db);

%% Channel and focusing metrics
delay_no_tr_ns = (0:numel(h_no_tr)-1).' / sample_rate * 1e9;
delay_tr_ns = (0:numel(h_effective_tr)-1).' / sample_rate * 1e9;
peak_to_avg_no_tr_db = 20*log10(max(abs(h_no_tr)) / mean(abs(h_no_tr)+eps));
peak_to_avg_tr_db = 20*log10(max(abs(h_effective_tr)) / mean(abs(h_effective_tr)+eps));

fprintf('\nChannel peak-to-average magnitude:\n');
fprintf('  Without TR: %.2f dB\n', peak_to_avg_no_tr_db);
fprintf('  With TR:    %.2f dB\n', peak_to_avg_tr_db);
fprintf('\n Es/N0(dB)      BER noTR      BER TR      EVM noTR(dB)    EVM TR(dB)\n');
for k = 1:numel(snr_db)
    fprintf(' %8.1f     %9.4g   %9.4g      %9.2f     %9.2f\n', ...
        snr_db(k), ber_no_tr(k), ber_tr(k), evm_no_tr_db(k), evm_tr_db(k));
end

%% Plots
fig = figure('Color','w','Position',[0 0 1000 500]);
stem(delay_no_tr_ns, abs(h_no_tr)/max(abs(h_no_tr)), 'filled', 'MarkerSize', 3);
grid on; xlabel('Relative delay (ns)'); ylabel('Normalized |h(t)|');
title('Five antennas, equal-weight transmission, no TR');
formatFigure(fig);
exportgraphics(fig,fullfile(out_dir,'01_channel_without_tr.png'),'Resolution',300);
savefig(fig,fullfile(out_dir,'01_channel_without_tr.fig'));

fig = figure('Color','w','Position',[0 0 1000 500]);
plot(delay_tr_ns,abs(h_effective_tr)/max(abs(h_effective_tr)),'LineWidth',2.5);
grid on; xlabel('Relative delay (ns)'); ylabel('Normalized magnitude');
title('Five-antenna TR-focused effective channel');
formatFigure(fig);
exportgraphics(fig,fullfile(out_dir,'02_channel_with_tr.png'),'Resolution',300);
savefig(fig,fullfile(out_dir,'02_channel_with_tr.fig'));

fig = figure('Color','w','Position',[0 0 1000 500]);
semilogy(snr_db,max(ber_no_tr,0.5/num_bits),'-o','LineWidth',2.5, ...
    'MarkerSize',9); hold on;
semilogy(snr_db,max(ber_tr,0.5/num_bits),'-s','LineWidth',2.5, ...
    'MarkerSize',9);
grid on; xlabel('E_s/N_0 (dB)'); ylabel('BER');
title('Single-carrier 16-QAM BER');
xlim([0 18]); ylim([0.5/num_bits 1]);
l = legend('Without TR','With TR','Location','southwest','Interpreter','none');
l.FontSize = 12;
formatFigure(fig);
set(gca,'XMinorGrid','off','YMinorGrid','off', ...
    'GridAlpha',0.24,'GridColor',[0.68 0.68 0.68]);
set(gca,'XMinorGrid','on','YMinorGrid','on', ...
    'MinorGridLineStyle',':','MinorGridAlpha',0.12, ...
    'MinorGridColor',[0.76 0.76 0.76],'LineWidth',2);
exportgraphics(fig,fullfile(out_dir,'03_ber_comparison.png'),'Resolution',300);
savefig(fig,fullfile(out_dir,'03_ber_comparison.fig'));

fig = figure('Color','w','Position',[0 0 1000 500]);
evm_no_tr_percent = 100 * 10.^(evm_no_tr_db/20);
evm_tr_percent = 100 * 10.^(evm_tr_db/20);
plot(snr_db,evm_no_tr_percent,'-o','LineWidth',2.5,'MarkerSize',9); hold on;
plot(snr_db,evm_tr_percent,'-s','LineWidth',2.5,'MarkerSize',9);
grid on; xlabel('E_s/N_0 (dB)'); ylabel('RMS EVM (%)');
title('Single-carrier 16-QAM EVM');
xlim([0 18]);
l = legend('Without TR','With TR','Location','southwest','Interpreter','none');
l.FontSize = 12;
formatFigure(fig);
exportgraphics(fig,fullfile(out_dir,'04_evm_comparison.png'),'Resolution',300);
savefig(fig,fullfile(out_dir,'04_evm_comparison.fig'));

levels = [-3 -1 1 3] / sqrt(10);
[ii, qq] = meshgrid(levels, levels);
ideal = ii(:) + 1j*qq(:);

fig = figure('Color','w','Position',[0 0 1000 500]);
scatter(real(rx_plot_no_tr), imag(rx_plot_no_tr), 9, 'filled'); hold on;
plot(real(ideal), imag(ideal), 'kx', 'LineWidth', 2, 'MarkerSize', 11);
xlim([-1.5 1.5]); ylim([-1.5 1.5]); pbaspect([1 1 1]);
grid on; xlabel('In-phase'); ylabel('Quadrature');
title(sprintf('Without TR, %g dB', display_snr_db));
formatFigure(fig,[0 0 700 700]);
exportgraphics(fig,fullfile(out_dir,'05_constellation_without_tr_18dB.png'), ...
    'Resolution',300);
savefig(fig,fullfile(out_dir,'05_constellation_without_tr_18dB.fig'));

fig = figure('Color','w','Position',[0 0 1000 500]);
scatter(real(rx_plot_tr), imag(rx_plot_tr), 9, 'filled'); hold on;
plot(real(ideal), imag(ideal), 'kx', 'LineWidth', 2, 'MarkerSize', 11);
xlim([-1.5 1.5]); ylim([-1.5 1.5]); pbaspect([1 1 1]);
grid on; xlabel('In-phase'); ylabel('Quadrature');
title(sprintf('With TR, %g dB', display_snr_db));
formatFigure(fig,[0 0 700 700]);
exportgraphics(fig,fullfile(out_dir,'06_constellation_with_tr_18dB.png'), ...
    'Resolution',300);
savefig(fig,fullfile(out_dir,'06_constellation_with_tr_18dB.fig'));

results = table(snr_db(:), ber_no_tr(:), ber_tr(:), evm_no_tr_db(:), evm_tr_db(:), ...
    'VariableNames', {'EsN0_dB','BER_without_TR','BER_with_TR','EVM_without_TR_dB','EVM_with_TR_dB'});
writetable(results, fullfile(out_dir, 'ara_5antenna_tr_ber_evm.csv'));
save(fullfile(out_dir, 'ara_5antenna_tr_ber_evm.mat'), ...
    'results','h_no_tr','h_tr_filters','h_effective_tr','h_channels', ...
    'selected_file','selected_meta', ...
    'sample_rate','symbol_rate','samples_per_symbol','rrc_rolloff');

fprintf('\nSaved outputs in:\n%s\n', out_dir);

%% Local functions
function [h_channels, meta] = loadAraMultiAntennaChannels(csi_file, bandwidth_hz, num_antennas)
    h5info(csi_file);
    bin = double(h5read(csi_file, '/ue_uplink/bin'));
    ch = double(h5read(csi_file, '/ue_uplink/ch'));
    hdl = double(h5read(csi_file, '/ue_uplink/hdl'));
    raw = double(h5read(csi_file, '/ue_uplink/data'));

    nrec = numel(bin);
    vals_per_rec = numel(raw) / nrec;
    if mod(vals_per_rec, 2) ~= 0, error('CSI data is not real/imag paired.'); end

    raw = reshape(raw, vals_per_rec, nrec);
    csi = raw(1:2:end,:) + 1j*raw(2:2:end,:);
    csi(abs(csi) > 1e6) = NaN;
    csi(csi == 0) = NaN;

    handles = unique(hdl(:))';
    ue_chans = unique(ch(:))';
    best = struct('score',-inf,'handle',NaN,'ue_channel',NaN);
    for hh = handles
        for cc = ue_chans
            rec_idx = (hdl == hh) & (ch == cc);
            if nnz(rec_idx) < 8, continue; end
            % Antennas are sparse across repeated raw records. Validate
            % antenna coverage only after records are averaged by CSI bin.
            score = nnz(rec_idx);
            if score > best.score
                best.score = score;
                best.handle = hh;
                best.ue_channel = cc;
            end
        end
    end
    if ~isfinite(best.score), error('No valid CSI records.'); end

    use = (hdl == best.handle) & (ch == best.ue_channel);
    bins = unique(bin(use));
    Hf_all = NaN(numel(bins),size(csi,1));
    for k = 1:numel(bins)
        idx = use & (bin == bins(k));
        Hf_all(k,:) = mean(csi(:,idx),2,'omitnan').';
    end

    coverage = sum(isfinite(Hf_all),1);
    antenna_energy = mean(abs(Hf_all).^2,1,'omitnan');
    eligible = find(coverage >= 0.8*numel(bins) & isfinite(antenna_energy));
    if numel(eligible) < num_antennas
        error('Too few BS antennas have sufficient bin coverage.');
    end
    [~, order] = sort(bins);
    bins = bins(order);
    bins = bins(:);
    [~,energy_order] = sort(antenna_energy(eligible),'descend');
    selected_antennas = eligible(energy_order(1:num_antennas));
    Hf = Hf_all(order,selected_antennas);
    for ant = 1:num_antennas
        valid = isfinite(Hf(:,ant));
        if nnz(valid) < 32, error('Too few valid bins on selected antenna.'); end
        Hf(:,ant) = interp1(bins(valid),real(Hf(valid,ant)),bins, ...
            'linear','extrap') + 1j*interp1(bins(valid), ...
            imag(Hf(valid,ant)),bins,'linear','extrap');
    end
    h_channels = ifft(ifftshift(Hf,1),[],1);

    meta = best;
    meta.num_bins = size(Hf,1);
    meta.num_antennas = num_antennas;
    meta.bs_ant_indices = selected_antennas;
    meta.duration_ns = size(h_channels,1) / bandwidth_hz * 1e9;
end

function [ber, evm_db, rx_symbols_for_plot] = simulateCase(tx_waveform, tx_symbols, tx_bits, rrc, channel, sps, eq_len, ~, snr_db, display_snr_db)
    clean = conv(tx_waveform, channel);
    sample_channel = conv(conv(rrc, channel), rrc);
    [~, peak] = max(abs(sample_channel));
    phase = mod(peak-1, sps) + 1;
    symbol_channel = sample_channel(phase:sps:end);

    ber = zeros(size(snr_db));
    evm_db = zeros(size(snr_db));
    rx_symbols_for_plot = [];
    [raw_peak_power,raw_peak_index] = max(abs(symbol_channel).^2);
    raw_sidelobe_energy = sum(abs(symbol_channel).^2)-raw_peak_power;
    raw_sir_db = 10*log10(raw_peak_power/max(raw_sidelobe_energy,eps));

    for k = 1:numel(snr_db)
        % Fixed transmit-referenced Es/N0 preserves five-antenna TR gain.
        % Both cases use identical total transmit energy and AWGN variance.
        symbol_energy = mean(abs(tx_symbols).^2);
        noise_var = symbol_energy * 10^(-snr_db(k)/10);
        noise = sqrt(noise_var/2) * (randn(size(clean)) + 1j*randn(size(clean)));
        rx_matched = conv(clean + noise, rrc);
        rx_sym_rate = rx_matched(phase:sps:end);

        if raw_sir_db >= 20
            % A clean TR focus needs only complex gain correction. A long
            % inverse filter would deepen channel nulls and amplify noise.
            w = 1;
            delay = raw_peak_index;
            gain = symbol_channel(raw_peak_index);
        else
            [w,delay,gain] = designAdaptiveMmse( ...
                symbol_channel,eq_len,noise_var);
        end
        equalized = conv(rx_sym_rate, w);

        idx = delay + (0:numel(tx_symbols)-1).';
        if idx(end) > numel(equalized)
            idx = idx(idx <= numel(equalized));
        end
        rx = equalized(idx) / gain;
        ref = tx_symbols(1:numel(rx));
        ref_bits = tx_bits(1:4*numel(rx));

        rx_bits = qam16Demod(rx);
        ber(k) = mean(rx_bits ~= ref_bits);
        evm_rms = sqrt(mean(abs(rx - ref).^2) / mean(abs(ref).^2));
        evm_db(k) = 20*log10(evm_rms + eps);

        if snr_db(k) == display_snr_db
            rx_symbols_for_plot = rx(1:min(6000,numel(rx)));
        end
    end
end

function symbols = qam16Mod(bits)
    groups = reshape(bits, 4, []).';
    i = (1 - 2*groups(:,1)) .* (3 - 2*groups(:,2));
    q = (1 - 2*groups(:,3)) .* (3 - 2*groups(:,4));
    symbols = (i + 1j*q) / sqrt(10);
end

function bits = qam16Demod(symbols)
    s = symbols * sqrt(10);
    groups = [real(s) < 0, abs(real(s)) < 2, imag(s) < 0, abs(imag(s)) < 2];
    bits = reshape(groups.', [], 1);
end

function h = rootRaisedCosine(alpha, span_symbols, sps)
    t = (-span_symbols*sps/2:span_symbols*sps/2).' / sps;
    h = zeros(size(t));
    tol = 100*eps;
    for k = 1:numel(t)
        if abs(t(k)) < tol
            h(k) = 1 - alpha + 4*alpha/pi;
        elseif alpha > 0 && abs(abs(t(k)) - 1/(4*alpha)) < tol
            h(k) = alpha/sqrt(2) * ((1+2/pi)*sin(pi/(4*alpha)) + ...
                (1-2/pi)*cos(pi/(4*alpha)));
        else
            h(k) = (sin(pi*t(k)*(1-alpha)) + 4*alpha*t(k)*cos(pi*t(k)*(1+alpha))) / ...
                (pi*t(k)*(1 - (4*alpha*t(k))^2));
        end
    end
    h = h / norm(h);
end

function A = makeConvolutionMatrix(channel, filter_length)
    channel = channel(:);
    A = zeros(numel(channel) + filter_length - 1, filter_length);
    for col = 1:filter_length
        A(col:col+numel(channel)-1, col) = channel;
    end
end

function [best_w,best_delay,best_gain] = designAdaptiveMmse(channel,max_length,noise_var)
    candidate_lengths = unique(min([1 2 4 8 12 16 24 32 48],max_length));
    [~,channel_peak] = max(abs(channel));
    best_mse = inf;
    best_w = 1;
    best_delay = channel_peak;
    best_gain = channel(channel_peak);

    for filter_length = candidate_lengths
        A = makeConvolutionMatrix(channel,filter_length);
        mu = noise_var;
        R = A'*A + mu*eye(filter_length);
        first_delay = channel_peak;
        last_delay = min(size(A,1),channel_peak+filter_length-1);
        for delay = first_delay:last_delay
            target = zeros(size(A,1),1);
            target(delay) = 1;
            w = R \ (A'*target);
            effective = A*w;
            gain = effective(delay);
            if abs(gain) < 1e-8, continue; end
            residual = effective;
            residual(delay) = 0;
            output_mse = (sum(abs(residual).^2) + ...
                noise_var*sum(abs(w).^2)) / abs(gain)^2;
            if output_mse < best_mse
                best_mse = output_mse;
                best_w = w;
                best_delay = delay;
                best_gain = gain;
            end
        end
    end
end

function formatFigure(fig,position)
    if nargin < 2, position = [0 0 1000 500]; end
    set(fig,'Position',position,'Color','White');
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        set(ax(k),'FontName','Times New Roman','FontSize',24, ...
            'LineWidth',3,'GridLineStyle','--','Box','on');
        grid(ax(k),'on');
    end
end
