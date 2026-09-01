%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  A02_interperformer_correlation
%
%  Inter-performer respiratory synchrony: lag-0 correlation + surrogate test
%
%  Project : Nagauta respiratory coordination
%
%  Description
%  -----------
%  Computes the lag-0 Pearson correlation (Fisher z-transformed) between
%  every pair of performers (dyad), separately for each musical section
%  and take. Dyad values are averaged across sections and takes, yielding
%  one Fisher z value per dyad (n = 15 dyads from 6 performers).
%
%  Two inferential approaches are reported:
%
%  1. Descriptive one-sample t-test
%     Tests dyad-averaged Fisher z against 0.
%     NOTE: the 15 dyads are not independent (each performer appears in
%     multiple dyads), so this test is reported for descriptive purposes
%     only. The circular-shift surrogate test below is the primary test.
%
%  2. Circular-shift surrogate test
%     For each iteration, every performer's signal is independently
%     circularly shifted by a random amount (>= min_shift_s). The same
%     shift is applied across all dyads for that performer, preserving
%     within-performer autocorrelation structure while destroying
%     between-performer simultaneity.
%     Test statistic: grand mean of all dyad Fisher z values, averaged
%     over dyads, sections, and takes.
%
%  Outputs: three figures
%    Figure 1 – section × dyad heatmap of Fisher z (mean across takes)
%    Figure 2 – dyad-level values per take + dyad-averaged summary
%    Figure 3 – surrogate null distribution with observed statistic
%
%  Folder structure assumed
%  ------------------------
%  <project_root>/
%    Mfile/        <- this file lives here
%    Mat_file/
%      Hexoskin_data.mat    (respiration, ECG, acceleration per performer and take)
%      SectionTimepoint.mat (section boundary times in ms)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
rng(20260205);   % fixed seed for reproducibility

%% ---- Path setup --------------------------------------------------------
code_dir    = fileparts(mfilename('fullpath'));
project_dir = fileparts(code_dir);
data_dir    = fullfile(project_dir, 'Mat_file');

%% ---- Parameters --------------------------------------------------------
fs_raw            = 128;    % original sampling rate (Hz)
fs_ds             = 32;     % target sampling rate after downsampling (Hz)
lowpass_cutoff_hz = 5;      % low-pass filter cutoff (Hz)

n_surrogate    = 1000;      % number of circular-shift surrogate iterations
min_shift_s    = 15;        % minimum circular shift (seconds)

downsample_factor = fs_raw / fs_ds;
ms_per_sample     = 1000   / fs_ds;
min_shift_samp    = round(min_shift_s * fs_ds);

%% ---- Experiment metadata -----------------------------------------------
take_field    = {'S1', 'S2'};
performers    = {'Uta','Shamisen1','Kotsuzumi','Fue','Taiko','Shamisen2'};
section_label = {'A','B','C','D','E','F','G','H'};

n_take      = numel(take_field);
n_performer = numel(performers);
n_section   = numel(section_label);

dyad_idx = nchoosek(1:n_performer, 2);   % [pf_a, pf_b] per dyad (15 rows)
n_dyad   = size(dyad_idx, 1);

dyad_label = strings(n_dyad, 1);
for dy = 1:n_dyad
    dyad_label(dy) = performers{dyad_idx(dy,1)} + "–" + performers{dyad_idx(dy,2)};
end

%% ---- Load data ---------------------------------------------------------
load(fullfile(data_dir, 'Hexoskin_data.mat'));    % struct 'Nagauta'
load(fullfile(data_dir, 'SectionTimepoint.mat')); % matrix 'Timepoint' (ms)

%% ---- Pre-extract z-scored section segments -----------------------------
%  segment{take, section}    : [L x n_performer] matrix, z-scored within section
%  valid_perf{take, section} : logical [n_performer x 1]
%
%  Pre-extracting segments ensures the observed and surrogate statistics
%  are computed on identical data.

segment    = cell(n_take, n_section);
valid_perf = cell(n_take, n_section);

for tk = 1:n_take

    resp_proc = cell(n_performer, 1);
    for pf = 1:n_performer
        raw = Nagauta.(take_field{tk}).(performers{pf}).respiration_thoracic;
        raw = raw(:);
        sig = detrend(raw);
        sig = lowpass(sig, lowpass_cutoff_hz, fs_raw);
        resp_proc{pf} = downsample(sig, downsample_factor);
    end

    for se = 1:n_section
        seg_start = max(1, round(Timepoint(se,   tk) / ms_per_sample));
        seg_end   =        round(Timepoint(se+1, tk) / ms_per_sample);

        seg_len = zeros(n_performer, 1);
        for pf = 1:n_performer
            e = min(numel(resp_proc{pf}), seg_end);
            seg_len(pf) = max(0, e - seg_start + 1);
        end

        has_data = seg_len > 0;
        if nnz(has_data) < 2
            segment{tk,se}    = [];
            valid_perf{tk,se} = false(n_performer, 1);
            continue;
        end

        L = min(seg_len(has_data));

        seg_mat = nan(L, n_performer);
        vperf   = false(n_performer, 1);

        for pf = 1:n_performer
            if seg_len(pf) < L || L < 10, continue; end
            col = resp_proc{pf}(seg_start : seg_start + L - 1);
            if any(isnan(col)) || std(col) == 0, continue; end
            seg_mat(:, pf) = (col - mean(col)) / std(col);
            vperf(pf)      = true;
        end

        segment{tk,se}    = seg_mat;
        valid_perf{tk,se} = vperf;
    end
end

%% ---- Observed statistic ------------------------------------------------
[obs_grand_mean_z, fisherz_dyad_section_take] = ...
    compute_dyad_grandmean(segment, valid_perf, dyad_idx, min_shift_samp, false);

fisherz_by_dyad = squeeze(mean(fisherz_dyad_section_take, [2 3], 'omitnan'));  % n_dyad x 1

%% ---- Circular-shift surrogate test -------------------------------------
surr_grand_mean = nan(n_surrogate, 1);
for b = 1:n_surrogate
    surr_grand_mean(b) = ...
        compute_dyad_grandmean(segment, valid_perf, dyad_idx, min_shift_samp, true);
end

null_mean = mean(surr_grand_mean);
null_lo   = prctile(surr_grand_mean,  2.5);
null_hi   = prctile(surr_grand_mean, 97.5);
pct_obs   = 100 * mean(surr_grand_mean < obs_grand_mean_z);

p_one_sided = (1 + sum(surr_grand_mean >= obs_grand_mean_z)) / (n_surrogate + 1);
p_two_sided = (1 + sum(abs(surr_grand_mean - null_mean) >= ...
                        abs(obs_grand_mean_z - null_mean))) / (n_surrogate + 1);

%% ---- Descriptive t-test ------------------------------------------------
[~, p_ttest, ci_ttest, stats_ttest] = ttest(fisherz_by_dyad, 0);
mean_z = mean(fisherz_by_dyad, 'omitnan');

%% ---- Print summary -----------------------------------------------------
se_z = std(fisherz_by_dyad) / sqrt(n_dyad);

fprintf('=== Inter-performer respiratory synchrony ===\n\n');
fprintf('Sample unit: dyad-averaged Fisher z, n = %d dyads (from %d performers)\n', ...
        n_dyad, n_performer);
fprintf('Averaged over %d sections x %d takes\n\n', n_section, n_take);

fprintf('Descriptive statistics (dyad-level Fisher z):\n');
fprintf('  Mean = %.4f, SD = %.4f, SE = %.4f\n', mean_z, std(fisherz_by_dyad), se_z);
fprintf('  Approx. mean r = tanh(mean z) = %.4f\n', tanh(mean_z));
fprintf('  Range: [%.4f, %.4f]\n\n', min(fisherz_by_dyad), max(fisherz_by_dyad));

fprintf('--- Descriptive one-sample t-test (dyad-averaged z vs. 0) ---\n');
fprintf('NOTE: dyads are not independent; treat as descriptive only.\n');
fprintf('t(%d) = %.3f, p = %.4f\n', stats_ttest.df, stats_ttest.tstat, p_ttest);
fprintf('95%% CI: [%.4f, %.4f]\n\n', ci_ttest(1), ci_ttest(2));

fprintf('--- Circular-shift surrogate test (primary test) ---\n');
fprintf('Surrogates = %d, min shift = %.0f s (%d samples)\n', ...
        n_surrogate, min_shift_s, min_shift_samp);
fprintf('Observed grand-mean z = %.4f (approx. r = %.4f)\n', ...
        obs_grand_mean_z, tanh(obs_grand_mean_z));
fprintf('Surrogate null: mean = %.4f, 95%% interval [%.4f, %.4f]\n', ...
        null_mean, null_lo, null_hi);
fprintf('Observed percentile within null = %.1f%%\n', pct_obs);
fprintf('p (one-sided, observed > chance) = %.4f\n', p_one_sided);
fprintf('p (two-sided)                    = %.4f\n\n', p_two_sided);

fprintf('--- Per-dyad Fisher z (averaged across sections and takes) ---\n');
fprintf('%-30s  Fisher z   Approx. r\n', 'Dyad');
for dy = 1:n_dyad
    fprintf('%-30s  %+.4f     %+.4f\n', dyad_label(dy), ...
            fisherz_by_dyad(dy), tanh(fisherz_by_dyad(dy)));
end

%% ---- Figure 1: section × dyad heatmap (mean across takes) --------------
fisherz_dyad_section = squeeze(mean(fisherz_dyad_section_take, 3, 'omitnan'));

figure(1);
set(gcf, 'Units','centimeters', 'Position',[2 2 17 10], ...
         'PaperUnits','centimeters', 'PaperPosition',[0 0 17 10]);
imagesc(fisherz_dyad_section);
clim([-0.4 0.4]);
colorbar;
set(gca, 'XTick',1:n_section, 'XTickLabel',section_label, ...
         'YTick',1:n_dyad,    'YTickLabel',dyad_label, ...
         'FontName','Arial',  'FontSize',10);
xlabel('Section', 'FontName','Arial', 'FontSize',12);
ylabel('Dyad',    'FontName','Arial', 'FontSize',12);
title('Lag-0 Fisher z by dyad and section (mean across takes)', ...
      'FontName','Arial', 'FontSize',11);
axis square;

%% ---- Figure 2: dyad-level values, take-wise + average ------------------
fisherz_dyad_take = squeeze(mean(fisherz_dyad_section_take, 2, 'omitnan'));

figure(2);
set(gcf, 'Units','centimeters', 'Position',[2 2 17 10], ...
         'PaperUnits','centimeters', 'PaperPosition',[0 0 17 10]);
hold on;
for dy = 1:n_dyad
    plot([1 2], fisherz_dyad_take(dy,:), '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1);
end
scatter(1.5*ones(n_dyad,1), fisherz_by_dyad, 40, 'k', 'filled');
yline(0, 'k--', 'LineWidth',1);
xlim([0.7 2.3]);
set(gca, 'XTick',[1 1.5 2], 'XTickLabel',{'Take 1','Mean','Take 2'}, ...
         'FontName','Arial', 'FontSize',10);
ylabel('Fisher z at lag 0', 'FontName','Arial', 'FontSize',12);
title('Dyad-level synchrony by take', 'FontName','Arial', 'FontSize',11);
grid on; box on; axis square;

%% ---- Figure 3: surrogate null with observed value ----------------------
figure(3);
set(gcf, 'Units','centimeters', 'Position',[2 2 17 10], ...
         'PaperUnits','centimeters', 'PaperPosition',[0 0 17 10]);
histogram(surr_grand_mean, 40, 'FaceColor',[0.75 0.75 0.75], 'EdgeColor','none');
hold on;
yl = ylim;
plot([obs_grand_mean_z obs_grand_mean_z], yl, 'r-',  'LineWidth',2, 'DisplayName','Observed');
plot([null_lo null_lo],                   yl, 'k--', 'LineWidth',1, 'DisplayName','95% null interval');
plot([null_hi null_hi],                   yl, 'k--', 'LineWidth',1, 'HandleVisibility','off');
xlabel('Grand-mean dyad Fisher z at lag 0', 'FontName','Arial', 'FontSize',12);
ylabel('Surrogate count',                   'FontName','Arial', 'FontSize',12);
title(sprintf('Circular-shift surrogate null  (p_{one-sided} = %.4f)', p_one_sided), ...
      'FontName','Arial', 'FontSize',11);
legend({'Surrogate null','Observed','95% null interval'}, 'Location','best');
set(gca, 'FontName','Arial', 'FontSize',10);
box on;

%% ======================================================================
%  Local functions
%  ======================================================================
function [grand_mean_z, fisherz_acc] = ...
        compute_dyad_grandmean(segment, valid_perf, dyad_idx, min_shift_samp, do_shift)
% Grand mean of lag-0 Fisher z across all dyads, sections, and takes.
%
% When do_shift is true, each valid performer's signal is independently
% circularly shifted by a random amount drawn uniformly from
% [min_shift_samp, L - min_shift_samp]. One shift per performer per
% take×section block, reused across all dyads involving that performer.

    [n_take, n_section] = size(segment);
    n_dyad = size(dyad_idx, 1);

    fisherz_acc = nan(n_dyad, n_section, n_take);

    for tk = 1:n_take
        for se = 1:n_section
            M = segment{tk, se};
            v = valid_perf{tk, se};
            if isempty(M), continue; end
            L = size(M, 1);

            if do_shift
                lo = max(1, min(min_shift_samp, floor(L / 4)));
                hi = max(lo + 1, L - lo);
                for pf = 1:size(M, 2)
                    if v(pf)
                        M(:, pf) = circshift(M(:, pf), randi([lo, hi]));
                    end
                end
            end

            for dy = 1:n_dyad
                pf_a = dyad_idx(dy, 1);
                pf_b = dyad_idx(dy, 2);
                if v(pf_a) && v(pf_b)
                    r0 = corr(M(:, pf_a), M(:, pf_b));
                    r0 = max(min(r0, 0.999999), -0.999999);
                    fisherz_acc(dy, se, tk) = atanh(r0);
                end
            end
        end
    end

    z_per_dyad   = squeeze(mean(fisherz_acc, [2 3], 'omitnan'));
    grand_mean_z = mean(z_per_dyad, 'omitnan');
end
