%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  A01_dtw_self_vs_other
%
%  Dynamic Time Warping (DTW) analysis: self-pairing vs. other-pairing
%
%  Project : Nagauta respiratory coordination 
%
%  Description
%  -----------
%  Computes the normalised DTW distance between the thoracic respiration
%  signals of the same performer recorded in two separate takes
%  (Take 1 / S1  vs. Take 2 / S2).
%
%  Self-pairing  : same performer across takes  (e.g. Uta-S1 vs. Uta-S2)
%  Other-pairing : different performers          (e.g. Uta-S1 vs. Shamisen1-S2)
%
%  The full 6x6 pairwise DTW distance matrix (averaged across sections)
%  is computed. Self-pairing uses the diagonal; other-pairing for each
%  S1 performer is the median of the five off-diagonal values in that row.
%  Values are averaged across sections, yielding one self-pairing and
%  one other-pairing value per performer (n = 6).
%
%  Statistical inference uses an exhaustive label-permutation test
%  (6! = 720 permutations of Take 2 performer labels), which respects the
%  dependent structure of the data. A paired t-test is additionally
%  reported for descriptive purposes.
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

%% ---- Path setup --------------------------------------------------------
code_dir    = fileparts(mfilename('fullpath'));
project_dir = fileparts(code_dir);
data_dir    = fullfile(project_dir, 'Mat_file');

%% ---- Parameters --------------------------------------------------------
fs_raw            = 128;    % original sampling rate (Hz)
fs_ds             = 32;     % target sampling rate after downsampling (Hz)
lowpass_cutoff_hz = 5;      % low-pass filter cutoff (Hz)
dtw_constraint_s  = 15;     % DTW Sakoe-Chiba band width (seconds)

downsample_factor = fs_raw / fs_ds;
dtw_constraint    = dtw_constraint_s * fs_ds;   % band width in samples
ms_per_sample     = 1000 / fs_ds;

%% ---- Experiment metadata -----------------------------------------------
take_field    = {'S1', 'S2'};
performers    = {'Uta','Shamisen1','Kotsuzumi','Fue','Taiko','Shamisen2'};
section_label = {'A','B','C','D','E','F','G','H'};

n_performer = numel(performers);
n_section   = numel(section_label);

performer_colors = [
    0.00, 0.45, 0.74;   % Uta
    0.85, 0.33, 0.10;   % Shamisen1
    0.47, 0.68, 0.19;   % Kotsuzumi
    0.49, 0.18, 0.56;   % Fue
    0.30, 0.75, 0.93;   % Taiko
    0.64, 0.08, 0.18];  % Shamisen2

%% ---- Load data ---------------------------------------------------------
load(fullfile(data_dir, 'Hexoskin_data.mat'));    % struct 'Nagauta'
load(fullfile(data_dir, 'SectionTimepoint.mat')); % matrix 'Timepoint' (ms)

%% ---- Pre-process respiration (low-pass + downsample) -------------------
resp_ds = cell(2, n_performer);
for tk = 1:2
    for pf = 1:n_performer
        raw = Nagauta.(take_field{tk}).(performers{pf}).respiration_thoracic;
        raw = raw(:);
        sig = lowpass(raw, lowpass_cutoff_hz, fs_raw);
        resp_ds{tk, pf} = downsample(sig, downsample_factor);
    end
end

%% ---- Compute full 6x6 pairwise DTW distance matrix --------------------
%  dist_mat(pf_s1, pf_s2, section) = normalised DTW distance
%  Storing all pairs enables the label-permutation test below.

dist_mat = nan(n_performer, n_performer, n_section);

for pf_s1 = 1:n_performer
    x_full = resp_ds{1, pf_s1};

    for pf_s2 = 1:n_performer
        y_full = resp_ds{2, pf_s2};

        for se = 1:n_section
            idx_start_x = max(1, round(Timepoint(se,   1) / ms_per_sample));
            idx_end_x   = min(numel(x_full), round(Timepoint(se+1, 1) / ms_per_sample));
            idx_start_y = max(1, round(Timepoint(se,   2) / ms_per_sample));
            idx_end_y   = min(numel(y_full), round(Timepoint(se+1, 2) / ms_per_sample));

            x = x_full(idx_start_x:idx_end_x);
            y = y_full(idx_start_y:idx_end_y);

            if std(x) == 0 || std(y) == 0, continue; end
            x = (x - mean(x)) / std(x);
            y = (y - mean(y)) / std(y);

            [dtw_raw, warp_path_x, ~] = dtw(x, y, dtw_constraint);
            dist_mat(pf_s1, pf_s2, se) = dtw_raw / numel(warp_path_x);
        end
    end
end

% Section-averaged distance matrix (used for permutation test)
dist_mat_avg = mean(dist_mat, 3, 'omitnan');   % [n_performer x n_performer]

%% ---- Aggregate to one value per performer ------------------------------
self_per_performer  = nan(n_performer, 1);
other_per_performer = nan(n_performer, 1);

for pf = 1:n_performer
    % Self: diagonal of section-averaged matrix
    self_per_performer(pf) = dist_mat_avg(pf, pf);

    % Other: median of the five off-diagonal values in row pf
    off_diag = dist_mat_avg(pf, setdiff(1:n_performer, pf));
    other_per_performer(pf) = median(off_diag, 'omitnan');
end

%% ---- Test statistic (shared by t-test and permutation) -----------------
% Mean difference: self minus other (negative = self < other = synchrony)
obs_stat = mean(self_per_performer - other_per_performer);

%% ---- Exhaustive label-permutation test (6! = 720) ----------------------
%  Permute the Take 2 performer labels. For each permutation, recompute the
%  self vs. other difference using the same aggregation as observed.
%  This respects the dependent structure of the 6-performer dataset.

all_perms  = perms(1:n_performer);   % 720 x 6
n_perm     = size(all_perms, 1);
perm_stats = nan(n_perm, 1);

for p = 1:n_perm
    perm = all_perms(p, :);

    self_p  = nan(n_performer, 1);
    other_p = nan(n_performer, 1);

    for pf = 1:n_performer
        % Under this permutation, Take 2 performer perm(pf) is matched to S1 performer pf
        self_p(pf) = dist_mat_avg(pf, perm(pf));

        % Other: the remaining 5 values in row pf under this permutation
        other_cols = dist_mat_avg(pf, perm(setdiff(1:n_performer, pf)));
        other_p(pf) = median(other_cols, 'omitnan');
    end

    perm_stats(p) = mean(self_p - other_p);
end

% Two-sided p: proportion of permutations with |stat| >= |observed|
p_perm_two  = mean(abs(perm_stats) >= abs(obs_stat));
% One-sided p: self < other (stat < 0)
p_perm_one  = mean(perm_stats <= obs_stat);

%% ---- Paired t-test (descriptive) ---------------------------------------
diff_vals = self_per_performer - other_per_performer;
[~, p_ttest, ci_ttest, stats_ttest] = ttest(self_per_performer, other_per_performer);
cohens_d = mean(diff_vals) / std(diff_vals);

%% ---- Descriptive statistics --------------------------------------------
self_mean  = mean(self_per_performer);
self_sd    = std(self_per_performer);
self_se    = self_sd / sqrt(n_performer);

other_mean = mean(other_per_performer);
other_sd   = std(other_per_performer);
other_se   = other_sd / sqrt(n_performer);

%% ---- Print summary -----------------------------------------------------
fprintf('--- DTW: self-pairing vs. other-pairing ---\n');
fprintf('Sample unit: normalised DTW distance per performer, n = %d performers\n\n', n_performer);

fprintf('Condition       Mean     SD       SE\n');
fprintf('Self-pairing    %.4f   %.4f   %.4f\n', self_mean,  self_sd,  self_se);
fprintf('Other-pairing   %.4f   %.4f   %.4f\n', other_mean, other_sd, other_se);

fprintf('\nObserved mean difference (self − other): %.4f\n', obs_stat);
fprintf('Complete separation: %s (max Self = %.4f, min Other = %.4f)\n', ...
        string(max(self_per_performer) < min(other_per_performer)), ...
        max(self_per_performer), min(other_per_performer));

fprintf('\nExhaustive label-permutation test (6! = %d permutations):\n', n_perm);
fprintf('  p (one-sided, self < other) = %.4f\n', p_perm_one);
fprintf('  p (two-sided)               = %.4f\n', p_perm_two);

fprintf('\nDescriptive paired t-test:\n');
fprintf('  t(%d) = %.3f, p = %.4f\n', stats_ttest.df, stats_ttest.tstat, p_ttest);
fprintf('  95%% CI of difference: [%.4f, %.4f]\n', ci_ttest(1), ci_ttest(2));
fprintf("  Cohen's dz = %.3f\n", cohens_d);

fprintf('\nPer-performer normalised DTW distance:\n');
fprintf('%-12s  Self     Other\n', 'Performer');
for pf = 1:n_performer
    fprintf('%-12s  %.4f   %.4f\n', performers{pf}, ...
            self_per_performer(pf), other_per_performer(pf));
end

%% ---- Figure ------------------------------------------------------------
figure('Units','centimeters', 'Position',[2 2 12 10], ...
       'PaperUnits','centimeters', 'PaperPosition',[0 0 12 10]);
hold on;

for pf = 1:n_performer
    plot([1 2], [self_per_performer(pf), other_per_performer(pf)], ...
         '-o', ...
         'Color',           performer_colors(pf,:), ...
         'MarkerFaceColor', performer_colors(pf,:), ...
         'MarkerEdgeColor', 'k', ...
         'LineWidth', 2, ...
         'MarkerSize', 10, ...
         'DisplayName', performers{pf});
end

xlim([0.5, 2.5]);
ylim([0,   0.4]);
set(gca, 'XTick', [1 2], 'XTickLabel', {'Self','Other'}, ...
         'FontName', 'Arial', 'FontSize', 10);
ylabel('Normalised DTW distance', 'FontName','Arial', 'FontSize',12);
legend('Location','northwest');
box on;
title(sprintf('Self vs. Other  (p_{perm} = %.4f)', p_perm_one), ...
      'FontName','Arial', 'FontSize',11);
