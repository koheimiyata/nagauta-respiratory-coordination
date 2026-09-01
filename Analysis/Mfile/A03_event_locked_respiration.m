%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  A03_event_locked_respiration
%
%  Event-locked respiration analysis with surrogate envelope
%
%  Project : Nagauta respiratory coordination
%
%  Description
%  -----------
%  For each event type in MusicalEvent.mat (default: Decel, Komi, Section),
%  this script:
%    1. Z-scores the full recording for each performer and take, then
%       extracts segments in a [-w1, +w2] second window around each
%       event onset.
%    2. Computes the grand-mean waveform pooled across performers and events.
%    3. Builds a 95% pointwise surrogate envelope using random event
%       times that avoid data boundaries and a guard region around
%       the original events.
%    4. Plots the observed waveform against the surrogate envelope,
%       with performer-wise averages shown as thin colored lines.
%
%  A separate figure is produced for each event type automatically.
%  To add or remove event types, edit the 'event_types' variable below.
%
%  Folder structure assumed
%  ------------------------
%  <project_root>/
%    Mfile/        <- this file lives here
%    Mat_file/
%      Hexoskin_data.mat (respiration, ECG, acceleration per performer and take)
%      MusicalEvent.mat  (event onset times in ms; struct with fields
%                         Decel/Komi/Section, each with sub-fields S1, S2)
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
lowpass_cutoff_hz = 5;    % low-pass filter cutoff (Hz)

pre_event_s    = 10;        % window before event onset (seconds)
post_event_s   = 10;        % window after  event onset (seconds)
n_surrogate    = 1000;      % number of surrogate iterations
guard_region_s = 10;        % exclude ± this many seconds around each event
                            %   when sampling surrogate times

downsample_factor = fs_raw / fs_ds;
pre_event_samp    = pre_event_s  * fs_ds;
post_event_samp   = post_event_s * fs_ds;
guard_region_samp = round(guard_region_s * fs_ds);
n_time            = pre_event_samp + post_event_samp + 1;
time_axis         = linspace(-pre_event_s, post_event_s, n_time);
ms_per_sample     = 1000 / fs_ds;

%% ---- Experiment metadata -----------------------------------------------
take_field  = {'S1', 'S2'};
performers  = {'Uta','Shamisen1','Kotsuzumi','Fue','Taiko','Shamisen2'};
event_types = {'Decel', 'Komi', 'Section'};   % edit to add/remove event types

n_take      = numel(take_field);
n_performer = numel(performers);

performer_colors = [
    0.00, 0.45, 0.74;   % Uta
    0.85, 0.33, 0.10;   % Shamisen1
    0.47, 0.68, 0.19;   % Kotsuzumi
    0.49, 0.18, 0.56;   % Fue
    0.30, 0.75, 0.93;   % Taiko
    0.64, 0.08, 0.18];  % Shamisen2

%% ---- Load data ---------------------------------------------------------
load(fullfile(data_dir, 'Hexoskin_data.mat'));  % struct 'Nagauta'
load(fullfile(data_dir, 'MusicalEvent.mat'));   % struct 'Event'

%% ---- Pre-process respiration -------------------------------------------
%  resp_ds{take, performer} = z-scored, downsampled signal (full recording)

resp_ds  = cell(n_take, n_performer);
take_len = zeros(n_take, 1);

for tk = 1:n_take
    for pf = 1:n_performer
        raw = Nagauta.(take_field{tk}).(performers{pf}).respiration_thoracic;
        raw = raw(:);
        sig = lowpass(raw, lowpass_cutoff_hz, fs_raw);
        sig = downsample(sig, downsample_factor);
        resp_ds{tk, pf} = zscore(sig);
    end
    take_len(tk) = numel(resp_ds{tk, 1});
end

%% ---- Loop over event types ---------------------------------------------
for ev = 1:numel(event_types)

    ev_type = event_types{ev};

    if ~isfield(Event, ev_type)
        fprintf('Event type "%s" not found in MusicalEvent.mat — skipping.\n', ev_type);
        continue;
    end

    fprintf('\nProcessing event type: %s\n', ev_type);

    %% -- Collect event onset indices (samples at fs_ds) per take ---------
    event_idx_take = cell(n_take, 1);
    for tk = 1:n_take
        if isfield(Event.(ev_type), take_field{tk})
            onset_ms = Event.(ev_type).(take_field{tk});
            event_idx_take{tk} = round(onset_ms(:) / ms_per_sample);
        else
            event_idx_take{tk} = [];
        end
    end

    %% -- Extract observed event-locked segments --------------------------
    trials_all          = [];
    trials_by_performer = cell(n_performer, 1);

    for tk = 1:n_take
        ev_idx = event_idx_take{tk};

        for pf = 1:n_performer
            sig = resp_ds{tk, pf};

            for k = 1:numel(ev_idx)
                c = ev_idx(k);
                if (c - pre_event_samp) < 1 || (c + post_event_samp) > numel(sig)
                    continue;
                end
                seg = sig((c - pre_event_samp):(c + post_event_samp));
                trials_all              = [trials_all, seg];
                trials_by_performer{pf} = [trials_by_performer{pf}, seg];
            end
        end
    end

    grand_mean = mean(trials_all, 2, 'omitnan');

    mean_by_performer = nan(n_time, n_performer);
    for pf = 1:n_performer
        if ~isempty(trials_by_performer{pf})
            mean_by_performer(:, pf) = mean(trials_by_performer{pf}, 2, 'omitnan');
        end
    end

    n_trials_total = size(trials_all, 2);
    fprintf('  Observed trials (all performers pooled): %d\n', n_trials_total);

    %% -- Surrogate distribution ------------------------------------------
    surr_grand_mean = nan(n_time, n_surrogate);

    for b = 1:n_surrogate
        surr_trials = [];

        for tk = 1:n_take
            L      = take_len(tk);
            ev_idx = event_idx_take{tk};
            n_ev   = numel(ev_idx);
            if n_ev == 0, continue; end

            valid = false(L, 1);
            valid((1 + pre_event_samp):(L - post_event_samp)) = true;

            for k = 1:n_ev
                c = ev_idx(k);
                valid(max(1,c-guard_region_samp) : min(L,c+guard_region_samp)) = false;
            end

            valid_idx = find(valid);
            if numel(valid_idx) < n_ev
                error(['Not enough valid surrogate positions for event type "%s", ' ...
                       'take %s.\nReduce guard_region_s or window length.'], ...
                       ev_type, take_field{tk});
            end

            surr_idx = valid_idx(randperm(numel(valid_idx), n_ev));

            for pf = 1:n_performer
                sig = resp_ds{tk, pf};
                for k = 1:n_ev
                    c   = surr_idx(k);
                    seg = sig((c - pre_event_samp):(c + post_event_samp));
                    surr_trials = [surr_trials, seg];
                end
            end
        end

        if ~isempty(surr_trials)
            surr_grand_mean(:, b) = mean(surr_trials, 2, 'omitnan');
        end
    end

    surr_env_lo = prctile(surr_grand_mean,  2.5, 2);
    surr_env_hi = prctile(surr_grand_mean, 97.5, 2);

    %% -- Quantitative summary --------------------------------------------
    % Pre-event baseline: mean of grand mean in [-w1, -1] s
    % Post-event response: mean in [0, +w2] s
    pre_mask  = time_axis >= -pre_event_s  & time_axis < 0;
    post_mask = time_axis >= 0             & time_axis <= post_event_s;

    pre_mean  = mean(grand_mean(pre_mask));
    post_mean = mean(grand_mean(post_mask));

    % Peak (signed maximum absolute deviation from zero) in post window
    [peak_val, peak_idx_rel] = max(abs(grand_mean(post_mask)));
    post_time_vec = time_axis(post_mask);
    peak_time     = post_time_vec(peak_idx_rel);
    peak_signed   = grand_mean(post_mask);
    peak_signed   = peak_signed(peak_idx_rel);

    % Proportion of time points outside the 95% surrogate envelope
    outside_mask = grand_mean > surr_env_hi | grand_mean < surr_env_lo;
    pct_outside_total = 100 * mean(outside_mask);
    pct_outside_post  = 100 * mean(outside_mask(post_mask));

    % Per-performer: peak in post-event window (mean across performers)
    performer_peaks = nan(n_performer, 1);
    for pf = 1:n_performer
        if ~all(isnan(mean_by_performer(:, pf)))
            post_sig = mean_by_performer(post_mask, pf);
            [~, pi]  = max(abs(post_sig));
            performer_peaks(pf) = post_sig(pi);
        end
    end

    fprintf('\n--- Quantitative summary: %s ---\n', ev_type);
    fprintf('Sample unit: z-scored thoracic respiration, n = %d trials\n', n_trials_total);
    fprintf('  (pooled across %d performers x %d takes)\n\n', n_performer, n_take);
    fprintf('Grand-mean respiration:\n');
    fprintf('  Pre-event mean  (%.0f to 0 s):  %.4f\n', -pre_event_s, pre_mean);
    fprintf('  Post-event mean (0 to +%.0f s): %.4f\n', post_event_s, post_mean);
    fprintf('  Change (post − pre): %.4f\n\n', post_mean - pre_mean);
    fprintf('Peak response (post-event window):\n');
    fprintf('  Peak amplitude: %.4f (at t = %.2f s)\n', peak_signed, peak_time);
    fprintf('  |Peak| = %.4f\n\n', peak_val);
    fprintf('Surrogate envelope (95%% pointwise):\n');
    fprintf('  %% time points outside envelope — full window: %.1f%%\n', pct_outside_total);
    fprintf('  %% time points outside envelope — post-event:  %.1f%%\n\n', pct_outside_post);
    fprintf('Per-performer signed peak in post-event window:\n');
    fprintf('%-12s  Peak z\n', 'Performer');
    for pf = 1:n_performer
        if ~isnan(performer_peaks(pf))
            fprintf('%-12s  %+.4f\n', performers{pf}, performer_peaks(pf));
        end
    end

    %% -- Figure ----------------------------------------------------------
    figure('Units','centimeters', 'Position',[2 2 17 10], ...
           'PaperUnits','centimeters', 'PaperPosition',[0 0 17 10]);
    hold on;

    fill([time_axis, fliplr(time_axis)], ...
         [surr_env_hi', fliplr(surr_env_lo')], ...
         [0.7 0.7 0.7], 'FaceAlpha',0.6, 'EdgeColor','none', ...
         'DisplayName','95% surrogate envelope');

    for pf = 1:n_performer
        if all(isnan(mean_by_performer(:, pf))), continue; end
        plot(time_axis, mean_by_performer(:, pf), ...
             'Color', performer_colors(pf,:), 'LineWidth',1, ...
             'DisplayName', performers{pf});
    end

    plot(time_axis, grand_mean, 'k', 'LineWidth',3, ...
         'DisplayName','Grand mean (observed)');

    xline(0, 'k--', 'LineWidth',1.25, 'HandleVisibility','off');

    xlim([-pre_event_s, post_event_s]);
    ylim([-1, 2]);
    set(gca, 'FontName','Arial', 'FontSize',10);
    xlabel('Time relative to event onset (s)',  'FontName','Arial', 'FontSize',12);
    ylabel('Respiration (z-scored, thoracic)',  'FontName','Arial', 'FontSize',12);
    title(sprintf('Event-locked respiration: %s  (n = %d trials)', ev_type, n_trials_total), ...
          'FontName','Arial', 'FontSize',11);
    grid on; box on; axis square;

    fprintf('  Surrogates: B = %d, guard region = ±%.0f s\n', n_surrogate, guard_region_s);

end % event type loop
