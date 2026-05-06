function analyze_DOPU_PhR_scatter_stats(data_dir, save_dir, target_Z, data_dir2, method1_name, method2_name)
% ANALYZE_DOPU_PHR_SCATTER_STATS  Clean correlation analysis: SS-DOPU vs PhR
%
% Usage:
%   analyze_DOPU_PhR_scatter_stats(data_dir, save_dir, target_Z)
%   analyze_DOPU_PhR_scatter_stats(data_dir, save_dir, target_Z, data_dir2)
%   analyze_DOPU_PhR_scatter_stats(data_dir, save_dir, target_Z, data_dir2, method1_name, method2_name)
%
%   - data_dir      (optional) : folder containing *_PhR_stats.mat and *_Metadata_stats.mat for Method 1
%   - save_dir      (optional) : output folder for PNG
%   - target_Z      (optional) : depth/frame index to analyze (default: 30)
%   - data_dir2     (optional) : folder for Method 2 data (if empty or omitted, single-method mode)
%   - method1_name  (optional) : display name for Method 1 (default: 'Method 1')
%   - method2_name  (optional) : display name for Method 2 (default: 'Method 2')
%
% Purpose:
% - Demonstrates the negative correlation between SD-DOPU and Phase Retardation
% - Shows that lower DOPU leads to artificially elevated PhR (noise-bias effect)
% - Produces a clean, publication-ready scatter plot with regression lines
% - **Supports dual-method comparison** with two colors + legend (backward compatible)
%
% Output file: <PhR_filename>_frame_<Z>_DOPU_PhR_Correlation.png

if nargin < 1 || isempty(data_dir), data_dir = 'D:\1-Liu Jian\yongxin.wang\Output\Graph\Mat'; end
if nargin < 2 || isempty(save_dir), save_dir = 'D:\1-Liu Jian\yongxin.wang\Output\Graph\Output'; end
if nargin < 3 || isempty(target_Z), target_Z = 50; end  % 与 tiff_frame=10 保持一致
if nargin < 4, data_dir2 = ''; end
% 方法名称 (用于图例和标注)
% 默认值：与全局稳定性脚本保持同步
% 方法1：DDG
% 方法2：SD-DOPU + Bayesian (之前为 SSDOPU++)
if nargin < 5 || isempty(method1_name), method1_name = 'DDG'; else
    % 根据调用者覆盖，generic 名称下方会自动发现真实名称
end
if nargin < 6 || isempty(method2_name), method2_name = 'SD-DOPU + Bayesian'; else
    % allow caller to override; fall back to auto-discovery below if generic
end

% === 单位换算参数 ===
phr_scale = 180/pi;  % rad → °（直接换算，不做范围截断）

if ~exist(save_dir, 'dir'), mkdir(save_dir); end

required_patterns = {'*_PhR_stats.mat', '*_Metadata_stats.mat'};

% Always treat the provided path as the ROOT and auto-discover method folders recursively.
if exist(data_dir, 'dir')
    discovered = auto_discover_method_dirs(data_dir, required_patterns);
    if isempty(discovered)
        error('No valid data folders found under root: %s', data_dir);
    end

    data_dir = discovered(1).data_dir;
    if strcmp(method1_name, 'Method 1')
        method1_name = discovered(1).method_name;
    end

    if isempty(data_dir2) && numel(discovered) >= 2
        data_dir2 = discovered(2).data_dir;
        if strcmp(method2_name, 'Method 2')
            method2_name = discovered(2).method_name;
        end
    end

    fprintf('Auto-discovered Method 1 dir: %s\n', data_dir);
    if ~isempty(data_dir2)
        fprintf('Auto-discovered Method 2 dir: %s\n', data_dir2);
    end
else
    error('Specified root directory does not exist: %s', data_dir);
end

if ~isempty(data_dir2) && exist(data_dir2, 'dir') && ~has_required_files_local(data_dir2, required_patterns)
    discovered2 = auto_discover_method_dirs(data_dir2, required_patterns);
    if ~isempty(discovered2)
        data_dir2 = discovered2(1).data_dir;
        if strcmp(method2_name, 'Method 2')
            method2_name = discovered2(1).method_name;
        end
        fprintf('Resolved Method 2 root to data dir: %s\n', data_dir2);
    else
        error('Method 2 path does not contain valid data folders: %s', data_dir2);
    end
end

has_method2 = ~isempty(data_dir2) && exist(data_dir2, 'dir');

% -------------------- 1) Load files — Method 1 --------------------
[vec_DOPU_1, vec_PhR_1, phR_files_1] = load_dopu_phr_local(data_dir, target_Z, phr_scale);
N_1 = numel(vec_DOPU_1);

% -------------------- 1b) Load files — Method 2 (if provided) --------------------
if has_method2
    [vec_DOPU_2, vec_PhR_2, ~] = load_dopu_phr_local(data_dir2, target_Z, phr_scale);
    N_2 = numel(vec_DOPU_2);
else
    vec_DOPU_2 = []; vec_PhR_2 = []; N_2 = 0;
end

% -------------------- 2) Statistical analysis — Method 1 --------------------
[r_1, p_val_1, fit_1, k_1, b_1, R2_1] = compute_corr_stats_local(vec_DOPU_1, vec_PhR_1);

% -------------------- 2b) Statistical analysis — Method 2 --------------------
if has_method2
    [r_2, p_val_2, fit_2, k_2, b_2, R2_2] = compute_corr_stats_local(vec_DOPU_2, vec_PhR_2);
end

% -------------------- 3) Plot --------------------
fig = figure('Color','w','Name','SS-DOPU vs PhR Correlation', 'Position', [100 100 700 600]);
hold on;

% Determine Y-range using unified percentile/margin logic
% This matches the logic used in analyze_phr_global_stability.m so
% that both plots share the same PhR bounds when using identical data.
if has_method2
    all_phr = [vec_PhR_1; vec_PhR_2];
else
    all_phr = vec_PhR_1;
end
[Ymin, Ymax] = compute_pr_axis_limits(all_phr);

% ---- Color scheme ----
% Method 1: Red tones
lineColor_1  = [0.85, 0.05, 0.05];     % strong red
pointColor_1 = [1.00, 0.65, 0.65];     % pale red
% Method 2: Blue tones
lineColor_2  = [0.05, 0.30, 0.85];     % strong blue
pointColor_2 = [0.65, 0.78, 1.00];     % pale blue

markerSize = 3;
markerAlpha = 0.45;
maxScatter = 5000;

% --- Scatter + regression: Method 1 ---
hLine_1 = plot_scatter_fit_local(vec_DOPU_1, vec_PhR_1, pointColor_1, lineColor_1, ...
    markerSize, markerAlpha, maxScatter, fit_1);

% --- Scatter + regression: Method 2 ---
if has_method2
    hLine_2 = plot_scatter_fit_local(vec_DOPU_2, vec_PhR_2, pointColor_2, lineColor_2, ...
        markerSize, markerAlpha, maxScatter, fit_2);
end

% Axis labels and title
xlabel('SD-DOPU (Polarization Uniformity)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Phase Retardation (°/μm)', 'FontSize', 20, 'FontWeight', 'bold');
if has_method2
    title('Correlation: SS-DOPU vs PhR', 'FontSize', 22, 'FontWeight', 'bold');
else
    title('Correlation Analysis: SS-DOPU vs Phase Retardation', 'FontSize', 22, 'FontWeight', 'bold');
end

set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 2.0, 'Box', 'on', ...
    'TickDir', 'in', 'TickLength', [0.015 0.015]);
xlim([0 1]);
ylim([Ymin, Ymax]);

% ---- Legend (with correlation coefficients) ----
if has_method2
    % put r values inside legend entries instead of separate textbox
    legend([hLine_1, hLine_2], ...
        {sprintf('%s (r=%.3f)', method1_name, r_1), ...
         sprintf('%s (r=%.3f)', method2_name, r_2)}, ...
        'Location', 'northeast', 'FontSize', 16, 'FontWeight', 'bold', 'Box', 'off');
else
    % single method: show r in legend as well for consistency
    legend(hLine_1, sprintf('%s (r=%.3f)', method1_name, r_1), ...
        'Location', 'northeast', 'FontSize', 16, 'FontWeight', 'bold', 'Box', 'off');
end

% remove annotation box, correlation now part of legend (cleaner layout)
% previous annotation code deleted

% ---- Save figure ----
[~,bn,~] = fileparts(phR_files_1(1).name);
bn = regexprep(bn, '_PhR_stats$', '');
outname = fullfile(save_dir, sprintf('%s_frame_%02d_DOPU_PhR_Correlation.png', bn, target_Z));
set(fig, 'PaperPositionMode', 'auto');
print(fig, outname, '-dpng', '-r300');

% ---- Console output ----
fprintf('\n======== SS-DOPU vs PhR Correlation Analysis ========\n');
fprintf('Frame: %d\n', target_Z);
fprintf('\n--- %s ---\n', method1_name);
fprintf('Total pixels: %d\n', N_1);
fprintf('  Pearson r = %.4f (p = %.2e)\n', r_1, p_val_1);
fprintf('  Linear fit: PhR = %.4f * DOPU + %.4f  (°/μm)\n', k_1, b_1);
fprintf('  R^2 = %.4f\n', R2_1);
print_interpretation_local(r_1);

if has_method2
    fprintf('\n--- %s ---\n', method2_name);
    fprintf('Total pixels: %d\n', N_2);
    fprintf('  Pearson r = %.4f (p = %.2e)\n', r_2, p_val_2);
    fprintf('  Linear fit: PhR = %.4f * DOPU + %.4f  (°/μm)\n', k_2, b_2);
    fprintf('  R^2 = %.4f\n', R2_2);
    print_interpretation_local(r_2);
end

fprintf('\nOutput saved: %s\n', outname);
fprintf('=====================================================\n\n');

end

%% ========================== Helper Functions ==========================

% shared percentile/margin routine used by multiple scripts
function [low, high] = compute_pr_axis_limits(values)
    p01 = prctile(values, 0.5);
    p99 = prctile(values, 99.5);
    margin = (p99 - p01) * 0.1;
    low = max(0, p01 - margin);
    high = p99 + margin;
end

function [vec_DOPU, vec_PhR, phR_files] = load_dopu_phr_local(data_dir, target_Z, phr_scale)
% Load PhR and DOPU data, convert rad → degrees, no upper clamp.
    if nargin < 3, phr_scale = 180/pi; end
    phR_files = dir(fullfile(data_dir, '*_PhR_stats.mat'));
    meta_files = dir(fullfile(data_dir, '*_Metadata_stats.mat'));
    if isempty(phR_files)
        error('No *_PhR_stats.mat found in %s', data_dir);
    end
    if isempty(meta_files)
        error('No *_Metadata_stats.mat found in %s', data_dir);
    end

    phR_path = fullfile(data_dir, phR_files(1).name);
    meta_path = fullfile(data_dir, meta_files(1).name);

    S1 = load(phR_path);
    S2 = load(meta_path);

    % Extract PhR
    if isfield(S1, 'PhR_phys_flat')
        PhR_all = S1.PhR_phys_flat;
    elseif isfield(S1, 'PhR_enface_avg')
        PhR_all = S1.PhR_enface_avg;
    else
        error('PhR variable not found in %s', phR_path);
    end

    % Extract DOPU
    if isfield(S2, 'dopu_splitSpectrum')
        DOPU_all = S2.dopu_splitSpectrum;
    elseif isfield(S2, 'DOPU')
        DOPU_all = S2.DOPU;
    else
        error('dopu_splitSpectrum variable not found in %s', meta_path);
    end

    % Extract target frame and convert PhR: rad → °/μm, then clamp
    PhR_frame = extract_frame_local(PhR_all, target_Z, 'PhR') * phr_scale;  % rad → °
    PhR_frame = max(PhR_frame, 0);  % 仅截断负値
    DOPU_frame = extract_frame_local(DOPU_all, target_Z, 'DOPU');

    % Reconcile sizes
    [PhR_frame, DOPU_frame] = reconcile_sizes_local(PhR_frame, DOPU_frame);

    % Mask
    mask = (PhR_frame > 0.0001) & isfinite(PhR_frame) & isfinite(DOPU_frame);
    vec_DOPU = DOPU_frame(mask);
    vec_PhR  = PhR_frame(mask);

    valid = isfinite(vec_DOPU) & isfinite(vec_PhR);
    vec_DOPU = vec_DOPU(valid);
    vec_PhR  = vec_PhR(valid);
end

function frame = extract_frame_local(data, target_Z, name)
    if ndims(data) == 3
        % 数据维度为 [nZ, nX, nY]，dim1 是深度 Z
        if size(data,1) < target_Z
            error('target_Z (%d) > available frames (%d) in %s data', target_Z, size(data,1), name);
        end
        frame = double(squeeze(data(target_Z, :, :)));  % [nX, nY]，En-face 切片
    elseif ismatrix(data)
        frame = double(data);
    else
        error('Unexpected %s array shape.', name);
    end
end

function [A, B] = reconcile_sizes_local(A, B)
    if isequal(size(A), size(B)), return; end
    [hA,wA] = size(A);
    [hB,wB] = size(B);
    if wA == wB && hB > hA
        B = B(1:hA, :);
        warning('Size mismatch resolved by cropping B to %dx%d.', hA, wA);
    elseif wA == wB && hA > hB
        A = A(1:hB, :);
        warning('Size mismatch resolved by cropping A to %dx%d.', hB, wA);
    elseif hA == hB && wB > wA
        B = B(:,1:wA);
        warning('Size mismatch resolved by cropping B columns to %d.', wA);
    elseif hA == hB && wA > wB
        A = A(:,1:wB);
        warning('Size mismatch resolved by cropping A columns to %d.', wB);
    else
        error('Size mismatch: A=%s, B=%s — automatic reconciliation not possible.', mat2str(size(A)), mat2str(size(B)));
    end
end

function [r, p_val, fit_coeff, k, b, R2] = compute_corr_stats_local(vec_DOPU, vec_PhR)
    try
        [r, p_val] = corr(vec_DOPU(:), vec_PhR(:), 'Type', 'Pearson', 'Rows', 'complete');
    catch
        C = corrcoef(vec_DOPU(:), vec_PhR(:));
        r = C(1,2);
        p_val = NaN;
    end
    fit_coeff = polyfit(vec_DOPU(:), vec_PhR(:), 1);
    k = fit_coeff(1);
    b = fit_coeff(2);
    R2 = r^2;
end

function hLine = plot_scatter_fit_local(vec_DOPU, vec_PhR, pointColor, lineColor, markerSize, markerAlpha, maxScatter, fit_coeff)
    N = numel(vec_DOPU);
    if N > maxScatter
        idx = randperm(N, maxScatter);
        scatter(vec_DOPU(idx), vec_PhR(idx), markerSize, pointColor, 'filled', ...
            'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeAlpha', 0);
    else
        scatter(vec_DOPU, vec_PhR, markerSize, pointColor, 'filled', ...
            'MarkerFaceAlpha', markerAlpha, 'MarkerEdgeAlpha', 0);
    end
    xf = linspace(0, 1, 400);
    hLine = plot(xf, polyval(fit_coeff, xf), 'Color', lineColor, 'LineWidth', 2.5);
end

function print_interpretation_local(r)
    fprintf('Interpretation:\n');
    if r < -0.05
        fprintf('  -> Negative correlation: SD-DOPU down leads to PhR up\n');
        fprintf('  -> Supports noise-bias hypothesis\n');
    elseif r < 0
        fprintf('  -> Weak negative trend detected\n');
    else
        fprintf('  -> No significant negative correlation\n');
    end
end

function tf = has_required_files_local(dir_path, required_patterns)
    tf = true;
    for i = 1:numel(required_patterns)
        if isempty(dir(fullfile(dir_path, required_patterns{i})))
            tf = false;
            return;
        end
    end
end





