% analyze_BScan_OA_ROI_PhR_Ascan_stats.m
% Independent B-Scan OA ROI agreement + PhR A-scan comparison (two methods)
% Notes:
% - B-Scan data size is [Z, X]
% - ROI mask size is [Z, X]
% - A-scan is fixed X, sampled along Z

clear; close all; clc;

% Ensure helper functions are on the MATLAB path so quColoring and others
% can be found. Try common relative locations (current working dir and
% script-relative ../function).
try
    addpath(genpath(fullfile(pwd, 'function')));
catch
end
try
    scriptPath = mfilename('fullpath');
    if ~isempty(scriptPath)
        scriptDir = fileparts(scriptPath);
        addpath(genpath(fullfile(scriptDir, '..', 'function')));
    else
        scriptFile = which('analyze_BScan_OA_ROI_PhR_Ascan_stats.m');
        if ~isempty(scriptFile)
            scriptDir = fileparts(scriptFile);
            addpath(genpath(fullfile(scriptDir, '..', 'function')));
        end
    end
catch
end
run_tic = tic;
fprintf('[%s] analyze_BScan_OA_ROI_PhR_Ascan_stats started\n', datestr(now, 'HH:MM:SS'));

% Use Times New Roman consistently for all figure text.
set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultLegendFontName', 'Times New Roman');

%% ==================== 1. Config ====================
input_root_dir = 'G:\1-Project\2023 王永鑫\Data\Rep\Output';      % root directory that contains both methods
method1_dir = '';            % optional: direct bscan_stats_mat folder
method2_dir = '';            % optional: direct bscan_stats_mat folder
output_dir = 'G:\1-Project\2023 王永鑫\Data\Rep\Graph_res';

method1_name = 'DDG';
method2_name = 'Bayesian optimization';

target_bscan_index = 35;     % IMPORTANT: B-Scan index, not Z-layer

oa_source = 'cum';           % 'cum' or 'local'
roi_mask_path = '';
roi_annotation_image = fullfile(input_root_dir, '光轴ROI区域.png');
roi_mode = 'from_annotation_image';   % 'load' / 'draw' / 'load_or_draw' / 'auto_default' / 'from_annotation_image'

ascan_annotation_image = fullfile(input_root_dir, '相位延迟要分析的Ascan.png');
ascan_mode = 'from_annotation_image'; % 'interactive' / 'from_annotation_image'

ascan_start_mode = 'boundary';   % 'boundary' / 'full'
ascan_depth_axis = 'relative';   % 'relative' / 'absolute'
ascan_axial_resolution_um_per_px = 3.0; % axial sampling (um/pixel) for A-scan plots

ascan_x = [];                % [] for interactive click; otherwise column index
ascan_radius_px = 0;         % 0 for single column; >0 for neighborhood averaging
phr_scale = 180/pi;          % rad -> deg
phr_display_range = [];      % [] = auto range from data (deg)
phr_colorbar_range = [0 5];  % fixed colorbar range (deg)
phr_colorbar_ticks = 0:1:5;  % fixed ticks for colorbar
oa_local_window_radius_px = 4;  % 9x9 window for local OA consistency

oa_mag_threshold = 1e-6;     % threshold for OA vector magnitude to reject low-signal pixels
use_tissue_mask = false;     % derive tissue mask from metadata if available

%% ==================== 2. Resolve method folders ====================
if isempty(output_dir)
    error('output_dir is empty.');
end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

[method1_dir, method2_dir] = resolve_method_dirs(input_root_dir, method1_dir, method2_dir);

%% ==================== 3. Load method bundles and extract target B-Scan ====================
fprintf('[%s] Loading Method 1 bundle...\n', datestr(now, 'HH:MM:SS'));
bundle1 = load_bscan_bundle(method1_dir, target_bscan_index, oa_source, phr_scale);
fprintf('[%s] Loading Method 2 bundle...\n', datestr(now, 'HH:MM:SS'));
bundle2 = load_bscan_bundle(method2_dir, target_bscan_index, oa_source, phr_scale);

fprintf('Method 1 dir: %s\n', method1_dir);
fprintf('Method 2 dir: %s\n', method2_dir);
fprintf('Target B-Scan index: %d\n', target_bscan_index);
fprintf('Method 1 files:\n  Metadata: %s\n  OA: %s\n  PhR: %s\n', bundle1.metadata_file, bundle1.oa_file, bundle1.phr_file);
fprintf('Method 2 files:\n  Metadata: %s\n  OA: %s\n  PhR: %s\n', bundle2.metadata_file, bundle2.oa_file, bundle2.phr_file);
fprintf('[%s] Bundles loaded, elapsed %.1f s\n', datestr(now, 'HH:MM:SS'), toc(run_tic));

OA1_deg = bundle1.OA_deg;
OA2_deg = bundle2.OA_deg;
PhR1_unit = bundle1.PhR_unit;
PhR2_unit = bundle2.PhR_unit;
PhR1_rad = bundle1.PhR_rad;
PhR2_rad = bundle2.PhR_rad;
base_img = bundle1.base_img;
topLine_bscan = bundle1.topLine_bscan;

if isempty(OA1_deg) || isempty(OA2_deg)
    error('OA slice extraction failed for target B-Scan %d.', target_bscan_index);
end
if isempty(PhR1_unit) || isempty(PhR2_unit)
    error('PhR slice extraction failed for target B-Scan %d.', target_bscan_index);
end

phr_display_range = resolve_phr_display_range(PhR1_unit, PhR2_unit, phr_display_range);
fprintf('PhR display range: [%.2f, %.2f] deg\n', phr_display_range(1), phr_display_range(2));

%% ==================== 4. ROI mask (same for both methods) ====================
fprintf('[%s] Preparing ROI mask...\n', datestr(now, 'HH:MM:SS'));
roi_mask = load_or_draw_roi_mask(roi_mask_path, roi_mode, base_img, output_dir, target_bscan_index, roi_annotation_image, size(OA1_deg));
if ~isequal(size(roi_mask), size(OA1_deg))
    error('ROI mask size mismatch. mask=[%d,%d], OA=[%d,%d]', size(roi_mask,1), size(roi_mask,2), size(OA1_deg,1), size(OA1_deg,2));
end

valid_mask = build_valid_mask(bundle1, bundle2, OA1_deg, OA2_deg, oa_mag_threshold, use_tissue_mask, target_bscan_index);
annotation_roi_mask = roi_mask;
stats_mask = annotation_roi_mask & valid_mask;

%% ==================== 5. OA local consistency (per method) ====================
fprintf('[%s] Computing OA local consistency...\n', datestr(now, 'HH:MM:SS'));
[oa_stats_m1, oa_local_m1] = compute_oa_local_consistency(OA1_deg, stats_mask, oa_local_window_radius_px);
[oa_stats_m2, oa_local_m2] = compute_oa_local_consistency(OA2_deg, stats_mask, oa_local_window_radius_px);
oa_stats_m1.method_name = method1_name;
oa_stats_m2.method_name = method2_name;
oa_neighbor_m1 = compute_neighbor_angle_change(OA1_deg, stats_mask);
oa_neighbor_m2 = compute_neighbor_angle_change(OA2_deg, stats_mask);

%% ==================== 6. OA figures ====================
outputs = struct();
outputs.oa_overlay_compare = fullfile(output_dir, sprintf('BScan%04d_OA_ROI_Overlay_Compare.png', target_bscan_index));
outputs.oa_crop_compare = fullfile(output_dir, sprintf('BScan%04d_OA_ROI_Crop_Compare.png', target_bscan_index));
outputs.oa_local_consistency = fullfile(output_dir, sprintf('BScan%04d_OA_LocalConsistency.png', target_bscan_index));
outputs.oa_local_consistency_legend_density = fullfile(output_dir, sprintf('BScan%04d_OA_LocalConsistency_LegendDensity.png', target_bscan_index));
outputs.oa_local_consistency_legend_metrics = fullfile(output_dir, sprintf('BScan%04d_OA_LocalConsistency_LegendMetrics.png', target_bscan_index));
outputs.oa_colorbar = fullfile(output_dir, sprintf('BScan%04d_OA_Colorbar.png', target_bscan_index));
outputs.oa_overlay_m1 = fullfile(output_dir, sprintf('BScan%04d_OA_ROI_Overlay_%s.png', target_bscan_index, safe_name(method1_name)));
outputs.oa_overlay_m2 = fullfile(output_dir, sprintf('BScan%04d_OA_ROI_Overlay_%s.png', target_bscan_index, safe_name(method2_name)));
outputs.oa_crop_m1 = fullfile(output_dir, sprintf('BScan%04d_OA_ROI_Crop_%s.png', target_bscan_index, safe_name(method1_name)));
outputs.oa_crop_m2 = fullfile(output_dir, sprintf('BScan%04d_OA_ROI_Crop_%s.png', target_bscan_index, safe_name(method2_name)));

plot_oa_roi_overlay_compare(bundle1.OA_rgb, bundle2.OA_rgb, annotation_roi_mask, method1_name, method2_name, outputs.oa_overlay_compare, outputs.oa_overlay_m1, outputs.oa_overlay_m2);
plot_oa_roi_crop_compare(bundle1.OA_rgb, bundle2.OA_rgb, annotation_roi_mask, method1_name, method2_name, outputs.oa_crop_compare, outputs.oa_crop_m1, outputs.oa_crop_m2);
plot_oa_local_consistency_figure(oa_local_m1.local_std_vals, oa_local_m2.local_std_vals, oa_stats_m1, oa_stats_m2, oa_neighbor_m1, oa_neighbor_m2, method1_name, method2_name, outputs.oa_local_consistency, outputs.oa_local_consistency_legend_density, outputs.oa_local_consistency_legend_metrics);
export_colorbar_only(outputs.oa_colorbar, hsv(256), [0 180], 0:45:180, '°');

%% ==================== 7. PhR A-scan selection ====================
fprintf('[%s] Resolving A-scan X...\n', datestr(now, 'HH:MM:SS'));
ascan_x = resolve_ascan_x(ascan_x, base_img, output_dir, target_bscan_index, ascan_mode, ascan_annotation_image, size(PhR1_unit));

%% ==================== 8. PhR A-scan stats ====================
fprintf('[%s] Computing PhR A-scan statistics...\n', datestr(now, 'HH:MM:SS'));
[ascan_stats, ascan_table] = compute_phr_ascan_stats_from_boundary(PhR1_unit, PhR2_unit, bundle1.topLine_bscan, bundle2.topLine_bscan, ascan_x, ascan_radius_px, method1_name, method2_name, ascan_start_mode, ascan_depth_axis);

%% ==================== 9. PhR A-scan figures ====================
outputs.ascan_plot = fullfile(output_dir, sprintf('BScan%04d_PhR_Ascan_BoundaryAligned.png', target_bscan_index));
outputs.ascan_plot_legend = fullfile(output_dir, sprintf('BScan%04d_PhR_Ascan_BoundaryAligned_Legend.png', target_bscan_index));
outputs.phr_overlay_compare = fullfile(output_dir, sprintf('BScan%04d_PhR_Ascan_Overlay_Compare.png', target_bscan_index));
outputs.phr_colorbar_horizontal = fullfile(output_dir, sprintf('BScan%04d_PhR_Colorbar_H.png', target_bscan_index));
outputs.phr_overlay_m1 = fullfile(output_dir, sprintf('BScan%04d_PhR_Ascan_Overlay_%s.png', target_bscan_index, safe_name(method1_name)));
outputs.phr_overlay_m2 = fullfile(output_dir, sprintf('BScan%04d_PhR_Ascan_Overlay_%s.png', target_bscan_index, safe_name(method2_name)));

plot_phr_ascan_boundary_aligned(ascan_table, outputs.ascan_plot, outputs.ascan_plot_legend, target_bscan_index, ascan_x, method1_name, method2_name, ascan_axial_resolution_um_per_px);
plot_phr_ascan_overlay_compare(PhR1_unit, PhR2_unit, ascan_x, bundle1.topLine_bscan, bundle2.topLine_bscan, phr_display_range, method1_name, method2_name, outputs.phr_overlay_compare, outputs.phr_overlay_m1, outputs.phr_overlay_m2);
export_colorbar_only(outputs.phr_colorbar_horizontal, parula(256), phr_colorbar_range, phr_colorbar_ticks, '°', 'south');

%% ==================== 10. Console summary ====================
print_summary(bundle1, bundle2, target_bscan_index, oa_stats_m1, oa_stats_m2, oa_neighbor_m1, oa_neighbor_m2, ascan_stats, ascan_x, outputs);
fprintf('[%s] Done. Total elapsed %.1f s\n', datestr(now, 'HH:MM:SS'), toc(run_tic));

%% ========================== Local Functions ==========================
function [dir1, dir2] = resolve_method_dirs(root_dir, dir1, dir2)
    if ~isempty(dir1) && ~exist(dir1, 'dir')
        error('method1_dir does not exist: %s', dir1);
    end
    if ~isempty(dir2) && ~exist(dir2, 'dir')
        error('method2_dir does not exist: %s', dir2);
    end

    if ~isempty(dir1) && ~isempty(dir2)
        return;
    end

    if isempty(root_dir) || ~exist(root_dir, 'dir')
        error('input_root_dir is empty or does not exist.');
    end

    discovered = discover_bscan_dirs(root_dir);
    if isempty(discovered)
        error('No B-Scan stats directories found under root: %s', root_dir);
    end

    if isempty(dir1)
        dir1 = discovered(1).data_dir;
        if strcmpi(dir1, dir2)
            error('Method 1 and Method 2 dirs resolved to the same path.');
        end
        fprintf('Auto Method 1 dir: %s\n', dir1);
    end

    if isempty(dir2)
        if numel(discovered) < 2
            error('Only one B-Scan stats directory found under root. Please set method2_dir manually.');
        end
        dir2 = discovered(2).data_dir;
        fprintf('Auto Method 2 dir: %s\n', dir2);
    end
end

function bundle = load_bscan_bundle(method_dir, target_bscan_index, oa_source, phr_scale)
    if nargin < 4 || isempty(phr_scale)
        phr_scale = 180/pi;
    end

    bundle = struct();
    bundle.method_dir = method_dir;
    bundle.target_bscan_index = target_bscan_index;

    bundle.metadata_file = local_find_exact_file(method_dir, '*_Metadata_bscan.mat');
    bundle.oa_file = local_find_exact_file(method_dir, '*_OA_bscan.mat');
    bundle.phr_file = local_find_exact_file(method_dir, '*_PhR_bscan.mat');

    if isempty(bundle.metadata_file)
        error('Missing *_Metadata_bscan.mat in %s', method_dir);
    end
    if isempty(bundle.oa_file)
        error('Missing *_OA_bscan.mat in %s', method_dir);
    end
    if isempty(bundle.phr_file)
        error('Missing *_PhR_bscan.mat in %s', method_dir);
    end

    local_tic = tic;
    fprintf('  [load] Reading %s for B-Scan %d...\n', method_dir, target_bscan_index);

    meta = local_load_mat_bundle(bundle.metadata_file);
    oa_mat = local_load_mat_bundle(bundle.oa_file);
    phr_mat = local_load_mat_bundle(bundle.phr_file);

    if ~isfield(meta, 'dopu_splitSpectrum')
        error('dopu_splitSpectrum not found in %s', bundle.metadata_file);
    end
    if ~isfield(meta, 'topLines')
        error('topLines not found in %s', bundle.metadata_file);
    end
    if ~isfield(oa_mat, 'cumLA_cfg1_avg')
        error('cumLA_cfg1_avg not found in %s', bundle.oa_file);
    end
    if ~isfield(phr_mat, 'PhR_c_cfg1_avg')
        error('PhR_c_cfg1_avg not found in %s', bundle.phr_file);
    end

    bundle.metadata = meta;
    bundle.topLines = meta.topLines;

    oa_cube = oa_mat.cumLA_cfg1_avg;
    phr_cube = phr_mat.PhR_c_cfg1_avg;

    bundle.bscan_count = size(oa_cube, 4);
    if target_bscan_index < 1 || target_bscan_index > bundle.bscan_count
        error('target_bscan_index=%d is out of range. OA cube has %d B-Scans.', target_bscan_index, bundle.bscan_count);
    end

    bundle.topLine_bscan = local_extract_topline(meta.topLines, target_bscan_index);
    [bundle.OA_deg, bundle.OA_mag, bundle.OA_vec] = extract_oa_slice(oa_cube, target_bscan_index, oa_source);
    bundle.OA_rgb = make_oa_display_rgb(bundle.OA_vec, bundle.topLine_bscan);

    [bundle.PhR_rad, bundle.PhR_unit] = extract_phr_slice(phr_cube, target_bscan_index, phr_scale);
    bundle.PRRrg = get_prrrg_from_meta(meta);
    bundle.PhR_rgb = make_phr_display_rgb(bundle.PhR_rad, bundle.PRRrg, bundle.topLine_bscan);

    bundle.base_img = bundle.OA_deg;
    bundle.oa_source = oa_source;
    bundle.phr_scale = phr_scale;
    bundle.OA_source_used = bundle.oa_source;
    bundle.PhR_source_used = 'PhR_c_cfg1_avg';
    fprintf('  [load] Finished %s in %.1f s (OA [%d,%d], PhR [%d,%d])\n', method_dir, toc(local_tic), size(bundle.OA_deg,1), size(bundle.OA_deg,2), size(bundle.PhR_unit,1), size(bundle.PhR_unit,2));
end

function S = local_load_mat_bundle(file_path)
    try
        m = matfile(file_path);
        vars = who(m);
        S = struct();
        for i = 1:numel(vars)
            name = vars{i};
            try
                S.(name) = m.(name);
            catch
                S = load(file_path);
                return;
            end
        end
    catch
        S = load(file_path);
    end
end

function file_path = local_find_exact_file(dir_path, pattern)
    d = dir(fullfile(dir_path, pattern));
    if isempty(d)
        file_path = '';
    else
        file_path = fullfile(d(1).folder, d(1).name);
    end
end

function topLine = local_extract_topline(topLines, bscan_index)
    topLines = double(topLines);
    if isvector(topLines)
        topLine = topLines(:);
        return;
    end
    if size(topLines, 2) >= bscan_index
        topLine = topLines(:, bscan_index);
    else
        topLine = topLines(:, end);
    end
end

function [OA_deg, OA_mag, OA_vec] = extract_oa_slice(oa_cube, bscan_index, oa_source)
    oa_cube = double(oa_cube);

    if ndims(oa_cube) ~= 4 || size(oa_cube, 3) < 3
        error('OA cube must be [Z,X,3,B]. Got size %s', mat2str(size(oa_cube)));
    end

    if strcmpi(oa_source, 'local')
        warning('oa_source=local requested, but only cumLA_cfg1_avg is available in this export. Using cum data.');
    end

    OA_vec = oa_cube(:,:,:,bscan_index);
    OA_mag = sqrt(sum(OA_vec.^2, 3));
    OA_deg = mod(0.5 * atan2d(OA_vec(:,:,2), OA_vec(:,:,1)), 180);
end

function [PhR_rad, PhR_unit] = extract_phr_slice(phr_cube, bscan_index, phr_scale)
    phr_cube = double(phr_cube);
    if ndims(phr_cube) ~= 3
        error('PhR cube must be [Z,X,B]. Got size %s', mat2str(size(phr_cube)));
    end
    PhR_rad = phr_cube(:,:,bscan_index);
    PhR_unit = PhR_rad * phr_scale;
end

function discovered = discover_bscan_dirs(root_dir)
    files = [dir(fullfile(root_dir, '**', '*_bscan*.mat')); dir(fullfile(root_dir, '**', '*_BScan*.mat')); dir(fullfile(root_dir, '**', '*_BSCAN*.mat'))];
    if isempty(files)
        discovered = struct('data_dir', {}, 'method_name', {}, 'session_name', {});
        return;
    end

    folders = unique(string({files.folder}));
    records = struct('data_dir', {}, 'method_name', {}, 'session_name', {});

    for i = 1:numel(folders)
        f = char(folders(i));
        if ~local_has_bscan_files(f)
            continue;
        end
        [method_name, session_name] = local_extract_names(root_dir, f);
        rec.data_dir = f;
        rec.method_name = method_name;
        rec.session_name = session_name;
        records(end+1) = rec; %#ok<AGROW>
    end

    if isempty(records)
        discovered = struct('data_dir', {}, 'method_name', {}, 'session_name', {});
        return;
    end

    [~, order] = sort(lower(string({records.method_name})));
    discovered = records(order);
end

function tf = local_has_bscan_files(dir_path)
    pats = {'**/*_bscan*.mat', '**/*_BScan*.mat', '**/*_BSCAN*.mat'};
    tf = false;
    for i = 1:numel(pats)
        if ~isempty(dir(fullfile(dir_path, pats{i})))
            tf = true;
            return;
        end
    end
end

function [method_name, session_name] = local_extract_names(root_dir, data_dir)
    root_clean = char(string(root_dir));
    data_clean = char(string(data_dir));

    if endsWith(root_clean, filesep)
        root_clean = root_clean(1:end-1);
    end

    rel = strrep(data_clean, [root_clean filesep], '');
    parts = strsplit(rel, filesep);

    session_name = '';
    method_name = 'Method';
    if isempty(parts)
        return;
    end

    session_name = parts{end};
    if numel(parts) >= 2
        method_name = parts{end-1};
    else
        [~, method_name] = fileparts(data_clean);
    end
end

function roi_mask = load_or_draw_roi_mask(roi_path, roi_mode, base_img, output_dir, bscan_index, roi_annotation_image, data_size)
    if isempty(roi_path)
        roi_path = fullfile(output_dir, sprintf('ROI_mask_BScan_%04d.mat', bscan_index));
    end

    can_load = exist(roi_path, 'file') == 2;

    if strcmpi(roi_mode, 'load')
        if ~can_load
            error('roi_mask_path not found: %s', roi_path);
        end
        tmp = load(roi_path);
        if isfield(tmp, 'roi_mask')
            roi_mask = tmp.roi_mask;
            fprintf('  ROI mask loaded from %s\n', roi_path);
            return;
        else
            error('roi_mask not found in %s', roi_path);
        end
    end

    if strcmpi(roi_mode, 'from_annotation_image')
        if isempty(roi_annotation_image) || ~exist(roi_annotation_image, 'file')
            error('roi_annotation_image not found: %s', roi_annotation_image);
        end
        roi_mask = extract_roi_mask_from_annotation(roi_annotation_image, data_size);
        fprintf('  ROI mask extracted from annotation: %s\n', roi_annotation_image);
        return;
    end

    if strcmpi(roi_mode, 'auto_default')
        roi_mask = auto_center_rect_roi(size(base_img));
        return;
    end

    if strcmpi(roi_mode, 'load_or_draw') && can_load
        tmp = load(roi_path);
        if isfield(tmp, 'roi_mask')
            roi_mask = tmp.roi_mask;
            fprintf('Loaded ROI mask: %s\n', roi_path);
            return;
        end
    end

    if strcmpi(roi_mode, 'draw') || strcmpi(roi_mode, 'load_or_draw')
        fig = figure('Color', 'w', 'Name', 'Draw ROI');
        imagesc(base_img); axis image; colormap(gca, gray); colorbar;
        title('Draw ROI on B-Scan (size [Z,X])', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        xlabel('X (column)'); ylabel('Z (row)');

        try
            roi_mask = draw_roi_mask(size(base_img));
            if isempty(roi_mask)
                roi_mask = auto_center_rect_roi(size(base_img));
                warning('ROI drawing returned empty mask. Falling back to automatic centered rectangle ROI.');
            end
            fprintf('ROI mask captured from drawing.\n');
        catch ME
            fprintf('  ROI drawing failed: %s. Falling back to automatic centered rectangle ROI.\n', ME.message);
            roi_mask = auto_center_rect_roi(size(base_img));
        end
        close(fig);
        return;
    end

    error('Invalid roi_mode: %s', roi_mode);
end

function roi_mask = draw_roi_mask(img_size)
    try
        h = drawpolygon('LineWidth', 1.5);
        wait(h);
        roi_mask = createMask(h);
    catch
        try
            roi_mask = roipoly;
        catch
            roi_mask = [];
        end
    end

    if ~isempty(roi_mask) && ~isequal(size(roi_mask), img_size)
        roi_mask = roi_mask(1:img_size(1), 1:img_size(2));
    end
end

function roi_mask = auto_center_rect_roi(img_size)
    z_size = img_size(1);
    x_size = img_size(2);

    rect_h = max(20, round(z_size * 0.4));
    rect_w = max(40, round(x_size * 0.4));

    z1 = max(1, floor((z_size - rect_h) / 2) + 1);
    x1 = max(1, floor((x_size - rect_w) / 2) + 1);
    z2 = min(z_size, z1 + rect_h - 1);
    x2 = min(x_size, x1 + rect_w - 1);

    roi_mask = false(z_size, x_size);
    roi_mask(z1:z2, x1:x2) = true;
    fprintf('  Using automatic centered ROI rectangle: Z[%d:%d], X[%d:%d]\n', z1, z2, x1, x2);
end

function [oa_stats, local_data] = compute_oa_local_consistency(OA_deg, roi_mask, window_radius_px)
    valid = roi_mask & isfinite(OA_deg);
    theta = deg2rad(2 * OA_deg);
    C = cos(theta);
    S = sin(theta);

    kernel = ones(2 * window_radius_px + 1);
    weight = conv2(double(valid), kernel, 'same');
    mean_C = conv2(C .* double(valid), kernel, 'same') ./ max(weight, eps);
    mean_S = conv2(S .* double(valid), kernel, 'same') ./ max(weight, eps);

    local_R = sqrt(mean_C.^2 + mean_S.^2);
    local_R(weight == 0) = NaN;
    local_R = min(max(local_R, eps), 1);

    local_std_deg = rad2deg(0.5 * sqrt(-2 * log(local_R)));
    local_std_deg(~valid) = NaN;

    local_std_vals = local_std_deg(valid);
    local_R_vals = local_R(valid);

    oa_stats = struct();
    oa_stats.valid_count = nnz(valid);
    if isempty(local_std_vals)
        oa_stats.median_local_variability = NaN;
        oa_stats.p90_local_variability = NaN;
    else
        oa_stats.median_local_variability = median(local_std_vals, 'omitnan');
        oa_stats.p90_local_variability = prctile(local_std_vals, 90);
    end
    if isempty(local_R_vals)
        oa_stats.mean_local_coherence_R = NaN;
    else
        oa_stats.mean_local_coherence_R = mean(local_R_vals, 'omitnan');
    end

    local_data = struct();
    local_data.local_std_deg = local_std_deg;
    local_data.local_R = local_R;
    local_data.local_std_vals = local_std_vals;
    local_data.local_R_vals = local_R_vals;
end

function stats = compute_neighbor_angle_change(OA_deg, roi_mask)
    valid = roi_mask & isfinite(OA_deg);

    dx = circular_diff_deg(OA_deg(:, 2:end), OA_deg(:, 1:end-1));
    dz = circular_diff_deg(OA_deg(2:end, :), OA_deg(1:end-1, :));

    dx_mask = valid(:, 2:end) & valid(:, 1:end-1);
    dz_mask = valid(2:end, :) & valid(1:end-1, :);

    vals = [abs(dx(dx_mask)); abs(dz(dz_mask))];

    stats = struct();
    if isempty(vals)
        stats.mean_neighbor_change = NaN;
        stats.median_neighbor_change = NaN;
        stats.p90_neighbor_change = NaN;
    else
        stats.mean_neighbor_change = mean(vals, 'omitnan');
        stats.median_neighbor_change = median(vals, 'omitnan');
        stats.p90_neighbor_change = prctile(vals, 90);
    end
end

function plot_oa_roi_overlay_compare(OA1_rgb, OA2_rgb, roi_mask, name1, name2, out_path, out_path_m1, out_path_m2)
    style = get_plot_style();
    fig = figure('Color', 'w', 'Position', [100 100 1200 430]);
    tiledlayout(1, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

    nexttile;
    imshow(OA1_rgb);
    hold on;
    draw_roi_boundary(roi_mask, style.boundary_linewidth);
    hold off;
    axis image off;
    title(name1, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);

    nexttile;
    imshow(OA2_rgb);
    hold on;
    draw_roi_boundary(roi_mask, style.boundary_linewidth);
    hold off;
    axis image off;
    title(name2, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);

    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);

    if nargin >= 7 && ~isempty(out_path_m1)
        export_single_roi_image(OA1_rgb, roi_mask, name1, out_path_m1, style);
    end
    if nargin >= 8 && ~isempty(out_path_m2)
        export_single_roi_image(OA2_rgb, roi_mask, name2, out_path_m2, style);
    end
end

function plot_oa_roi_crop_compare(OA1_rgb, OA2_rgb, roi_mask, name1, name2, out_path, out_path_m1, out_path_m2)
    style = get_plot_style();
    margin_px = 30;
    [r0, r1, c0, c1] = roi_crop_bounds(roi_mask, margin_px);

    img1 = OA1_rgb(r0:r1, c0:c1, :);
    img2 = OA2_rgb(r0:r1, c0:c1, :);
    roi_crop = roi_mask(r0:r1, c0:c1);

    fig = figure('Color', 'w', 'Position', [100 100 1200 430]);
    tiledlayout(1, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

    nexttile;
    imshow(img1);
    hold on;
    draw_roi_boundary(roi_crop, style.boundary_linewidth);
    hold off;
    axis image off;
    title(name1, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);

    nexttile;
    imshow(img2);
    hold on;
    draw_roi_boundary(roi_crop, style.boundary_linewidth);
    hold off;
    axis image off;
    title(name2, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);

    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);

    if nargin >= 7 && ~isempty(out_path_m1)
        export_single_roi_image(img1, roi_crop, name1, out_path_m1, style);
    end
    if nargin >= 8 && ~isempty(out_path_m2)
        export_single_roi_image(img2, roi_crop, name2, out_path_m2, style);
    end
end

function plot_oa_local_consistency_figure(local_std1, local_std2, oa_stats_m1, oa_stats_m2, neigh_m1, neigh_m2, name1, name2, out_path, legend_density_path, legend_metrics_path)
    colors = get_plot_colors();
    style = get_plot_style();
    fig = figure('Color', 'w', 'Position', [100 100 1320 600]);

    deg = char(176);
    legend_colors = [colors.ddg; colors.bayes];

    % Use manual axes positions so large labels have fixed reserved space.
    ax1 = axes('Parent', fig, 'Position', [0.12 0.22 0.38 0.66]);
    hold(ax1, 'on');
    if ~isempty(local_std1)
        histogram(ax1, local_std1, 'Normalization', 'pdf', 'FaceColor', colors.ddg, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        [x1, y1] = compute_kde(local_std1);
        plot(ax1, x1, y1, '-', 'Color', colors.ddg, 'LineWidth', style.line_width);
    end
    if ~isempty(local_std2)
        histogram(ax1, local_std2, 'Normalization', 'pdf', 'FaceColor', colors.bayes, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        [x2, y2] = compute_kde(local_std2);
        plot(ax1, x2, y2, '--', 'Color', colors.bayes, 'LineWidth', style.line_width);
    end
    hold(ax1, 'off');

    xlabel(ax1, sprintf('Variability (%s)', deg), 'FontName', 'Times New Roman', 'FontSize', style.font_label, 'FontWeight', 'bold');
    ylabel(ax1, 'Probability density', 'FontName', 'Times New Roman', 'FontSize', style.font_label, 'FontWeight', 'bold');
    title(ax1, 'Local OA variability', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    grid(ax1, 'on');
    set(ax1, 'GridColor', colors.grid, 'FontName', 'Times New Roman', 'FontSize', style.font_tick, ...
        'LineWidth', style.axis_linewidth, 'TickDir', 'in', 'FontWeight', 'bold');

    export_legend_only({'-', '--'}, {'none', 'none'}, legend_colors, {name1, name2}, legend_density_path, 'northeast', style.font_legend);

    ax2 = axes('Parent', fig, 'Position', [0.62 0.22 0.34 0.66]);
    metrics = [oa_stats_m1.median_local_variability, oa_stats_m1.p90_local_variability, neigh_m1.mean_neighbor_change; ...
               oa_stats_m2.median_local_variability, oa_stats_m2.p90_local_variability, neigh_m2.mean_neighbor_change]';

    b = bar(ax2, metrics, 'grouped');
    b(1).FaceColor = colors.ddg;
    b(2).FaceColor = colors.bayes;
    b(1).EdgeColor = 'none';
    b(2).EdgeColor = 'none';

    set(ax2, 'XTick', 1:3, ...
        'XTickLabel', {'', '', ''}, ...
        'XTickLabelRotation', 0, ...
        'XLim', [0.5 3.5]);

    ylabel(ax2, sprintf('Degrees (%s)', deg), 'FontName', 'Times New Roman', 'FontSize', style.font_label, 'FontWeight', 'bold');
    title(ax2, 'Core local consistency metrics', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);

    grid(ax2, 'on');
    set(ax2, 'GridColor', colors.grid, 'FontName', 'Times New Roman', 'FontSize', style.font_tick, ...
        'LineWidth', style.axis_linewidth, 'TickDir', 'in', 'FontWeight', 'bold');

    % Freeze Y limits before drawing manual labels, otherwise MATLAB may expand
    % the axis to include the text below the x-axis.
    yl = ylim(ax2);
    ylim(ax2, yl);
    yr = yl(2) - yl(1);

    label_y1 = yl(1) - 0.045 * yr;
    label_y2 = yl(1) - 0.115 * yr;

    label_top = {'Median', 'P90', 'Mean'};
    label_bottom = {'local \sigma', 'local \sigma', 'neighbor \Delta'};

    text(ax2, 1:3, repmat(label_y1, 1, 3), label_top, ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.font_tick, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Clipping', 'off');

    text(ax2, 1:3, repmat(label_y2, 1, 3), label_bottom, ...
        'FontName', 'Times New Roman', ...
        'FontSize', style.font_tick, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'Clipping', 'off');

    export_legend_only({'none', 'none'}, {'s', 's'}, legend_colors, {name1, name2}, legend_metrics_path, 'northwest', style.font_legend);

    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);
end

function [xs, ys] = compute_kde(vals)
    vals = vals(isfinite(vals));
    if numel(vals) < 3
        xs = vals(:)';
        ys = zeros(size(xs));
        return;
    end

    if exist('ksdensity', 'file') == 2
        [ys, xs] = ksdensity(vals);
        return;
    end

    xmin = min(vals);
    xmax = max(vals);
    if xmin == xmax
        xs = xmin;
        ys = 1;
        return;
    end

    xs = linspace(xmin, xmax, 200);
    bw = 1.06 * std(vals) * numel(vals)^(-1/5);
    if ~isfinite(bw) || bw <= 0
        bw = (xmax - xmin) / 30;
    end

    diffs = (xs(:) - vals(:)')./ bw;
    ys = mean(exp(-0.5 * diffs.^2), 2) ./ (bw * sqrt(2 * pi));
    ys = ys(:)';
end

function ascan_x = resolve_ascan_x(ascan_x, base_img, output_dir, bscan_index, ascan_mode, ascan_annotation_image, data_size)
    if isempty(ascan_x) && strcmpi(ascan_mode, 'from_annotation_image')
        if isempty(ascan_annotation_image) || ~exist(ascan_annotation_image, 'file')
            error('ascan_annotation_image not found: %s', ascan_annotation_image);
        end
        ascan_x = extract_ascan_from_annotation(ascan_annotation_image, data_size);
    elseif isempty(ascan_x)
        fig = figure('Color', 'w', 'Name', 'Select A-scan X');
        imagesc(base_img);
        colormap(gca, gray);
        title('Click to choose A-scan X on B-Scan', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        axis image; colorbar; xlabel('X (column)'); ylabel('Z (row)');
        fprintf('  Please click one X position in the figure window...\n');
        [x, ~] = ginput(1);
        close(fig);
        if isempty(x)
            error('No A-scan X selected.');
        end
        ascan_x = round(x);
    end

    fprintf('A-scan X resolved: %d\n', ascan_x);
end

function [ascan_stats, T] = compute_phr_ascan_stats_from_boundary(PhR1_unit, PhR2_unit, topLine1, topLine2, ascan_x, radius_px, name1, name2, start_mode, depth_axis)
    [Z, X] = size(PhR1_unit);
    if ascan_x < 1 || ascan_x > X
        error('ascan_x out of bounds: %d (X size=%d)', ascan_x, X);
    end

    if nargin < 9 || isempty(start_mode)
        start_mode = 'boundary';
    end
    if nargin < 10 || isempty(depth_axis)
        depth_axis = 'relative';
    end

    x0 = max(1, ascan_x - radius_px);
    x1 = min(X, ascan_x + radius_px);

    if strcmpi(start_mode, 'boundary')
        [v1, v2, z_vals] = extract_boundary_aligned_ascan(PhR1_unit, PhR2_unit, topLine1, topLine2, x0, x1);
        if strcmpi(depth_axis, 'absolute')
            z_vals = z_vals + (max(1, round(max([topLine1(ascan_x), topLine2(ascan_x)]))) - 1);
        end
    else
        if radius_px > 0
            v1 = mean(PhR1_unit(:, x0:x1), 2, 'omitnan');
            v2 = mean(PhR2_unit(:, x0:x1), 2, 'omitnan');
        else
            v1 = PhR1_unit(:, ascan_x);
            v2 = PhR2_unit(:, ascan_x);
        end
        z_vals = (1:Z).';
    end

    diff = v2 - v1;
    valid = isfinite(v1) & isfinite(v2);

    ascan_stats = struct();
    ascan_stats.method1_name = name1;
    ascan_stats.method2_name = name2;
    ascan_stats.valid_count = nnz(valid);
    ascan_stats.method1_mean = mean(v1(valid));
    ascan_stats.method1_std = std(v1(valid));
    ascan_stats.method1_median = median(v1(valid));
    ascan_stats.method2_mean = mean(v2(valid));
    ascan_stats.method2_std = std(v2(valid));
    ascan_stats.method2_median = median(v2(valid));
    ascan_stats.mae = mean(abs(diff(valid)));
    ascan_stats.rmse = sqrt(mean(diff(valid).^2));
    ascan_stats.mean_diff = mean(diff(valid));
    ascan_stats.max_abs_diff = max(abs(diff(valid)));
    if nnz(valid) >= 2
        ascan_stats.pearson_r = corr(v1(valid), v2(valid));
    else
        ascan_stats.pearson_r = NaN;
    end

    T = table(z_vals, v1, v2, diff, 'VariableNames', {'z', 'PhR_method1', 'PhR_method2', 'diff'});
end

function [v1, v2, rel_z] = extract_boundary_aligned_ascan(PhR1_unit, PhR2_unit, topLine1, topLine2, x0, x1)
    [Z, ~] = size(PhR1_unit);
    cols = x0:x1;
    n_cols = numel(cols);

    z0_1 = round(topLine1(cols));
    z0_2 = round(topLine2(cols));
    z0_1 = max(1, min(Z, z0_1));
    z0_2 = max(1, min(Z, z0_2));
    z0 = max([z0_1(:); z0_2(:); 1]);

    max_len = Z - z0 + 1;
    v1_stack = NaN(max_len, n_cols);
    v2_stack = NaN(max_len, n_cols);

    for i = 1:n_cols
        c = cols(i);
        zc = max([z0_1(i), z0_2(i), z0]);
        if zc <= Z
            v1_stack(1:(Z - zc + 1), i) = PhR1_unit(zc:Z, c);
            v2_stack(1:(Z - zc + 1), i) = PhR2_unit(zc:Z, c);
        end
    end

    v1 = mean(v1_stack, 2, 'omitnan');
    v2 = mean(v2_stack, 2, 'omitnan');
    rel_z = (0:(max_len - 1)).';
end

function plot_phr_ascan_boundary_aligned(T, out_path, legend_path, bscan_idx, ascan_x, name1, name2, axial_res_um_per_px)
    colors = get_plot_colors();
    style = get_plot_style();
    fig = figure('Color', 'w', 'Position', [100 100 900 900]);
    tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    marker_step = 15;
    if nargin < 8 || isempty(axial_res_um_per_px) || axial_res_um_per_px <= 0
        axial_res_um_per_px = 1.0;
    end
    idx = 1:marker_step:numel(T.z);
    deg = char(176);
    mu = char(181);
    unit_label = sprintf('%s/%sm', deg, mu);
    depth_um = T.z * axial_res_um_per_px;
    phr1_rate = T.PhR_method1 / axial_res_um_per_px;
    phr2_rate = T.PhR_method2 / axial_res_um_per_px;

    nexttile;
    plot(depth_um, phr1_rate, '-', 'Color', colors.ddg, 'LineWidth', style.line_width, 'Marker', 'o', ...
        'MarkerIndices', idx, 'MarkerSize', style.marker_size, 'MarkerFaceColor', colors.ddg);
    hold on;
    plot(depth_um, phr2_rate, '--', 'Color', colors.bayes, 'LineWidth', style.line_width, 'Marker', 's', ...
        'MarkerIndices', idx, 'MarkerSize', style.marker_size, 'MarkerFaceColor', colors.bayes);
    hold off;
    xlabel(sprintf('Depth from boundary (%sm)', mu), 'FontName', 'Times New Roman', 'FontSize', style.font_label, 'FontWeight', 'bold');
    ylabel(sprintf('Local Phase Retardation (%s)', unit_label), 'FontName', 'Times New Roman', 'FontSize', style.font_label, 'FontWeight', 'bold');
    title(sprintf('PhR A-scan (BScan=%d, X=%d, boundary-aligned)', bscan_idx, ascan_x), ...
        'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    legend_colors = [colors.ddg; colors.bayes];
    export_legend_only({'-', '--'}, {'o', 's'}, legend_colors, {name1, name2}, legend_path, 'northeast', style.font_legend, style.marker_size + 1);
    grid on;
    set(gca, 'GridColor', colors.grid, 'FontName', 'Times New Roman', 'FontSize', style.font_tick, 'LineWidth', style.axis_linewidth, 'TickDir', 'in', 'FontWeight', 'bold');

    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);
end

function plot_phr_ascan_overlay_compare(PhR1_unit, PhR2_unit, ascan_x, topLine1, topLine2, disp_range, name1, name2, out_path, out_path_m1, out_path_m2)
    style = get_plot_style();
    fig = figure('Color', 'w', 'Position', [100 100 1200 430]);
    tiledlayout(1, 2, 'Padding', 'loose', 'TileSpacing', 'loose');

    nexttile;
    rgb1 = make_phr_display_rgb_unit(PhR1_unit, disp_range, topLine1);
    imshow(rgb1);
    hold on;
    xline(ascan_x, 'w-', 'LineWidth', style.line_width);
    hold off;
    axis image off;
    title(name1, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);

    nexttile;
    rgb2 = make_phr_display_rgb_unit(PhR2_unit, disp_range, topLine2);
    imshow(rgb2);
    hold on;
    xline(ascan_x, 'w-', 'LineWidth', style.line_width);
    hold off;
    axis image off;
    title(name2, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);

    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);

    if nargin >= 10 && ~isempty(out_path_m1)
        export_single_phr_overlay(PhR1_unit, disp_range, topLine1, ascan_x, name1, out_path_m1, style);
    end
    if nargin >= 11 && ~isempty(out_path_m2)
        export_single_phr_overlay(PhR2_unit, disp_range, topLine2, ascan_x, name2, out_path_m2, style);
    end
end

function colors = get_plot_colors()
    colors.ddg = [0.000, 0.447, 0.741];
    colors.bayes = [0.835, 0.369, 0.000];
    colors.diff = [0.13, 0.13, 0.13];
    colors.ref = [0.70, 0.70, 0.70];
    colors.grid = [0.88, 0.88, 0.88];
end

function style = get_plot_style()
    font_scale = 1.2;
    style.font_title = round(20 * font_scale);
    style.font_label = round(18 * font_scale);
    style.font_tick = round(16 * font_scale);
    style.font_legend = round(16 * font_scale);
    style.axis_linewidth = 1.6;
    style.line_width = 2.2;
    style.boundary_linewidth = 2.8;
    style.marker_size = 5;
    style.top_margin = 0.18;
end

function export_legend_only(line_styles, markers, colors, labels, out_path, location, font_size, marker_size)
    if nargin < 7 || isempty(font_size)
        font_size = 10;
    end
    if nargin < 8 || isempty(marker_size)
        marker_size = 8;
    end

    n_items = numel(labels);
    fig = figure('Color', 'w', 'Position', [100 100 360 220]);
    ax = axes(fig);
    hold(ax, 'on');
    h = gobjects(n_items, 1);

    for i = 1:n_items
        h(i) = plot(ax, nan, nan, 'LineStyle', line_styles{i}, 'Marker', markers{i}, 'Color', colors(i,:), 'LineWidth', 2.0);
        if ~strcmpi(markers{i}, 'none')
            set(h(i), 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:), 'MarkerSize', marker_size);
        end
    end

    hold(ax, 'off');
    axis(ax, 'off');
    set(ax, 'Position', [0 0 1 1]);
    legend(ax, h, labels, 'Location', location, 'Box', 'off', 'FontSize', font_size, 'FontName', 'Times New Roman');
    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);
end

function [r0, r1, c0, c1] = roi_crop_bounds(roi_mask, margin_px)
    [rows, cols] = find(roi_mask);
    if isempty(rows)
        r0 = 1; r1 = size(roi_mask, 1);
        c0 = 1; c1 = size(roi_mask, 2);
        return;
    end
    r0 = max(1, min(rows) - margin_px);
    r1 = min(size(roi_mask, 1), max(rows) + margin_px);
    c0 = max(1, min(cols) - margin_px);
    c1 = min(size(roi_mask, 2), max(cols) + margin_px);
end

function draw_roi_boundary(roi_mask, line_width)
    if nargin < 2 || isempty(line_width)
        line_width = 2.0;
    end
    bnd = bwboundaries(roi_mask, 'noholes');
    for k = 1:numel(bnd)
        plot(bnd{k}(:,2), bnd{k}(:,1), 'w', 'LineWidth', line_width);
    end
end

function export_single_roi_image(img, roi_mask, title_str, out_path, style)
    fig = figure('Color', 'w', 'Position', [100 100 600 430]);
    imshow(img);
    hold on;
    if ~isempty(roi_mask)
        draw_roi_boundary(roi_mask, style.boundary_linewidth);
    end
    hold off;
    axis image off;
    title(title_str, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);
    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);
end

function export_single_phr_overlay(PhR_unit, disp_range, topLine, ascan_x, title_str, out_path, style)
    fig = figure('Color', 'w', 'Position', [100 100 600 430]);
    rgb = make_phr_display_rgb_unit(PhR_unit, disp_range, topLine);
    imshow(rgb);
    hold on;
    xline(ascan_x, 'w-', 'LineWidth', style.line_width);
    hold off;
    axis image off;
    title(title_str, 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', style.font_title);
    apply_top_margin(gca, style.top_margin);
    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);
end

function apply_top_margin(ax, top_margin)
    if nargin < 2 || isempty(top_margin)
        top_margin = 0.15;
    end
    if isa(ax.Parent, 'matlab.graphics.layout.TiledChartLayout')
        if ~isempty(ax.Title) && isgraphics(ax.Title)
            ax.Title.Units = 'normalized';
            tpos = ax.Title.Position;
            tpos(2) = min(tpos(2), 1 - top_margin);
            ax.Title.Position = tpos;
        end
        return;
    end

    ax.Units = 'normalized';
    pos = ax.Position;
    desired_top = 1 - top_margin;
    if desired_top > pos(2)
        pos(4) = max(0.05, desired_top - pos(2));
        ax.Position = pos;
    end
    if ~isempty(ax.Title) && isgraphics(ax.Title)
        ax.Title.Units = 'normalized';
        tpos = ax.Title.Position;
        tpos(2) = min(tpos(2), 1 - top_margin);
        ax.Title.Position = tpos;
    end
end

function export_colorbar_only(out_path, cmap, clim, ticks, label_str, location)
    style = get_plot_style();
    if nargin < 6 || isempty(location)
        location = 'eastoutside';
    end
    if strcmpi(location, 'south')
        fig = figure('Color', 'w', 'Position', [100 100 900 120]);
    else
        fig = figure('Color', 'w', 'Position', [100 100 220 420]);
    end
    ax = axes(fig);
    h_img = imagesc(ax, linspace(clim(1), clim(2), 256).');
    set(h_img, 'AlphaData', 0); % hide background image, keep colorbar only
    set(ax, 'Visible', 'off', 'Color', 'none');
    colormap(ax, cmap);
    caxis(ax, clim);
    cb = colorbar(ax, 'Location', location);
    cb.Ticks = ticks;
    cb.FontSize = style.font_tick;
    cb.LineWidth = style.axis_linewidth;
    cb.Label.String = label_str;
    cb.Label.FontSize = style.font_label;
    cb.Label.FontWeight = 'bold';
    exportgraphics(fig, out_path, 'Resolution', 600);
    close(fig);
end

function disp_range = resolve_phr_display_range(v1, v2, disp_range)
    if nargin < 3 || isempty(disp_range)
        vals = [v1(:); v2(:)];
        vals = vals(isfinite(vals));
        if isempty(vals)
            disp_range = [0 1];
        else
            disp_range = [min(vals), max(vals)];
            if disp_range(1) == disp_range(2)
                disp_range(2) = disp_range(1) + eps;
            end
        end
    end
end

function diff_deg = circular_diff_deg(ang2_deg, ang1_deg)
    diff_deg = mod(ang2_deg - ang1_deg + 90, 180) - 90;
end

function mean_deg = circular_mean_pair_deg(ang1_deg, ang2_deg)
    a = deg2rad(ang1_deg(:) * 2);
    b = deg2rad(ang2_deg(:) * 2);
    mean_x = cos(a) + cos(b);
    mean_y = sin(a) + sin(b);
    mean_deg = mod(0.5 * atan2d(mean_y, mean_x), 180);
    mean_deg = reshape(mean_deg, size(ang1_deg));
end

function ascan_x = extract_ascan_from_annotation(img_path, data_size)
    [img, map] = imread(img_path);
    if ~isempty(map)
        img = ind2rgb(img, map);
    end
    white = extract_white_mask(img, false);

    [h, w] = size(white);
    col_counts = sum(white, 1);
    if max(col_counts) == 0
        error('No white pixels detected in %s', img_path);
    end

    cc = bwconncomp(white);
    stats = regionprops(cc, 'BoundingBox', 'Area', 'PixelIdxList');
    best_idx = 0;
    best_score = -inf;
    for i = 1:numel(stats)
        bb = stats(i).BoundingBox;
        height = bb(4);
        width = bb(3);
        if width <= 0
            continue;
        end
        aspect = height / width;
        score = stats(i).Area * aspect;
        if aspect >= 3 && height >= 0.3 * h && score > best_score
            best_score = score;
            best_idx = i;
        end
    end

    if best_idx > 0
        [~, cols] = ind2sub([h, w], stats(best_idx).PixelIdxList);
        col_center = median(cols);
    else
        [~, col_center] = max(col_counts);
    end

    ascan_x_img = round(col_center);
    ascan_x = map_coord(ascan_x_img, w, data_size(2));
    fprintf('  A-scan X extracted from annotation: img_x=%d -> data_x=%d\n', ascan_x_img, ascan_x);
end

function roi_mask = extract_roi_mask_from_annotation(img_path, data_size)
    [img, map] = imread(img_path);
    if ~isempty(map)
        img = ind2rgb(img, map);
    end
    white = extract_white_mask(img, false);
    white = bwareaopen(white, 20);

    if nnz(white) == 0
        error('No white pixels detected in %s', img_path);
    end

    cc = bwconncomp(white);
    stats = regionprops(cc, 'Area', 'PixelIdxList');
    [~, idx] = max([stats.Area]);
    mask_main = false(size(white));
    mask_main(stats(idx).PixelIdxList) = true;

    [rows, cols] = find(mask_main);
    k = convhull(cols, rows);
    roi_mask_img = poly2mask(cols(k), rows(k), size(white,1), size(white,2));
    roi_mask_img = imfill(roi_mask_img, 'holes');

    roi_mask = imresize(roi_mask_img, data_size, 'nearest');
end

function white = extract_white_mask(img, clear_border)
    if size(img, 3) == 1
        img = repmat(img, 1, 1, 3);
    end
    img = double(img);

    mx = max(img, [], 3);
    mn = min(img, [], 3);
    white = mx > 200 & mn > 170 & (mx - mn) < 70;
    white = bwareaopen(white, 5);

    if nargin >= 2 && clear_border
        white = imclearborder(white);
    end
end

function out = map_coord(coord_img, img_size, data_size)
    out = round(coord_img * data_size / img_size);
    out = max(1, min(data_size, out));
end

function valid_mask = build_valid_mask(bundle1, bundle2, OA1_deg, OA2_deg, oa_mag_threshold, use_tissue_mask, bscan_index)
    valid_mask = isfinite(OA1_deg) & isfinite(OA2_deg);
    if isfield(bundle1, 'OA_mag') && isfield(bundle2, 'OA_mag')
        valid_mask = valid_mask & (bundle1.OA_mag > oa_mag_threshold) & (bundle2.OA_mag > oa_mag_threshold);
    end

    if use_tissue_mask
        tissue_mask = derive_tissue_mask(bundle1.metadata, bundle2.metadata, size(OA1_deg), bscan_index);
        valid_mask = valid_mask & tissue_mask;
    end
end

function tissue_mask = derive_tissue_mask(meta1, meta2, data_size, bscan_index)
    tissue_mask = true(data_size);

    candidates = {'Strus', 'strus', 'dopu_splitSpectrum', 'DOPU_splitSpectrum', 'I', 'intensity'};
    mask1 = extract_mask_from_meta(meta1, candidates, data_size, bscan_index);
    mask2 = extract_mask_from_meta(meta2, candidates, data_size, bscan_index);

    if ~isempty(mask1)
        tissue_mask = tissue_mask & mask1;
    end
    if ~isempty(mask2)
        tissue_mask = tissue_mask & mask2;
    end
end

function mask = extract_mask_from_meta(meta, fields, data_size, bscan_index)
    mask = [];
    for i = 1:numel(fields)
        f = fields{i};
        if isfield(meta, f)
            raw = meta.(f);
            if ndims(raw) == 3
                if size(raw, 3) >= bscan_index
                    raw = raw(:,:,bscan_index);
                else
                    raw = raw(:,:,1);
                end
            end
            if ~isequal(size(raw), data_size)
                raw = imresize(raw, data_size, 'nearest');
            end
            raw = double(raw);
            thresh = graythresh(mat2gray(raw));
            mask = mat2gray(raw) > thresh;
            return;
        end
    end
end

function rgb = make_oa_display_rgb(OA_vec, topLine)
    rgb = quColoring(OA_vec, 440);
    rgb = apply_topline_mask_rgb(rgb, topLine);
end

function rgb = apply_topline_mask_rgb(rgb, topLine)
    if isempty(topLine)
        return;
    end
    topLine = round(double(topLine(:)));
    [Z, X, ~] = size(rgb);
    if numel(topLine) ~= X
        return;
    end
    for x = 1:X
        z1 = max(1, min(Z, topLine(x)));
        rgb(1:z1, x, :) = 0;
    end
end

function PRRrg = get_prrrg_from_meta(meta)
    PRRrg = [0.06 0.4];
    if isfield(meta, 'PRRrg')
        v = meta.PRRrg;
        if numel(v) == 2
            PRRrg = double(v(:)).';
            return;
        end
    end
    if isfield(meta, 'config_params') && isstruct(meta.config_params)
        if isfield(meta.config_params, 'polarization') && isfield(meta.config_params.polarization, 'PRRrg')
            v = meta.config_params.polarization.PRRrg;
            if numel(v) == 2
                PRRrg = double(v(:)).';
            end
        end
    end
end

function rgb = make_phr_display_rgb(PhR_rad, PRRrg, topLine)
    if nargin < 2 || isempty(PRRrg)
        PRRrg = [0.06 0.4];
    end
    scaled = mat2gray(PhR_rad, PRRrg);
    rgb = ind2rgb(uint8(scaled * 255), parula(256));

    if nargin >= 3 && ~isempty(topLine)
        rgb = apply_topline_mask_rgb(rgb, topLine);
    end
end

function rgb = make_phr_display_rgb_unit(PhR_unit, disp_range, topLine)
    if nargin < 2 || isempty(disp_range)
        disp_range = [0 5];
    end
    scaled = mat2gray(PhR_unit, disp_range);
    rgb = ind2rgb(uint8(scaled * 255), parula(256));

    if nargin >= 3 && ~isempty(topLine)
        rgb = apply_topline_mask_rgb(rgb, topLine);
    end
end


function name_out = safe_name(name_in)
    name_out = regexprep(name_in, '[^A-Za-z0-9]+', '');
    if isempty(name_out)
        name_out = 'Method';
    end
end

function print_summary(bundle1, bundle2, bscan_index, oa_stats_m1, oa_stats_m2, oa_neighbor_m1, oa_neighbor_m2, ascan_stats, ascan_x, outputs)
    deg = char(176);
    mu = char(181);
    unit_str = deg;

    fprintf('\nOA local consistency:\n');
    fprintf('%s:\n', oa_stats_m1.method_name);
    fprintf('  Median local variability = %.3f %s\n', oa_stats_m1.median_local_variability, deg);
    fprintf('  P90 local variability    = %.3f %s\n', oa_stats_m1.p90_local_variability, deg);
    fprintf('  Mean neighbor change     = %.3f %s\n', oa_neighbor_m1.mean_neighbor_change, deg);
    fprintf('  Mean local coherence R   = %.3f\n', oa_stats_m1.mean_local_coherence_R);

    fprintf('\n%s:\n', oa_stats_m2.method_name);
    fprintf('  Median local variability = %.3f %s\n', oa_stats_m2.median_local_variability, deg);
    fprintf('  P90 local variability    = %.3f %s\n', oa_stats_m2.p90_local_variability, deg);
    fprintf('  Mean neighbor change     = %.3f %s\n', oa_neighbor_m2.mean_neighbor_change, deg);
    fprintf('  Mean local coherence R   = %.3f\n', oa_stats_m2.mean_local_coherence_R);

    fprintf('\nInterpretation:\n');
    fprintf('  Lower variability / neighbor change means better local OA consistency.\n');
    fprintf('  Higher local coherence R means better local OA consistency.\n');

    fprintf('\nPhR A-scan:\n');
    fprintf('  X = %d (BScan %d)\n', ascan_x, bscan_index);
    fprintf('  Unit = %s\n', unit_str);
    fprintf('  MAE = %.3f %s\n', ascan_stats.mae, unit_str);
    fprintf('  RMSE = %.3f %s\n', ascan_stats.rmse, unit_str);
    fprintf('  Pearson r = %.4f\n', ascan_stats.pearson_r);

    fprintf('\nFigures saved to:\n');
    disp(outputs);
end
