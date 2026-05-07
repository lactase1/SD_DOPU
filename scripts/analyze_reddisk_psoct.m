% analyze_reddisk_psoct.m
% 稳定重构版：统一 Mask 驱动 + 文献级绘图 (对标 Fig 10)
%
% 修复记录 (2026-02-21):
%   - 统一 Mask 驱动：画线、提取、统计全部从同一个 region_masks 矩阵获取
%   - 修复 Section 5 两个版本代码的错误合并（重复变量、语法碎片、未闭合 if）
%   - 修复 draw_half_color_wheel 重复定义（合并为单一函数）
%   - 补全 draw_stability_violin 函数签名及初始化代码
%   - 直方图极坐标不再旋转 Global_Ang_Offset（OA 是组织固有角度 0-180 deg）
%   - 图 (b)(c) 保留原图尺寸和坐标，不裁剪不居中

clear; close all; clc;

%% ==================== 1. 用户配置区域 ====================
% --- 数据路径配置 ---
input_root_dir = 'G:\1-Project\2023 王永鑫\Data\05_1310_redDisk';
input_dir2 = '';  % 留空则自动搜索第二个方法

output_dir = 'G:\1-Project\2023 王永鑫\Data\05_1310_redDisk\Graph_res';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

method1_name = 'DDG';
method2_name = 'Bayesian optimization';

% --- 几何与 ROI 配置 ---
Disk_Center  = [252, 252]; % 圆心坐标 [X, Y] (图像列, 行)
Radius_Inner = 60;
Radius_Outer = 180;

% --- 【新增】ROI 绘图基准方法 ---
% 1 = 使用 Method 1 结果作为基底（默认，与历史行为一致）
% 2 = 使用 Method 2 结果作为基底（仅在 has_method2 为 true 时有效）
roi_base_method = 2;  % 可改为 2 切换基准

% 全局零度偏移量（度）。
% 控制扇区 R1 (0 deg) 在图像中的起始方向。
% 135 = 副对角线方向（左下角），角度沿图像顺时针递增：
%   0 deg = 左下, 45 deg = 左, 90 deg = 左上, 135 deg = 上, 180 deg = 右上
Global_Ang_Offset = 135;

% 【新增】PS-OCT 系统固有的零点偏振相位误差补偿（用于对齐物理目标）
System_OA_Offset = 72; % 根据小提琴误差图，中心偏移约在 72 度

% 定义扇区：{名称, [起始角, 终止角], 目标角度}
% 角度在 Theta_map 坐标系中，顺时针从左下角开始
Regions = {
    'R1', [0,   30],   0;
    'R2', [30,  60],  30;
    'R3', [60,  90],  60;
    'R4', [90,  120], 90;
    'R5', [120, 150], 120;
    'R6', [150, 180], 150;
};

%% ==================== 2. 自动路径解析 ====================
fprintf('开始自动解析输入目录...\n');
required_patterns = {'*_OA_stats.mat', '*_PhR_stats.mat'};

if isempty(input_root_dir) || ~exist(input_root_dir, 'dir')
    error('input_root_dir is empty or does not exist.');
end

discovered = auto_discover_method_dirs(input_root_dir, required_patterns);
if isempty(discovered)
    error('在根目录下未找到可用 OA/PhR 数据目录: %s', input_root_dir);
end

dir_method1 = discovered(1).data_dir;
fprintf('自动识别 Method 1: %s\n', dir_method1);

has_method2 = false;
if ~isempty(input_dir2) && exist(input_dir2, 'dir')
    dir_method2 = input_dir2;
    has_method2 = true;
    fprintf('使用指定的 Method 2: %s\n', dir_method2);
elseif numel(discovered) >= 2
    dir_method2 = discovered(2).data_dir;
    has_method2 = true;
    fprintf('自动识别 Method 2: %s\n', dir_method2);
else
    fprintf('仅发现一组数据，进入单方法模式。\n');
end

%% ==================== 3. 数据加载与统一掩膜生成 ====================
fprintf('\n加载数据中...\n');
[OA_deg_1, PhR_1, ColorMap_1, H, W] = load_data_robust(dir_method1);

if has_method2
    [OA_deg_2, PhR_2, ColorMap_2, H2, W2] = load_data_robust(dir_method2);
    if H ~= H2 || W ~= W2
        warning('两种方法的图像尺寸不同: Method1=%dx%d, Method2=%dx%d', H, W, H2, W2);
    end
end

% --- 构建极坐标网格 ---
[xg, yg] = meshgrid(1:W, 1:H);
Xrel = xg - Disk_Center(1);
Yrel = yg - Disk_Center(2);
R_map = sqrt(Xrel.^2 + Yrel.^2);

% Theta_map: 空间极坐标角度（含 Global_Ang_Offset）
% 使用 atan2d(Yrel, Xrel)（不翻转 Y）使角度在图像坐标系中沿顺时针方向递增
% 减去 Global_Ang_Offset 将 0 deg 对齐到副对角线（左下方向）
Theta_map = mod(atan2d(Yrel, Xrel) - Global_Ang_Offset, 360);

% --- 统一生成所有区域的掩膜矩阵 ---
% region_masks(r,c) = 区域编号 (1~N)，0 表示不属于任何 ROI
num_regions = size(Regions, 1);
region_masks = zeros(H, W);           % 每个像素归属的区域编号
global_roi_mask = false(H, W);        % 全局 ROI 并集
stats_1 = cell(num_regions, 1);
stats_2 = cell(num_regions, 1);
targets = zeros(num_regions, 1);

r_mask = (R_map >= Radius_Inner) & (R_map <= Radius_Outer);

for i = 1:num_regions
    ang_range = Regions{i, 2};
    targets(i) = Regions{i, 3};

    % 角度掩膜（处理跨 0 deg / 360 deg 的情况）
    if ang_range(1) <= ang_range(2)
        ang_mask = (Theta_map >= ang_range(1)) & (Theta_map <= ang_range(2));
    else
        % 跨越 0 deg 的情况（如 [350, 10]）
        ang_mask = (Theta_map >= ang_range(1)) | (Theta_map <= ang_range(2));
    end

    curr_mask = ang_mask & r_mask;
    region_masks(curr_mask) = i;
    global_roi_mask = global_roi_mask | curr_mask;

    % 提取 OA 值（统一从掩膜驱动）
    oa1 = mod(OA_deg_1(curr_mask), 180);
    stats_1{i} = oa1(~isnan(oa1));

    if has_method2
        oa2 = mod(OA_deg_2(curr_mask), 180);
        stats_2{i} = oa2(~isnan(oa2));
    end
end

fprintf('掩膜生成完毕：共 %d 个区域, 总 ROI 像素 = %d\n', num_regions, nnz(global_roi_mask));

%% ==================== 4. 绘图：Fig 10 (a), (b), (c) ====================
% 所有绘图和数据提取均基于同一个 region_masks / global_roi_mask
fprintf('生成 ROI 及区域提取图...\n');

% --- 根据 roi_base_method 选择基准色图与名称 ---
if roi_base_method == 2 && ~has_method2
    warning('Requested ROI base method 2 but Method 2 data not available; falling back to Method 1.');
    roi_base_method = 1;
end
if roi_base_method == 1
    cmap_base = ColorMap_1;
    base_name = method1_name;
else
    cmap_base = ColorMap_2;
    base_name = method2_name;
end
fprintf('使用 "%s" 作为 ROI/区域提取基准\n', base_name);

% --- Fig 10(a): En-face 图 + ROI 白线边界（基于选定方法） ---
fig_a = figure('Name', 'Fig 10(a) ROI Overlay', 'Color', 'w');
imshow(cmap_base); hold on;
for i = 1:num_regions
    % 用 region_masks 的轮廓来画 ROI 边界（与数据提取完全一致）
    bw_i = (region_masks == i);
    bnd = bwboundaries(bw_i, 'noholes');
    for b = 1:numel(bnd)
        plot(bnd{b}(:,2), bnd{b}(:,1), 'w', 'LineWidth', 3.5);
    end
    % 放置标签
    [row_idx, col_idx] = find(bw_i);
    if ~isempty(row_idx)
        text(mean(col_idx), mean(row_idx), Regions{i,1}, ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 16, ...
            'HorizontalAlignment', 'center');
    end
end
title(sprintf('(a) %s enface with ROI', base_name));
saveas(fig_a, fullfile(output_dir, sprintf('Fig10_a_ROI_Overlay_%s.png', base_name)));

% --- Fig 10(b): 基准方法的提取区域（保留原图尺寸） ---
fig_b = figure('Name', 'Fig 10(b) Base Method Extracted', 'Color', 'k');
masked_base = bsxfun(@times, cmap_base, cast(global_roi_mask, class(cmap_base)));
imshow(masked_base); title(base_name, 'Color', 'w');
saveas(fig_b, fullfile(output_dir, sprintf('Fig10_b_%s_Extracted.png', base_name)));

% --- Fig 10(c): 另一种方法的提取区域（如果有） ---
if has_method2
    % 非基准方法作为 Fig10(c)
    if roi_base_method == 1
        other_cmap = ColorMap_2;
        other_name = method2_name;
    else
        other_cmap = ColorMap_1;
        other_name = method1_name;
    end
    fig_c = figure('Name', 'Fig 10(c) Other Method Extracted', 'Color', 'k');
    masked_other = bsxfun(@times, other_cmap, cast(global_roi_mask, class(other_cmap)));
    imshow(masked_other); title(other_name, 'Color', 'w');
    saveas(fig_c, fullfile(output_dir, sprintf('Fig10_c_%s_Extracted.png', other_name)));
end

%% ==================== 5. 绘图：Fig 10(d) 极坐标直方图 ====================
% OA 角度是组织固有属性 (0-180 deg)，直方图无需旋转 Global_Ang_Offset
fprintf('生成 Fig 10(d) 分布直方图...\n');
fig_d = figure('Name', 'Fig 10(d) Polar Histograms', 'Color', 'w', ...
    'Position', [100, 100, 1200, 600]);

% 线型/标记区分（适配黑白打印）
m1_style = struct('LineStyle', '-',  'Marker', 'o', 'Color', [0 0 0]);
m2_style = struct('LineStyle', '--', 'Marker', '^', 'Color', [0.5 0.5 0.5]);

edges = 0:5:180;
theta_centers = deg2rad((edges(1:end-1) + edges(2:end)) / 2);
rows_d = 2; cols_d = ceil(num_regions / 2);

for i = 1:num_regions
    ax = subplot(rows_d, cols_d, i);
    hold(ax, 'on'); axis(ax, 'equal');
    set(ax, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Color', 'w');

    % 绘制半圆色轮背景（标准 0-180 deg 方向）
    draw_half_color_wheel(ax);

    % 计算直方图
    counts1 = histcounts(stats_1{i}, edges, 'Normalization', 'probability');
    max_rho = max(counts1);

    if has_method2
        counts2 = histcounts(stats_2{i}, edges, 'Normalization', 'probability');
        max_rho = max([max_rho, max(counts2)]);
    end

    if max_rho <= 0, max_rho = 1e-6; end
    sfac = 0.9 / max_rho;

    % 绘制 Method 1（黑线 + 圆点）
    c1_scaled = counts1 * sfac;
    [x1, y1] = pol2cart(theta_centers, c1_scaled);
    x1 = [0, x1, 0]; y1 = [0, y1, 0];
    plot(ax, x1, y1, 'LineStyle', m1_style.LineStyle, 'Marker', m1_style.Marker, ...
        'Color', m1_style.Color, 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'MarkerFaceColor', m1_style.Color);

    % 绘制 Method 2（灰虚线 + 方块）
    if has_method2
        c2_scaled = counts2 * sfac;
        [x2, y2] = pol2cart(theta_centers, c2_scaled);
        x2 = [0, x2, 0]; y2 = [0, y2, 0];
        plot(ax, x2, y2, 'LineStyle', m2_style.LineStyle, 'Marker', m2_style.Marker, ...
            'Color', m2_style.Color, 'LineWidth', 1.5, 'MarkerSize', 4, ...
            'MarkerFaceColor', m2_style.Color);
    end

    title(ax, Regions{i,1}, 'FontSize', 12, 'FontWeight', 'bold');
    xlim(ax, [-1.1 1.1]); ylim(ax, [-0.1 1.1]);

    % 角度刻度标签
    for tg = [0 30 60 90 120 150 180]
        [tx, ty] = pol2cart(deg2rad(tg), 1.05);
        text(ax, tx, ty, sprintf('%d', tg), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end

    % —— 保存当前子图到单独文件 ——
    fig_single = figure('Visible','off','Color','w','Position',[100,100,400,400]);
    ax2 = axes('Parent', fig_single);
    hold(ax2, 'on'); axis(ax2, 'equal');
    set(ax2, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Color', 'w');
    draw_half_color_wheel(ax2);
    % repeat plotting for this one region
    plot(ax2, x1, y1, 'LineStyle', m1_style.LineStyle, 'Marker', m1_style.Marker, ...
        'Color', m1_style.Color, 'LineWidth', 1.5, 'MarkerSize', 4, ...
        'MarkerFaceColor', m1_style.Color);
    if has_method2
        plot(ax2, x2, y2, 'LineStyle', m2_style.LineStyle, 'Marker', m2_style.Marker, ...
            'Color', m2_style.Color, 'LineWidth', 1.5, 'MarkerSize', 4, ...
            'MarkerFaceColor', m2_style.Color);
    end
    title(ax2, Regions{i,1}, 'FontSize', 12, 'FontWeight', 'bold');
    xlim(ax2, [-1.1 1.1]); ylim(ax2, [-0.1 1.1]);
    for tg = [0 30 60 90 120 150 180]
        [tx, ty] = pol2cart(deg2rad(tg), 1.05);
        text(ax2, tx, ty, sprintf('%d', tg), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
    outname_single = fullfile(output_dir, sprintf('Fig10_d_region_%s.png', Regions{i,1}));
    saveas(fig_single, outname_single);
    close(fig_single);
end

sgtitle('Fig 10(d) OA Distribution Histograms', 'FontWeight', 'bold', 'FontSize', 14);
% 使用真实线型+标记的图例（避免文本符号显示异常）
ax_leg = axes('Parent', fig_d, 'Position', [0 0 1 1], 'Visible', 'off');
hold(ax_leg, 'on');
h_leg1 = plot(ax_leg, NaN, NaN, 'LineStyle', m1_style.LineStyle, 'Marker', m1_style.Marker, ...
    'Color', m1_style.Color, 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerFaceColor', m1_style.Color);
if has_method2
    h_leg2 = plot(ax_leg, NaN, NaN, 'LineStyle', m2_style.LineStyle, 'Marker', m2_style.Marker, ...
        'Color', m2_style.Color, 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', m2_style.Color);
    legend([h_leg1, h_leg2], {method1_name, method2_name}, ...
        'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'FontSize', 16, 'FontWeight', 'bold', 'Box', 'off');
else
    legend(h_leg1, method1_name, 'Location', 'southoutside', ...
        'FontSize', 16, 'FontWeight', 'bold', 'Box', 'off');
end

saveas(fig_d, fullfile(output_dir, 'Fig10_d_Polar_Histograms.png'));
fprintf('所有处理完毕！结果存放在: %s\n', output_dir);

%% ==================== 6. 绘图：小提琴稳定性评估图 ====================
fprintf('生成稳定性小提琴图...\n');
draw_stability_violin(stats_1, targets, Regions(:,1), method1_name, output_dir, 'Fig11_Stability_Method1');

if has_method2
    draw_stability_violin(stats_2, targets, Regions(:,1), method2_name, output_dir, 'Fig11_Stability_Method2');
end

%% ==================== 7. 诊断验证图 ====================
fprintf('生成掩膜诊断图...\n');
fig_diag = figure('Name', 'Diagnostic: Region Mask', 'Color', 'w', ...
    'Position', [100 100 800 400]);

subplot(1,2,1);
imagesc(region_masks); axis image; colorbar;
title('Region Mask (1-6)');
colormap(gca, [0 0 0; lines(num_regions)]);
caxis([0 num_regions]);

subplot(1,2,2);
imshow(ColorMap_1); hold on;
cmap_regions = lines(num_regions);
for i = 1:num_regions
    bw_i = (region_masks == i);
    bnd = bwboundaries(bw_i, 'noholes');
    for b = 1:numel(bnd)
        plot(bnd{b}(:,2), bnd{b}(:,1), 'Color', cmap_regions(i,:), 'LineWidth', 3.5);
    end
end
title('Overlay Verification');
legend(Regions{:,1}, 'Location', 'best', 'Box', 'off');

saveas(fig_diag, fullfile(output_dir, 'Diagnostic_Region_Mask.png'));
fprintf('诊断图已保存。\n');

%% ==================== 辅助功能函数 ====================

function [OA_deg, PhR, OA_Color_Map, H, W] = load_data_robust(dir_path)
    % 注意：这里通过 evalin 获取主工作区的 System_OA_Offset 变量
    try 
        sys_offset = evalin('base', 'System_OA_Offset'); 
    catch
        sys_offset = 0; 
    end

    oa_files = dir(fullfile(dir_path, '*_OA_stats.mat'));
    phr_files = dir(fullfile(dir_path, '*_PhR_stats.mat'));
    if isempty(oa_files) || isempty(phr_files)
        error('路径 %s 下缺失必要的 mat 文件', dir_path);
    end

    OA_data = load(fullfile(dir_path, oa_files(1).name));
    PhR_data = load(fullfile(dir_path, phr_files(1).name));

    if isfield(OA_data, 'OA_enface_avg'), OA_raw = OA_data.OA_enface_avg;
    else, OA_raw = OA_data.OA_enface_single; end
    if isfield(PhR_data, 'PhR_enface_avg'), PhR = PhR_data.PhR_enface_avg;
    else, PhR = PhR_data.PhR_enface_single; end

    OA_arr = squeeze(OA_raw);
    Q = double(OA_arr(:,:,1)); U = double(OA_arr(:,:,2));
    
    % 【核心修正 1：应用系统的固有偏振补偿】
    OA_deg_raw = 0.5 * atan2d(U, Q);
    OA_deg = mod(OA_deg_raw - sys_offset, 180);

    % 【核心修正 2：抛弃假彩色图，用补偿后的真实数据自己生成HSV图】
    % 不再读取 OA_data.OA_enface_rgb
    H_val = OA_deg / 180;    
    S_val = ones(size(H_val));
    V_val = ones(size(H_val));
    OA_Color_Map = hsv2rgb(cat(3, H_val, S_val, V_val));
    
    [H, W] = size(OA_deg);
end

function draw_half_color_wheel(ax)
% DRAW_HALF_COLOR_WHEEL  在坐标轴 ax 上绘制标准 0-180 deg 半圆 HSV 色轮
    [TH, R] = meshgrid(linspace(0, pi, 180), linspace(0, 1, 100)');
    [XB, YB] = pol2cart(TH, R);
    C_bg = hsv2rgb(cat(3, TH / pi, ones(size(TH)), ones(size(TH))));
    surf(ax, XB, YB, zeros(size(XB)), C_bg, ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    % 同心圆参考线
    for rt = 0.2:0.2:1
        [xcg, ycg] = pol2cart(linspace(0, pi, 180), rt);
        plot(ax, xcg, ycg, 'Color', [1 1 1 0.4], 'LineWidth', 0.5);
    end
    % 径向参考线
    for thl = deg2rad([30 60 90 120 150])
        [xl, yl] = pol2cart(thl, [0 1]);
        plot(ax, xl, yl, 'Color', [1 1 1 0.4], 'LineWidth', 0.5);
    end
end

function draw_stability_violin(stats, targets, region_names, method_label, output_dir, file_prefix)
% DRAW_STABILITY_VIOLIN  绘制 OA 角度误差的小提琴图
    nS = numel(stats);
    x = 1:nS;
    maxWidth = 0.35;

    % 配色
    palette = lines(nS);

    fig = figure('Name', ['Violin: ' method_label], 'Color', 'w', ...
        'Position', [100 100 900 500]);
    hold on;

    % 零误差参考线
    yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

    for i = 1:nS
        % 计算误差：映射到 [-90, 90] 区间
        vals = stats{i} - targets(i);
        vals = mod(vals + 90, 180) - 90;
        n = numel(vals);
        if n == 0, continue; end

        % 核密度估计 (KDE)
        xi_query = linspace(-90, 90, 200);
        try
            [f, xi] = ksdensity(vals, xi_query);
        catch
            [cts, edges_v] = histcounts(vals, 60, 'Normalization', 'pdf');
            xi = (edges_v(1:end-1) + edges_v(2:end)) / 2;
            f = smooth(cts, 5)';
        end

        if max(f) <= eps, continue; end
        f_norm = f / max(f) * maxWidth;
        xv = [x(i) - f_norm, fliplr(x(i) + f_norm)];
        yv = [xi, fliplr(xi)];

        % 小提琴斑块
        baseColor = palette(i,:);
        patch(xv, yv, baseColor, 'FaceAlpha', 0.55, ...
            'EdgeColor', baseColor * 0.5, 'LineWidth', 1.2);

        % 散点抖动
        displayN = min(300, n);
        idx = randperm(n, displayN);
        scatter_vals = vals(idx);
        jitter = zeros(displayN, 1);
        for j = 1:displayN
            [~, closest_idx] = min(abs(xi - scatter_vals(j)));
            jitter(j) = (rand - 0.5) * 2 * (f_norm(closest_idx) * 0.7);
        end
        scatter(x(i) + jitter, scatter_vals, 10, ...
            'MarkerFaceColor', baseColor * 0.7, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.35);

        % 四分位距 & 中位数
        med = median(vals);
        q1 = prctile(vals, 25);
        q3 = prctile(vals, 75);
        plot([x(i), x(i)], [q1, q3], '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 4);
        plot(x(i), med, 'o', 'MarkerSize', 7, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

        % 统计文本
        theta_circ = deg2rad(2 * vals);
        C = mean(cos(theta_circ));
        S = mean(sin(theta_circ));
        Rval = sqrt(C^2 + S^2);
        mean_ang = mod(rad2deg(atan2(S, C) / 2), 180);
        txt = sprintf('n=%d\n\\mu_c=%.1f^{\\circ}\nR=%.2f', n, mean_ang, Rval);
        text(x(i) + 0.22, -65, txt, 'FontSize', 8, ...
            'Color', [0.2 0.2 0.2], 'Interpreter', 'tex');
    end

    xlim([0.4 nS + 0.6]);
    ylim([-90 90]);
    set(gca, 'XTick', 1:nS, 'XTickLabel', region_names, ...
        'TickDir', 'in', 'LineWidth', 2.0, 'TickLength', [0.015 0.015]);
    ylabel('Angle Error (deg)', 'FontWeight', 'bold');
    title(sprintf('[%s] OA Stability Evaluation', method_label), 'FontWeight', 'bold');
    grid on; box on;
    set(gca, 'Layer', 'top');
    saveas(fig, fullfile(output_dir, [file_prefix '.png']));
end
