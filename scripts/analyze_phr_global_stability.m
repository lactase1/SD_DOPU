% analyze_phr_global_stability.m
%
% PS-OCT PhR 整体稳定性分析 (Global Stability) — 支持双方法对比
% 核心功能：
% 1. 读取 3D PhR 数据，提取指定深度 (En-face) 切片
% 2. 假设全图均为有效数据 (无 Mask)，计算整体均值与标准差
% 3. 绘制 SCI 级 "直方图 + KDE 密度曲线"，两种方法用不同颜色对比
% 4. 兼容单方法模式 (不提供 input_dir2 时退化为原始行为)
%
% 输出：Analysis_Global_PhR.csv 和 PhR_Global_Histogram.png

clear; close all; clc;

%% --------------------- 1. 用户配置 (Configuration) ---------------------
% 输入目录 — Method 1 (not used — script discovers from root)
% (root-only discovery — set input_root_dir below)

% 根目录（必需）：脚本将仅从此根目录递归查找方法数据
% 示例：D:\1-Liu Jian\yongxin.wang\Output\Graph\Mat
input_root_dir = 'G:\1-Project\2023 王永鑫\Data\05_1310_redDisk';

% 输入目录 — Method 2 (留空 '' 则退化为单方法模式)
input_dir2 = '';

% 方法名称 (用于图例和标注)
method1_name = 'DDG';
method2_name = 'Bayesian optimization'; % 以前叫 SSDOPU++

% 输出目录
output_dir = 'G:\1-Project\2023 王永鑫\Data\05_1310_redDisk\Graph_res';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% === 切片参数 ===
Target_Depth_Index = 50;  % 与 tiff_frame=10 保持一致

% === 单位换算参数 ===
phr_scale = 180/pi;  % rad → °（直接换算，不做范围截断）

%% --------------------- 2. 数据加载与处理 ---------------------
required_patterns = {'*_PhR_stats.mat'};

% Use the configured root directory and always auto-discover method folders from it.
if isempty(input_root_dir) || ~exist(input_root_dir, 'dir')
    error('input_root_dir is empty or does not exist. Set it to the Mat root directory (e.g. D:\\...\\Output\\Graph\\Mat).');
end

discovered = auto_discover_method_dirs(input_root_dir, required_patterns);
if isempty(discovered)
    error('在根目录下未找到可用的 PhR 数据目录: %s', input_root_dir);
end
input_dir = discovered(1).data_dir;
if strcmp(method1_name, 'Method 1')
    method1_name = discovered(1).method_name;
end
if isempty(input_dir2) && numel(discovered) >= 2
    input_dir2 = discovered(2).data_dir;
    if strcmp(method2_name, 'Method 2')
        method2_name = discovered(2).method_name;
    end
end
fprintf('自动识别 Method 1: %s\n', input_dir);
if ~isempty(input_dir2)
    fprintf('自动识别 Method 2: %s\n', input_dir2);
end

if ~isempty(input_dir2) && exist(input_dir2, 'dir') && ~has_required_files_local(input_dir2, required_patterns)
    discovered2 = auto_discover_method_dirs(input_dir2, required_patterns);
    if ~isempty(discovered2)
        input_dir2 = discovered2(1).data_dir;
        if strcmp(method2_name, 'Method 2')
            method2_name = discovered2(1).method_name;
        end
        fprintf('Method 2 根目录已解析到: %s\n', input_dir2);
    else
        error('Method 2 路径下未找到可用 PhR 数据目录: %s', input_dir2);
    end
end

has_method2 = ~isempty(input_dir2) && exist(input_dir2, 'dir');

% --- Method 1 ---
fprintf('=== %s ===\n', method1_name);
[PhR_Values_1, PhR_Frame_1] = load_phr_slice(input_dir, Target_Depth_Index, phr_scale);

% --- Method 2 (if provided) ---
if has_method2
    fprintf('\n=== %s ===\n', method2_name);
    [PhR_Values_2, PhR_Frame_2] = load_phr_slice(input_dir2, Target_Depth_Index, phr_scale);
else
    PhR_Values_2 = [];
    PhR_Frame_2 = [];
end

%% --------------------- 3. 统计计算 ---------------------
Global_Mean_1 = mean(PhR_Values_1);
Global_Std_1  = std(PhR_Values_1);
Pixel_Count_1 = length(PhR_Values_1);

fprintf('\n=== %s 统计结果 ===\n', method1_name);
fprintf('Mean: %.4f °/μm\n', Global_Mean_1);
    fprintf('Std : %.4f °/μm (越小越稳定)\n', Global_Std_1);

if has_method2
    Global_Mean_2 = mean(PhR_Values_2);
    Global_Std_2  = std(PhR_Values_2);
    Pixel_Count_2 = length(PhR_Values_2);
    
    fprintf('\n=== %s 统计结果 ===\n', method2_name);
    fprintf('Mean: %.4f °/μm\n', Global_Mean_2);
    fprintf('Std : %.4f °/μm (越小越稳定)\n', Global_Std_2);
end

% 保存 CSV
if has_method2
    T = table([Target_Depth_Index; Target_Depth_Index], ...
              [Global_Mean_1; Global_Mean_2], ...
              [Global_Std_1; Global_Std_2], ...
              [Pixel_Count_1; Pixel_Count_2], ...
              {method1_name; method2_name}, ...
              'VariableNames', {'Depth_Layer', 'PhR_Mean_deg_per_um', 'PhR_Std_deg_per_um', 'Pixel_Count', 'Method'});
else
    T = table(Target_Depth_Index, Global_Mean_1, Global_Std_1, Pixel_Count_1, ...
        'VariableNames', {'Depth_Layer', 'PhR_Mean_deg_per_um', 'PhR_Std_deg_per_um', 'Pixel_Count'});
end
csv_path = fullfile(output_dir, 'Analysis_Global_PhR.csv');
writetable(T, csv_path);
fprintf('统计表格已保存: %s\n', csv_path);

% ------------------ 额外输出：局部 std 热图 (每种方法分别生成) ------------------
generate_local_std_map(PhR_Frame_1, output_dir, method1_name, 'PhR_std_map');
if has_method2
    generate_local_std_map(PhR_Frame_2, output_dir, method2_name, 'PhR_std_map_method2');
end

%% --------------------- 4. 高级绘图 — 双方法对比 ---------------------
% ---- 配色 ----
% Method 1: 蓝色系 (保持与原版一致)
Color_Hist_1 = [75, 139, 190] / 255;   % 淡青蓝
Color_KDE_1  = [63, 104, 221] / 255;   % 深蓝
% Method 2: 橙红色系
Color_Hist_2 = [230, 159, 100] / 255;  % 淡橙
Color_KDE_2  = [210, 80, 50] / 255;    % 深红橙
% 公共
Color_Mean_Line = [154, 153, 153] / 255;

% 线型区分（适配黑白打印）
m1_style = struct('LineStyle', '-');
m2_style = struct('LineStyle', '--');

% 创建画布
fig = figure('Name', 'PhR Global Stability', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on;

% --- Method 1: 直方图 + KDE ---
h1 = histogram(PhR_Values_1, 100, 'Normalization', 'pdf');
h1.FaceColor = Color_Hist_1;
h1.EdgeColor = 'k';
h1.FaceAlpha = 0.5;
h1.LineStyle = m1_style.LineStyle;

[f1, xi1] = ksdensity(PhR_Values_1);
hKDE_1 = plot(xi1, f1, 'Color', Color_KDE_1, 'LineWidth', 2.5, 'LineStyle', m1_style.LineStyle);

% ±1σ 区间 (Method 1)
sigma_low_1 = Global_Mean_1 - Global_Std_1;
sigma_high_1 = Global_Mean_1 + Global_Std_1;
idx_shade_1 = xi1 >= sigma_low_1 & xi1 <= sigma_high_1;
if any(idx_shade_1)
    x_sh = xi1(idx_shade_1);
    y_sh = f1(idx_shade_1);
    fill([x_sh, fliplr(x_sh)], [y_sh, zeros(1,numel(y_sh))], Color_KDE_1, ...
        'EdgeColor','none', 'FaceAlpha', 0.12);
end

% --- Method 2: 直方图 + KDE (if available) ---
if has_method2
    h2 = histogram(PhR_Values_2, 100, 'Normalization', 'pdf');
    h2.FaceColor = Color_Hist_2;
    h2.EdgeColor = 'k';
    h2.FaceAlpha = 0.5;
    h2.LineStyle = m2_style.LineStyle;
    
    [f2, xi2] = ksdensity(PhR_Values_2);
    hKDE_2 = plot(xi2, f2, 'Color', Color_KDE_2, 'LineWidth', 2.5, 'LineStyle', m2_style.LineStyle);
    
    % ±1σ 区间 (Method 2)
    sigma_low_2 = Global_Mean_2 - Global_Std_2;
    sigma_high_2 = Global_Mean_2 + Global_Std_2;
    idx_shade_2 = xi2 >= sigma_low_2 & xi2 <= sigma_high_2;
    if any(idx_shade_2)
        x_sh2 = xi2(idx_shade_2);
        y_sh2 = f2(idx_shade_2);
        fill([x_sh2, fliplr(x_sh2)], [y_sh2, zeros(1,numel(y_sh2))], Color_KDE_2, ...
            'EdgeColor','none', 'FaceAlpha', 0.12);
    end
end

% === 图表美化 ===
if has_method2
    title('Global Phase Retardation Stability', 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Arial');
else
    title('Global Phase Retardation Stability', 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Arial');
end
xlabel('Phase Retardation (°/μm)', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Probability Density', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');

grid off;
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.FontWeight = 'bold';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.4;
ax.TickDir = 'in';
ax.LineWidth = 2.0;
ax.TickLength = [0.015 0.015];
ax.Box = 'off';

% X 轴范围 (考虑双方法数据) — use shared limit logic
if has_method2
    all_vals = [PhR_Values_1; PhR_Values_2];
else
    all_vals = PhR_Values_1;
end
[low_lim, high_lim] = compute_pr_axis_limits(all_vals);
xlim([low_lim, high_lim]);

% === 标注框 ===
if has_method2
    str_stats = {
        sprintf('--- %s ---', method1_name);
        sprintf('Mean: %.3f °/μm', Global_Mean_1);
        sprintf('Std : %.3f °/μm', Global_Std_1);
        sprintf('N   : %d px', Pixel_Count_1);
        '';
        sprintf('--- %s ---', method2_name);
        sprintf('Mean: %.3f °/μm', Global_Mean_2);
        sprintf('Std : %.3f °/μm', Global_Std_2);
        sprintf('N   : %d px', Pixel_Count_2)
    };
else
    str_stats = {
        sprintf('Global Statistics');
        sprintf('Mean: %.3f °/μm', Global_Mean_1);
        sprintf('Std : %.3f °/μm', Global_Std_1);
        sprintf('N   : %d px', Pixel_Count_1)
    };
end

% move statistics box down toward lower right corner but keep it
% slightly above the bottom edge so it's not flush with the axis
% adjust values until the position feels natural

dim = [0.6 0.39 0.25 0.20];
annotation('textbox', dim, 'String', str_stats, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'none', ...
    'FontName', 'Arial', ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'LineWidth', 0.5, ...
    'Interpreter', 'none');

% === 图例 ===
if has_method2
    % put legend in the upper right corner inside the axes
    hLeg = legend([hKDE_1, hKDE_2], {method1_name, method2_name}, ...
        'Location', 'northeast', 'FontSize', 16, 'FontWeight', 'bold', 'Box', 'off');
    % if it still overlaps the annotation box, tweak dimensions there instead
end

% === 保存图片 ===
set(gcf, 'PaperPositionMode', 'auto');
img_filename = fullfile(output_dir, 'PhR_Global_Histogram.png');
print(fig, img_filename, '-dpng', '-r300');
fprintf('高清图表已保存: %s\n', img_filename);


%% ========================== Helper Functions ==========================

function [PhR_Values, PhR_Frame] = load_phr_slice(dir_path, depth_idx, phr_scale)
% Load PhR data, extract a single en-face slice, convert rad → degrees.
    if nargin < 3, phr_scale = 180/pi; end
    fprintf('正在搜索 PhR 数据文件...\n');
    list = dir(fullfile(dir_path, '*_PhR_stats.mat'));
    if isempty(list)
        error('未在目录中找到 *_PhR_stats.mat 文件: %s', dir_path);
    end
    phr_filename = fullfile(list(1).folder, list(1).name);
    fprintf('加载文件: %s\n', list(1).name);
    
    data = load(phr_filename);
    
    if isfield(data, 'PhR_phys_flat')
        PhR_3D = data.PhR_phys_flat;
    else
        error('文件中未找到变量 PhR_phys_flat');
    end
    
    [dim1, ~, ~] = size(PhR_3D);
    if depth_idx > dim1
        error('目标层数 (%d) 超出了数据深度范围 (%d)', depth_idx, dim1);
    end
    
    PhR_Frame = squeeze(PhR_3D(depth_idx, :, :)) * phr_scale;  % rad → °
    PhR_Frame = max(PhR_Frame, 0);  % 仅截断负値
    
    if ismatrix(PhR_Frame)
        fprintf('成功提取第 %d 层切片，尺寸: %d x %d\n', depth_idx, size(PhR_Frame,1), size(PhR_Frame,2));
    else
        error('提取的切片不是 2D 矩阵，请检查数据维度。');
    end
    
    valid_mask = ~isnan(PhR_Frame);
    PhR_Values = PhR_Frame(valid_mask);
    
    if isempty(PhR_Values)
        error('该层数据为空或全部为 NaN，请检查 depth_idx 是否正确。');
    end
end

function generate_local_std_map(PhR_Frame, output_dir, method_label, file_prefix)
% 计算并保存局部 std 热图
    LocalStdWindow = 31;
    ws = LocalStdWindow;
    kernel = ones(ws, ws);
    
    valid_mask = ~isnan(PhR_Frame);
    count_map = conv2(double(valid_mask), kernel, 'same');
    sum_map = conv2(double(PhR_Frame .* valid_mask), kernel, 'same');
    sumsq_map = conv2(double((PhR_Frame.^2) .* valid_mask), kernel, 'same');
    
    local_mean = sum_map ./ max(count_map, eps);
    local_mean_sq = sumsq_map ./ max(count_map, eps);
    local_var = max(local_mean_sq - local_mean.^2, 0);
    PhR_std_map = sqrt(local_var);
    PhR_std_map(count_map == 0) = NaN;
    
    save(fullfile(output_dir, [file_prefix '.mat']), 'PhR_std_map', 'LocalStdWindow');
    std_vals = PhR_std_map(:);
    std_vals = std_vals(isfinite(std_vals));
    if ~isempty(std_vals)
        clim_upper = prctile(std_vals, 95);
        clim_upper = max(clim_upper, nanmedian(std_vals) + 2*nanstd(std_vals));
        fh = figure('Visible','off','Color','w','Position',[100,100,700,600]);
        imagesc(PhR_std_map); axis image off; colormap(jet); colorbar;
        caxis([0 clim_upper]);
        title(sprintf('[%s] Local PhR std (window=%d px) — median=%.3f °/μm', method_label, ws, nanmedian(std_vals)), 'FontSize', 12);
        print(fh, fullfile(output_dir, [file_prefix '.png']), '-dpng', '-r300');
        close(fh);
        fprintf('已保存局部 std 热图 (%s) -> %s\n', method_label, fullfile(output_dir, [file_prefix '.png']));
    else
        fprintf('局部 std 计算无有效像素 (%s)，跳过热图保存。\n', method_label);
    end
end

function [low, high] = compute_pr_axis_limits(values)
% Shared helper used by both global stability and scatter scripts.
    p01 = prctile(values, 0.5);
    p99 = prctile(values, 99.5);
    margin = (p99 - p01) * 0.1;
    low = max(0, p01 - margin);
    high = p99 + margin;
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
