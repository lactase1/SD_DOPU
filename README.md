# SSDOPU PS-OCT Pipeline

![MATLAB](https://img.shields.io/badge/MATLAB-R2020b+-blue.svg)
![Python](https://img.shields.io/badge/Python-3.7+-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

PS-OCT 数据处理流。输入 `.oct` 数据，输出结构图、DOPU、光轴/双折射与相位延迟结果，支持 DICOM、TIFF 与统计 `.mat` 导出。

## Overview

- 分光谱、空间、组合三种 DOPU 计算路径
- LA/OA、PhR 偏振参数计算与彩色编码
- 分批流式处理与并行计算，适配大体量数据
- 展平 En-face 与后巩膜边界展平（可选）
- 批处理输出组织与统计数据导出

## Quick Start

以下步骤用于最小化配置并直接产出结果。

### 1) 环境

- MATLAB R2020b 或更高版本
- 推荐安装 Parallel Computing Toolbox
- Python 3.7+（仅用于 scripts 中辅助脚本）

### 2) 数据

- 必需：待处理 `.oct` 文件（同一目录）
- 建议：表面边界 `.mat`（`topLines`，通常与 `.oct` 同名）
- 可选：后巩膜边界 `.mat`（放在 `scaler_mat/`，或在配置中指定路径）
- 说明：缺失边界文件时程序可自动估计边界，但精度通常低于人工/离线分割结果

### 3) 修改输入输出路径

编辑 `split_dopu_psoct_main.m`：

```matlab
data_path   = 'D:\path\to\oct_data';
output_base = 'D:\path\to\output';
```

### 4) 设置基础参数

编辑 `function/config_params.m`：

```matlab
params.tiff.saveDicom = 1;                   % 主输出
params.tiff.make_tiff = 0;                   % 首次建议关闭以加快验证
params.processing.max_frames = 0;            % 0=全帧；调试可设 20
params.parallel.maxWorkers = 0;              % 0=自动选择
params.processing.enable_flatten_enface = 1; % 生成展平与 En-face
params.processing.hasSeg = 1;                % 1=使用已有表面边界；0=自动估计
params.files.use_sclera_boundary = 0;        % 1=启用后巩膜边界；0=禁用
```

### 5) 运行

在 MATLAB 当前目录切换到项目根目录后执行：

```matlab
split_dopu_psoct_main
```

### 6) 输出检查

`output_base` 典型结构：

```text
output_base/
  run_YYYYMMDD_HHMMSS/      # 多文件批处理时创建
    sample_xxx/
      dcm/
      tiff/                 # make_tiff=1 时生成
      *_OA_stats.mat
      *_PhR_stats.mat
      *_Metadata_stats.mat
```

## Run Modes

### 单文件与多文件

- 程序自动扫描 `data_path` 下全部 `.oct`
- 单文件直接输出到 `output_base`
- 多文件输出到 `run_时间戳` 子目录

### 快速试跑

```matlab
params.processing.max_frames = 20;
params.tiff.make_tiff = 0;
```

### 开启 TIFF

```matlab
params.tiff.make_tiff = 1;
params.tiff.tiff_frame = 10;
```

## Repository Structure

### Root

- `split_dopu_psoct_main.m`：主入口，负责扫描、调度、并行池与整体流程
- `README.md`：项目文档

### function/

- `config_params.m`：参数中心
- `overwrite_config_params.m`：参数覆写
- `calculateSplitSpectrumDOPU.m`：分光谱 DOPU
- `calculateSpatialDOPU.m`：空间 DOPU
- `calculateCombinedDOPU.m`：组合 DOPU
- `calLAPhRALL.m`：LA/PhR 核心计算
- `calOAC.m`：光轴计算
- `cumulativeQUV.m`：累积 Q/U/V
- `filterVectorFieldAdaptive.m`：自适应向量滤波
- `vWinAvgFiltOpt.m`：窗口平均滤波优化
- `surf_seg.m`：表面分割
- `FreeSpace_PSOCT_Optimized.m`：流程优化组件
- `convert_dcm_to_tiff_local.m`：DCM 转 TIFF
- `dicomwrite_with_lowquality.m`：高/低质量 DCM 写出
- `quColoring.m`：偏振彩色编码
- `print_progress.m`：进度显示

### scripts/

- `organize_results_by_oct.py`：按样本整理结果
- `dcm_to_gif.py`：DCM 转 GIF
- `reduce_images.py`：图像降采样/压缩
- `generate_structure_pngs.m`：结构图 PNG 批量导出
- `extract_dcm_frame.m`：提取 DCM 指定帧
- `parameter_sweep_script.m`：参数扫描
- `run_planefit_sweep.m`：平面拟合参数扫描
- `analyze_DOPU_PhR_scatter_stats.m`：DOPU/PhR 散点统计
- `analyze_phr_global_stability.m`：全局稳定性分析
- `analyze_reddisk_psoct.m`：特定数据分析
- `auto_discover_method_dirs.m`：自动发现方法目录
- `generate_hsv_oa.m`：OA HSV 可视化
- `find_first_match.m`：匹配查找工具
- `circular_stats_deg.m`：角度统计
- `data_helper.py`、`tiff.py`：Python 数据与 TIFF 工具
- `polarpatch.m`：可视化工具

## Configuration

所有参数位于 `function/config_params.m`。

### 常用参数

- `params.processing.max_frames`：处理帧数上限，`0` 为全量
- `params.parallel.maxWorkers`：并行 worker 上限，`0` 为自动
- `params.parallel.batchSize`：分批大小
- `params.dopu.ss_mode`：`overlap9` 或 `nonoverlap5`
- `params.tiff.saveDicom` / `params.tiff.make_tiff`：输出开关
- `params.processing.enable_flatten_enface`：展平与 En-face 开关
- `params.files.use_sclera_boundary`：后巩膜边界流程开关

### 后巩膜边界

- 首选自动查找：`<data_path>/scaler_mat/*.mat`
- 回退路径：`params.files.sclera_boundary_path`
- 变量名可自动识别，或设置 `params.files.sclera_boundary_var`

## Output

每个样本目录通常包含：

- `dcm/`：体数据、DOPU、LA/PhR、展平体、En-face
- `tiff/`：可选 TIFF 输出
- `*_OA_stats.mat`：OA 统计数据
- `*_PhR_stats.mat`：PhR 统计数据
- `*_Metadata_stats.mat`：结构图、DOPU、边界元数据

## FAQ

### 内存不足

- 降低 `params.parallel.maxWorkers`
- 降低 `params.parallel.batchSize`
- 临时限制 `params.processing.max_frames`

### 处理速度慢

- 关闭 TIFF：`params.tiff.make_tiff = 0`
- 先用小帧数验证流程
- 在硬件允许下提高并行 worker

### 无输出或结果异常

- 检查 `split_dopu_psoct_main.m` 中路径配置
- 确认输入目录存在 `.oct` 文件
- 确认 `params.tiff.saveDicom = 1`

## License

MIT License

## Author

Yongxin Wang
