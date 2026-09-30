# 统一斑图初值读取约定

反向扫描的斑图初值与动力学噪声种子不是同一个概念。噪声重复实验仍然运行；取消的是额外的斑图初值训练。

## 标准斑图种子

- 文件：`results/evolution_er_N{N}_K{K}_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat`。
- 同一 N、K 的 WS、BA、ER、混合配置以及不同 alpha、beta、噪声参数共用此文件，不因参数改变重新训练。
- `load_standard_pattern_seed(cfg)` 仅读取。文件缺失、长度不等于 `2*N*K`、终态存在非有限值或未形成明显斑图时，明确报错。
- 实验入口可以调用 `load_standard_pattern_seed(cfg, true)`：仅当 N 不等于基准 200 且文件缺失时，允许生成一次该 N、K 的标准 ER 种子。训练参数固定为 alpha=0.05、beta=0.005、sigma=100、noise=0；不继承目标实验参数。
- K 改变也必须匹配状态尺寸。N=200 的 K 扫描不会自动生成种子，需要提前准备相应 K 的文件，不能截取或拼接五层状态。
- `pattern_evolution(..., false)` 仍保留为显式的正向演化/种子准备入口，但不再被上述实验脚本隐式调用来训练 N=200 的种子。

## 已调整的实验入口

- `fig2c_topology_mean_std`：取消逐配置、逐耦合比参考态训练，只计算谱阈值；所有反向搜索使用同一个读取的标准种子。使用 `scan_sharedseed_*.mat` 和带 `sharedseed` 标记的图片，旧版训练初值结果不会被混用或自动覆盖。
- `fig1a`、`fig1cd`、`fig1b_seed`、`fig3a/b/c`、`fig3c_permuted`、`fig4a/b/c`、`fig5a/b/c`：统一使用上述读取规则。
- `sigma_beta_alpha`、`sigma_alpha_beta`、`sigma_beta_noise_er`、`sigma_beta_connectivity_er`、`sigma_eta_connectivity_er`、两个 `phase_diagram_*`：统一读取。
- `find_thresholds`、`sweep_param`、`sweep_param_ws`：保持已有函数参数和阈值/积分设置，改为统一读取；`find_thresholds` 仍支持显式传入 `cfg.y_seed`，但会检查长度及有限值。
- `fig_hetero_comparison`：读取其 K=3 对应的标准文件，不再误用 K=5 状态。
- `plot_A_time_evolution`：不再生成参数专属 `A_time_turing_seed_*`；读取标准初值，结果使用新的 `A_time_evolution_sharedseed_cache_*` 名称，旧缓存保留。

## 吸引概率分析的目标参考态

`analyze_attraction_probability_vs_c` 的插值终点需要是目标 sigma 下的斑图，不应将 sigma=100 的标准种子直接冒充目标态。

按严格读取约定，该脚本只使用已有缓存中的有效参考态，或通过 `load_existing_pattern_reference` 读取匹配的 `evolution_*` 文件；不自动训练、不延拓。匹配依据为 MAT 中保存的 cfg（N、K、alpha、beta、sigma、noise 和网络矩阵），而不是精度有限的文件名。找不到目标态就报错。

旧 `code/` 目录的普通扫描已是读取预生成种子的模式，本次未改动旧目录。求解器 `solve_multiplex` 的动力学、噪声、收敛与早停逻辑均未修改。
