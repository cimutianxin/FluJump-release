# figdata/fig5_depth — 深度方向翻转 + 规模普适性（Results 3.4，Figure 5）

全部数值为 09-19 host 标签修复后口径；150M 逐层 AUC 取自
`validation_exp/depth_reversal/output/g1_layer_geometry.json`
（**非** `ESM_clf/jump_exp/output/layer_sweep.json`，后者为修复前口径）。

- `panela_per_layer_auc.csv`：三规模（150M/650M/3B）× 双标签逐层 AUC。
  `auc_h7`/`auc_h5` 为全 holdout raw-logit AUC；`auc_test_h13` 为同分布 H1+H3 test AUC；
  `rel_depth = layer / N_tr`，N_tr = 30/33/36（transformer 层数，L0=embedding 输出不入分母）。
- `panelb_scale_summary.csv`：cluster target-val 五 seed 选层与 test AUC（tab:scale 口径），
  末层全 holdout AUC、全 holdout 最优层/AUC。
- `panelc_gap_150m.csv`：150M 亚型内符号 gap g_H7(ℓ)=w·(μ⁺−μ⁻)（H7 口径）逐层，双标签，
  取自 g5_scale_geometry.json 的 150m_ref（与 g1 decomp 的 within+ − within− 逐位一致）。
- `paneld_aligned_stats_150m.csv`：150M 三统计量逐层（jump 标签）：gap_H7、gap_H5
  （由 g1 decomp 重建，复算 ρ(gap_H5,gap_H7)=−0.702/−0.773 与存储一致）、
  participation ratio、AA 身份解码 margin/acc（V1）。
- `summary_metrics.json`：caption 全部 headline 数值（ρ、变号边界、cos、PR/margin 最小值）。
- `h5_layer_scan.csv`（10-06，S7）：三规模 × 双标签 H5 逐层 holdout AUC 单列导出
  （列 scale, label, layer, rel_depth, h5_auc），自 `panela_per_layer_auc.csv` 提取，供 Fig 5 备选面板/补充材料。
- `cos_w_best_w_last.json`：10-03 导出的 150M 真末层（L30）口径 cos（本目录既有文件，
  prepare 只读不改）；650M/3B 的 cos 取自 g5 的 `cos_w_bestlayer_w_last`。

由 figs/fig5_depth/prepare_figdata.py 生成（本地重打包 + 断言，无新实验）。
