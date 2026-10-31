# figdata/ — 绘图数据同步目录

服务器端实验**不直接出图**，只把绘图所需数据保存到本目录，经 git 同步到本地后统一绘图。

## 规则

- 每张图 / 每组图一个子文件夹，命名与本地 `figs/` 下的子文件夹**一一对应**：
  `figdata/<图名>/`（本仓库，服务器写） ↔ `figs/<图名>/`（本地，绘图）。
- 数据文件**只保存文本格式**（csv / tsv / json 等）。仓库 `.gitignore` 忽略 `.npy` / `.npz` / `.pt` 等二进制文件，保存为这类格式将无法通过 git 同步。
- **例外（10-04 起）**：结构渲染等**位图资产**（PNG）可放 `figdata/<图名>/assets/` 并正常 commit——此类是渲染产物而非绘图数据，本地 `figs/<图名>/` 组 EPS 时直接引用。
- 数据文件应自洽：包含绘图全部所需数值与列名说明，不依赖服务器本地路径。
- 写入后 `git add figdata && git commit && git push`；本地 `git pull` 后在 `figs/<图名>/` 中绘图。

## 现有子目录

| 目录 | 内容 | 来源脚本 |
|------|------|----------|
| `fig2_transfer/` | MLP head 迁移对比数据 | （历史） |
| `fig4_boundary/` | panel (b) 几何散点逐点数据（panelb_geometry_points.csv/meta.json）；panel (e) 1HGG 结构映射（structure_colors.csv + panele_1hgg_render.pml） | `ESM_clf/jump_exp/export_geometry_points.py`、`validation_exp/reversal_sites/export_structure_colors.py` |
| `fig5_depth/` | 150M cos(w_best, w_last)（target-val best 层 vs L30，主拟合 + 5 seed） | `validation_exp/depth_reversal/export_cos_w_150m.py` |
| `fig6_surveillance/` | 2024/2025/2026 三协议 train logit 分布（6 个 CSV） | `validation_exp/temporal_validation/run_multiyear_forward.py`、`run_2026_workflow_eval.py` |

注：`fig4_boundary/assets/` 下的渲染 PNG 按上方例外规则正常 commit（v1: panele_1hgg_{side,top}.png；v2: panele_1hgg_{side,top}_v2.png）。
