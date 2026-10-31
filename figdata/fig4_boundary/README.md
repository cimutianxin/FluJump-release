# figdata/fig4_boundary — 失效边界 = H7/H10 亚支（Results 3.3，Figure 4）

- `panelb_geometry_points.csv`（既有，10-04 服务器导出）：11,060 株逐株 PC1/PC2/proj_w
  （L28、jump_human、raw logit 口径；meta 的 dev_note：旧 geometry_analysis.json 有
  ≤0.1 标签位移偏差，以本 csv 为准）。
- `panelb_unit_transform.json`（本脚本派生）：raw logit → 单位 w 投影（rel. train 质心）
  的线性映射，由 meta 的 8 个校验均值解出（max|残差|<2e-6）。caption 的 −5.2/−3.3 与
  "negative half-space" 以变换后/原始 logit 两个口径分别核对（H7 负半空间比例 0.994）。
- `panelc_profiles.csv`：H10/H4 × 双标签 × 31 层逐层 isolate/cluster AUC
  （方向检验；符号检验口径 = jump 的 cluster_auc，已复算核对 **20/31、6/31**——
  10-06 起为 09-19 词边界修复后标签口径，修复前为 24/31、3/31）。
- `panelc_preregistered.csv`：预注册层结果（jump 17/28；jump_human 13/28；含 bootstrap
  CI 与 P(AUC>0.5)）。H4 无 jump_human 阳性簇，其 jump_human cluster AUC 未定义（全 null）。
- `panelc_signtest_rank.json`：符号检验计数/p（10-06 修复后口径：H10 20/31 p=0.150、
  H4 6/31 p=8.8e-4）、人源 H10N8 簇（H10_197，江西 2013）
  rank 313→9/321（jump_human L13 翻转）、H10/H4 评估规模。
- `paneld_columns.csv`：逐位点驱动散点（jump L28 格，与结构图同格；region/共享核心标记
  来自 structure_colors.csv，d_pos 与 g3 输出逐列一致）。**10-06 起为 S5 词边界修复后
  口径**（H10 阳性簇 11）：仅 d_pos_H10 列随阳性集漂移（max |Δ| 11.9），
  d_pos_H7/region/core 标记不变。
- `paneld_overlap_summary.csv`：6 个 layer×label 格的 ρ(H7,H10)/p/top-20 Jaccard 及
  H4/H5 对照（H4 无 jump_human 阳性 → 对照仅 jump 三格）。修复后口径：ρ(H7,H10)
  0.631–0.718（旧 0.632–0.718）、Jaccard 0.176–0.333（旧上限 0.379）、H4 对照
  0.368–0.420；jump_human 三格标签未变、数值与修复前完全一致。
- `panele_summary.json` + `structure_colors.csv` + `assets/panele_1hgg_*_v2.png`
  （10-06 重渲染）：共享核心 17 列（成员不变，0 落 RBS）、JS 百分位 52.8–85.4
  （区间不变；jump L13 格 74.0→78.7，其共享集 11→10 列、列 669 跌出 H10 top-20）、
  亚型身份对照探针 ρ 0.11–0.15（不变）。

panelb/paneld/panele 由 figs/fig4_boundary/prepare_figdata.py 生成（本地重打包 + 断言，
无新实验）；panelc 与本次 paneld/panele 口径刷新在服务器侧完成（S5/S8，日志
notes/10-06-04_h10h4_label_audit.md、notes/10-06-*_reversal_sites_relabel.md）。
