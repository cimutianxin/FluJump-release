# figdata/fig6_surveillance — 前瞻监测排序（Results 3.5）

- `panela_summary.csv`：三年×四方法前向 AUC（probe=cluster 含 CI95/perm p；VirHostPRED=cluster 单侧 p；naive/CNN 为 isolate+cluster 双口径行，10-06 起 Fig 6a 全方法统一用 cluster 行）。源：forward_multiyear.json / forward_2026_hardening.json / forward_2026.json / hostpred_baseline_forward.json。
- `panelb_2026_ranking.csv`：2026 前向 201 条逐株排序（probe jump_human 头 score/rank；10 阳性名次 1–9,11）。源：forward_2026_ranking.csv。
- `panelc_2024_cluster_ranking.csv`：2024 簇级排序 136 簇（score=簇内 mean logit，双标签；pos_ever=全期标签；first_jump=2024 首跳标记）。源：probe_scores_forward.csv ⋈ forward_2024_ranking.csv（prepare 脚本内断言复现 0.907/名次 3/10/13/20/28/36、(H5,70) 3 vs 55）。
- `paneld_2026_tiers.csv`：2026 逐株 jump_human logit/trainpct/tier（中+=8/10 TP、0/191 FP）。源：workflow_2026_tiers.csv。
- `train_logits_{2024,2025,2026}_{jump,jump_human}.csv`：三年 train logit 分布（10-03 服务器导出，见 notes/10-03-23_fig6_train_logits_export.md）。

本地重打包脚本：figs/fig6_surveillance/prepare_figdata.py（含 caption 数值断言）。
