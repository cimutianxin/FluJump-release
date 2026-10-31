# figdata/fig3_anchor — 方向的生物学锚定（Results 3.2）

- `panela_binding_anchor.json`：binding 监督逐层 AUC + jump/jump_human 零样本 AUC（双口径）。
- `panelb_ism_deltalogit.csv`：ISM 饱和突变逐突变 Δlogit（59,186 行；已删 mut_seq 整序列列）。
  口径：论文数值用 **bg_idx < 72**（72 个 H1+H3 train 背景；bg_idx 72–81 为 10 条 H5 背景，仅 G228S/Q226L 有）
  与 **dlogit_label_is_jump_human** 列（caption 的 5 替换与 0.80/1.01 均出自该列）。
  control 行的 site 列为空。
由 figs/fig3_anchor/prepare_figdata.py 从 validation_exp 已同步输出生成（本地重打包，无新实验）。
