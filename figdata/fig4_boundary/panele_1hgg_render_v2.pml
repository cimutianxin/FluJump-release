# fig4_boundary panel (e) v2：1HGG 共享核心驱动位点结构映射（视角/强调优化版）
# 与 v1（panele_1hgg_render.pml）同一数据源：核心列与着色值实时读自
# figdata/fig4_boundary/structure_colors.csv（残基列表不手抄）
# v2 变更：
#   1) 俯视图 zoom 到头部（以 RBS 为中心 +25A 缓冲，裁掉 stem），RBS 改蓝色大球更醒目
#   2) 两个视角都收紧边界（v1 留白偏多）
#   3) 侧面图保持 v1 已定稿视角
# 运行（仓库根目录）: /root/miniconda3/envs/pymol/bin/pymol -cq panele_1hgg_render_v2.pml

load validation_exp/reversal_sites/data/1HGG.pdb, hgg
remove hetatm
hide everything, hgg
show cartoon, hgg
set cartoon_fancy_helices, on
color grey85, hgg
bg_color white
set ray_opaque_background, on
set ray_shadows, on
set antialias, 2
set ray_trace_mode, 1

# ── RBS 参照（蓝：sticks + 大球，醒目优先；残基表同 v1，MARKER_SITES 口径）──
select RBS, (chain A+C+E) and (resi 134+135+136+137+138+186+190+193+194+221+222+224+225+226+227+228)
show sticks, RBS
show spheres, RBS and name CA
set stick_radius, 0.35, RBS
set sphere_scale, 1.0, RBS
color deepblue, RBS

# ── 共享核心列（白→红渐变，sticks + CA spheres；数据驱动读 csv）──
python
import csv
core_resi = {"A+C+E": [], "B+D+F": []}
with open("figdata/fig4_boundary/structure_colors.csv") as f:
    for row in csv.DictReader(f):
        if row["in_shared_core"] != "1":
            continue
        chains = "A+C+E" if row["chain"] == "HA1" else "B+D+F"
        core_resi[chains].append(row["h3_num"])
        cmd.alter("(chain %s) and resi %s" % (chains, row["h3_num"]),
                  "b=%.4f" % (-float(row["d_pos_mean"])))
for chains, resis in core_resi.items():
    cmd.select("CORE_" + chains.replace("+", "_"),
               "(chain %s) and (resi %s)" % (chains, "+".join(resis)))
cmd.create("CORE", "CORE_A_C_E or CORE_B_D_F")
python end
show sticks, CORE
show spheres, CORE and name CA
set stick_radius, 0.35, CORE
set sphere_scale, 0.6, CORE
spectrum b, white red, CORE

# ── 视角 1：侧面（同 v1 视角，zoom 收紧 5→4）──
orient hgg
turn y, 90
turn x, -90
zoom hgg, 4
png figdata/fig4_boundary/assets/panele_1hgg_side_v2.png, width=2400, height=2400, dpi=300, ray=1

# ── 视角 2：面向 RBS 俯视（zoom 到头部：RBS + 25A 缓冲，stem 出画）──
turn x, 90
zoom RBS, 25
png figdata/fig4_boundary/assets/panele_1hgg_top_v2.png, width=2400, height=2400, dpi=300, ray=1

quit
