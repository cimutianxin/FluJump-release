# fig4_boundary panel (e)：1HGG 共享核心驱动位点结构映射
# 由 validation_exp/reversal_sites/export_structure_colors.py 生成（勿手改）
# 白→红渐变 = 推负强度（-d_pos_mean）；蓝 = RBS 参照（130-loop/190-helix/220-loop）
# 运行: /root/miniconda3/envs/pymol/bin/pymol -cq panele_1hgg_render.pml

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

# ── RBS 参照（蓝，sticks）──
select RBS, (chain A+C+E) and (resi 134+135+136+137+138+186+190+193+194+221+222+224+225+226+227+228)
show sticks, RBS
color deepblue, RBS
set stick_radius, 0.35, RBS

# ── 共享核心列（白→红渐变，sticks + CA spheres）──
alter (chain A+C+E) and resi 7, b=445.4301
alter (chain A+C+E) and resi 8, b=-220.9725
alter (chain A+C+E) and resi 9, b=325.1551
alter (chain A+C+E) and resi 20, b=180.4125
alter (chain A+C+E) and resi 21, b=-197.2007
alter (chain A+C+E) and resi 23, b=658.4346
alter (chain A+C+E) and resi 93, b=556.5577
alter (chain A+C+E) and resi 116, b=822.3290
alter (chain A+C+E) and resi 117, b=929.4715
alter (chain A+C+E) and resi 124, b=843.8376
alter (chain A+C+E) and resi 129, b=874.1136
alter (chain A+C+E) and resi 196, b=780.9492
alter (chain A+C+E) and resi 263, b=676.1779
alter (chain A+C+E) and resi 282, b=476.9202
alter (chain A+C+E) and resi 285, b=695.7790
alter (chain A+C+E) and resi 312, b=371.7668
alter (chain B+D+F) and resi 83, b=434.8267
select CORE, (chain A+C+E) and (resi 7+8+9+20+21+23+93+116+117+124+129+196+263+282+285+312) or (chain B+D+F) and (resi 83)
show sticks, CORE
show spheres, CORE and name CA
set stick_radius, 0.35, CORE
set sphere_scale, 0.6, CORE
spectrum b, white red, CORE

# ── 视角 1：侧面（HA 主轴竖直，头部/RBS 在上）──
orient hgg
turn y, 90
turn x, -90
zoom hgg, 5
png figdata/fig4_boundary/assets/panele_1hgg_side.png, width=2400, height=2400, dpi=300, ray=1

# ── 视角 2：面向 RBS 俯视（自头部沿主轴向下看；+90 转向头端）──
turn x, 90
zoom hgg, 6
png figdata/fig4_boundary/assets/panele_1hgg_top.png, width=2400, height=2400, dpi=300, ray=1

quit
