# Lohse (2016) Thouless pump 经验参考

对 `reference/lohes_nature_physics.pdf`（Lohse et al., Nature Physics 12, 350 (2016)）
的定量整理：连续超晶格哈密顿量、Rice–Mele 映射（J1/J2/Δ）、能隙、绝热时间、
均匀 vs gap-adaptive 驱动的实时模拟。**看 `index.html`。**

## Files

- `index.html`: 主报告（哈密顿量、J1/J2/Δ、能隙、绝热时间、模拟曲线）。
- `make_lohse_reference.py`: 一条命令重算全部数字与图（约 10–15 分钟）：
  `python3 make_lohse_reference.py`
- `lohse_reference_data.npz` / `lohse_reference_summary.json`: 计算结果。
- `gap_adaptive_lohse_schedule.csv`: 论文约定的 gap-adaptive 相位 schedule（`s,phi`）。
- `gap_adaptive_lohse_maincpp_schedule.csv`: 换算到 main.cpp 约定
  （`phi^m(s) = 3pi/2 - phi(s)`，递减，泵 +a/cycle —— 与 3/3 案例的 -a 相反）。
- `make_lohse_initial_basis.py`: 生成 ECG 初态
  `initial_state/Vs10Ers_Vl5Ers/initial_lohse_N1_K16.csv`
  （band-0 Wannier 拟合，fidelity 0.9989，path-pad 铺向 +x）。
- 对应的 C++ 配置：`pumpconfig/lohse_10_5.cpp`（T=500pi code units，已 smoke test：
  初始能量与本目录 Python 预测吻合到 1e-7）。
- `figs/`: 全部由脚本生成。

## System

Vs = 10 E_rs, Vl = 20 E_rl = 5 E_rs, dl = 2 ds（87Rb, λl = 1534 nm）。
单位 ds = 1, E_rs = 1。算法与 `../Vs3Vl3_3_3/` 参考数据同源
（平面波 Bloch、投影位置算符 Wannier、分裂步 TDSE、dφ/dt ∝ g² schedule）。
