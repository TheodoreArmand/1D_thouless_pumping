# Vs3/Vl3 3/3 Rice-Mele reference data

这个目录只放当前 3/3 Thouless pump 案例会用到的 Rice-Mele 参考数据。

## Files

- `gap_adaptive_vs3vl3_full_depth_schedule.csv`: 当前 C++ 与 grid 配置共同读取的 canonical full-depth schedule，列为 `s,phi,g,inv_g2`，其中 `s=t/T`。
- `make_vs3vl3_full_depth_schedule.py`: 从 full-depth \(V_s=V_l=3E_R\) Bloch gap 可复现地生成上述 CSV。
- `gap_adaptive_vs3vl3_maincpp_schedule.csv`: 历史 half-depth-derived schedule；只用于复现旧 N=1 成功 case 和旧 N=2 数据，不再作为默认值。
- `gap_adaptive_Vs3Vl3_convention_note.md`: 当前 full-depth 协议、旧协议和数据 provenance 的约定说明。
- `gap_adaptive_grid_Vs3Vl3.npz`: full-depth Hamiltonian + full-depth schedule 的 N=1 网格数据。
- `gap_adaptive_grid_Vs3Vl3_maincpp.npz`: 历史 half-depth Hamiltonian 网格数据。
- `gapadaptive_grid_ref_T160pi.npz`: 历史 full-depth Hamiltonian + half-depth-derived schedule 参考数据。
- `grid_rm_curves_Vs3Vl3.npz`: Vs3/Vl3 Rice-Mele 曲线数据。
- `figs/`: 从旧工程搬来的 3/3 参考图和已有动画。

## Source

这些文件来自：

- `../periodic/rice_mele_reference/`
- `../periodic/rice_mele_grid/`

这里只归档 3/3 案例。其他 `Vs/Vl` 扫描结果、旧脚本和 `__pycache__` 没有搬进来。
