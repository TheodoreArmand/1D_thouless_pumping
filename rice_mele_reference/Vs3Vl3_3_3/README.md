# Vs3/Vl3 3/3 Rice-Mele reference data

这个目录只放当前 3/3 Thouless pump 案例会用到的 Rice-Mele 参考数据。

## Files

- `gap_adaptive_vs3vl3_maincpp_schedule.csv`: 当前 C++ 配置读取的相位 schedule，列为 `s,phi`，其中 `s=t/T`。
- `gap_adaptive_Vs3Vl3_convention_note.md`: Vs3/Vl3 参考数据的约定说明。
- `gap_adaptive_grid_Vs3Vl3.npz`: Vs3/Vl3 gap-adaptive 网格数据。
- `gap_adaptive_grid_Vs3Vl3_maincpp.npz`: 和当前 `main.cpp` 相位约定对齐后的 gap-adaptive 网格数据。
- `gapadaptive_grid_ref_T160pi.npz`: `T=160*pi` 使用过的 gap-adaptive 参考数据。
- `grid_rm_curves_Vs3Vl3.npz`: Vs3/Vl3 Rice-Mele 曲线数据。
- `figs/`: 从旧工程搬来的 3/3 参考图和已有动画。

## Source

这些文件来自：

- `../periodic/rice_mele_reference/`
- `../periodic/rice_mele_grid/`

这里只归档 3/3 案例。其他 `Vs/Vl` 扫描结果、旧脚本和 `__pycache__` 没有搬进来。
