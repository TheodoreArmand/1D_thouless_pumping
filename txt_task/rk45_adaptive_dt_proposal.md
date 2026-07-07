# 自适应时间步进提案：RK45、以及比 RK45 更好的方案

*2026-07-07。纯提案，未改任何代码。针对 `lohse_10_5` 配置（V_s^m=5, V_l^m=2.5 E_rs，T=500π，
目前 fixed RK4 dt=2e-3 ≈ 785k 步），但结论对 legacy_prb_3_3 同样成立。*

---

## 0. TL;DR

| 方案 | 改动量 | 预期加速 | 风险 | 建议 |
|---|---|---|---|---|
| (0) dt scan 选定更大的 fixed dt | 1 行（`cfg.dt`） | 视 scan 结果，可能 2–5× | 无（scan 定量验证） | **先做**（正在跑，`out/dtscan_lohse/`） |
| (1) RK45 自适应 | ~100 行，3 个文件 | 1–3×（相对 scan 最优 fixed dt 可能只有 1× 出头） | 低 | 可选，价值主要是安全性 |
| (2) "dt 跟随 gap" 的手工 schedule | — | ~0 | 概念错位 | **不建议**（§3 解释为什么） |
| (3) 全局相位规范预处理 | ~20 行 + 诊断 | 未知，可能数倍 | 低 | 先做 1 小时诊断实验 |
| (4) u-split 指数积分器 | 新模块 ~200 行 | **10–100×**（期望） | 中（需系统验证） | 真正值得投入的方案 |

---

## 1. 问题framing：wall time 的三个因子

```
wall time = (T / dt) × (每步秒数)
```

- **T**（物理周期）已经被 gap-adaptive schedule 压到物理下限：φ̇ ∝ g² 把绝热要求从
  均匀驱动的 ~13000π（code units）压到 ~500π。这一项没有余量了。
- **每步秒数**由 K、param_dim、OMP 线程数决定（本机 8 核 ≈ 2.5 s/step；
  Vs3Vl3 的 summary 记录 0.59 s/step，应是更大的节点——上 SLURM 是免费的 4×）。
- **dt** 是剩下的数值旋钮，也是本提案的主题。当前 2e-3 是从 Vs3Vl3 继承的猜测值，
  没有做过收敛性标定。

---

## 2. 这个体系里，究竟是什么在限制 dt？

实时 ECG-TDVP 的 RK4 局部误差有三个来源，量级差别很大：

### 2a. 动力学相位旋转（主导，全周期近似恒定）

每个高斯的复系数 u_k 以 ~e^{−iE_eff t} 旋转，E_eff 的散布由基组决定
（窄高斯 B=2.6 的动能 ~B/2≈1.3，势能幅度 ±(V_s+V_l)≈±2.3 code units）。
smoke test 观测到 **|dz| ≈ 3.3**：参数矢量以 ~3/时间单位 的速率运动——
这与 gap 无关、与 schedule 无关、全周期近似恒定。
**这就是 dt 的真正约束。** RK4 对相位旋转的局部误差 ~ (ω·dt)⁵/120，
ω≈3–5，所以 dt 每放大 2×，误差放大 32×——能放大到哪里由正在跑的 scan 用数据回答。

### 2b. 驱动速率 φ̇（已经被 schedule 吸收，不构成约束）

从 `gap_adaptive_lohse_maincpp_schedule.csv` 实测（T=500π）：

```
|dphi/dt|:  min = 9.6e-5   mean = 4.0e-3   max = 0.155   (code units)
max/min = 1617  (= g_max²/g_min²)
每步相位增量 max|dphi|·dt:  dt=0.002 → 3.1e-4 rad;  dt=0.01 → 1.6e-3 rad
```

即使在 mid-cycle 驱动最快的地方、用最大的候选 dt=0.01，
哈密顿量每步也只变化 1.6 毫弧度——**驱动分辨率在任何地方都不是 dt 的约束**。

### 2c. TDVP 度规瞬态（偶发、局部）

C 矩阵接近奇异的瞬间（rank 跳变、path-pad 高斯被激活、两个高斯中心交叉）
会让 ż 短暂变大。fixed dt 对这类瞬态只能"硬扛"（表现为 raw_resid 尖峰）；
这是自适应步长真正有用武之地的地方。

---

## 3. 回应你的直觉："gap 小的时候用小 dt，gap 大的时候用大 dt"

这个直觉对 **φ（泵浦相位）** 是完全正确的——而且**已经实现了**：
φ̇ ∝ g² 的 schedule 正是"gap 小 → 相位走得慢，gap 大 → 相位走得快"，
它就住在 CSV 里，main.cpp 逐行插值执行。

但对 **dt（积分步长）** 这个直觉不成立，原因是一个抵消效应：

> gap 小的地方，schedule 已经把 φ̇ 压小了 1617 倍，哈密顿量随时间几乎不变；
> 而 dt 的真正约束（§2a 的相位旋转）与 gap 无关、到处一样。
> 所以"误差 vs 时间"曲线在整个周期里**近似是平的**——
> 这正是 g² schedule 的设计目标：让每单位时间的激发率恒定。

如果强行让 dt ∝ g²（跟着 gap 走）：crossing 处 Δφ/step ∝ g⁴，把算力浪费在
分辨一个几乎不动的哈密顿量上；mid-cycle 处 dt 放大 1617× 会立刻撞上相位旋转的
稳定性墙。clip 之后整条曲线退化回常数 dt——所以这个旋钮转不出东西来。

**你对 RK45 的怀疑因此是对的，但结论反了一半**：RK45 确实不会实现"dt 跟随 gap"
——因为它跟随的是*实测局部误差*，而这里的误差不由 gap 主导。
这不是 RK45 的缺陷，而是它会自动给出的正确答案。
附带的好处：跑一次 RK45，输出 dt(t) 的轨迹/直方图，
"哪里能大步走"就从猜测变成实验事实。

---

## 4. RK45 详细设计（如果要做）

### 4.1 现状盘点

- `src/realtime_tdvp.cpp` **已有** Dormand–Prince 5(4) 的 7-stage 实现
  （k1…k7，约 L257–292），但只支持 frozen H（整步用同一个哈密顿量），
  且只被 `realtime_tdvp_evolution`（非泵浦驱动器）使用。
- 时间依赖的 stage 求值模式已有样板：`realtime_tdvp_step_rk4_time_dependent`
  （k1 用 H(t)，k2/k3 用 H(t+dt/2)，k4 用 H(t+dt)）。
- `RealtimeEvolutionConfig` 已有 rk45_abs_tol / rel_tol / dt_min / dt_max 字段；
  `run_report` 的 summary.txt 已经预留了 `adaptive_*` 输出行（目前恒为 0）。
  ——也就是说，管线的三分之二早就修好了，缺的只是中间一段。

### 4.2 数学

Dormand–Prince 5(4)，stage 时刻 c = (0, 1/5, 3/10, 4/5, 8/9, 1, 1)：

```
k_i = f( t + c_i·dt,  z + dt·Σ_j a_ij·k_j ),   f(t,z) = −i C(z)⁺ g(z; H(t))
z5  = z + dt·Σ b_i·k_i          （5 阶解，用它推进）
z4  = z + dt·Σ b*_i·k_i         （4 阶嵌入解，只用于误差估计）
err = sqrt( mean_j |z5_j − z4_j|² / (atol + rtol·max(|z_j|,|z5_j|))² )
```

每个 stage 调一次 `compute_rhs_dz`（内含 C 组装 + SVD），
用 `terms_at(t + c_i·dt)` 取 stage 哈密顿量——与现有 RK4 时间依赖版完全同构。

### 4.3 步长控制（PI 控制器，工业标准）

```
接受 (err ≤ 1):  dt ← dt · clip( 0.9 · err^(−0.7/5) · err_prev^(0.4/5), 0.2, 5.0 )
拒绝 (err > 1):  dt ← dt / 2, 重试（连续拒绝 > max_reject 次则 abort）
额外拒绝条件（TDVP 特有）: 非有限参数 / min Re(B) ≤ 0 / raw_resid 突增 10×
dt 全程 clip 到 [dt_min, dt_max]
```

### 4.4 成本核算（诚实版）

- DP45 每步 7 stages，FSAL（k7 = 下一步的 k1）后有效 **6 次 rhs/step**，
  fixed RK4 是 4 次 → 每步贵 1.5×。
- 净加速 = (dt 平均增益) / 1.5。若 scan 显示 fixed dt=0.006 已经安全，
  RK45 平均 dt 大概也在 0.006 附近晃——净收益可能只有零点几倍。
- **真正的收益是**：(a) 度规瞬态处自动缩步，避免整条轨迹为最坏瞬间买单；
  (b) 换配置（更深晶格、更大 K、N>1）不用重新猜 dt；(c) dt(t) 本身是诊断数据。

### 4.5 代码接触面（约 100 行，3 个文件 + 配置）

| 文件 | 改动 |
|---|---|
| `src/realtime_tdvp.cpp` | 新函数 `realtime_tdvp_step_rk45_time_dependent`（~80 行：7 个 stage、嵌入误差、返回 err；照抄现有 frozen-H DP45 的系数表 + rk4_time_dependent 的 terms_at 模式） |
| `src/realtime_tdvp.hpp` | 1 个声明 + `RealtimeStepResult` 加一个 `err_norm` 字段 |
| `main.cpp` | 演化循环加 accept/reject（~20 行）：目前 `dt_step = min(cfg.dt, T−t)` 改为受控 dt；trace 已有 `used_dt` 列，无需改 I/O |
| `pumpconfig/pump_config.hpp` | 5 个字段：`use_rk45, rk45_abs_tol, rk45_rel_tol, dt_min, dt_max`（默认关闭 → legacy 行为逐位不变） |
| `src/run_report.cpp` | 把已预留的 `adaptive_*` 行接到真实计数上 |

### 4.6 验证计划

1. 关闭自适应（use_rk45=0）：逐位复现 legacy fixed-RK4 轨迹（回归保护）。
2. dt-scan 同款窗口（T_test=1.5）：RK45 的自动 dt 应落在 scan 判定的安全 dt 附近，
   末态参数与 dt=0.002 参考一致到容差水平。
3. 全周期对照（可先跑 T=160π 的 Vs3Vl3，41 h 有已知答案）：
   ΔP、norm/energy drift、solver 诊断都不劣于 fixed dt。

---

## 5. 比 RK45 更好的方案（你要的 "better proposal"）

### 5a. 全局相位规范预处理（半天，先做诊断）

TDVP 参数里最快的运动很可能是个**没有物理意义的自由度**：全局相位
（所有 u_k 同乘 e^{−iĒt}，Ē=⟨H⟩≈−0.727）。把它解析地吸收掉
（演化 ũ_k = u_k·e^{+iĒt}，或每步做规范投影），剩下的参数运动就只有
慢得多的相对运动，dt 的天花板直接抬高。

**先花 1 小时做诊断**（不改产品代码，写个一次性分析）：
把 smoke test 的 dz 向量投影到全局相位方向 iz 上，看 |dz_phase|/|dz| 占比。
- 占比 ~90%+ → 这个方案性价比爆表，先于 RK45 做；
- 占比小 → 说明 P_perp 度规已经处理了它，跳过。

### 5b. u-split 指数积分器（真正的 10–100×，~1 周）

观察：固定 (A,B,R) 时，u 的运动是**线性薛定谔方程**
（K 维 Galerkin 投影：`i S u̇ = H_eff u`，S/H_eff 是 K×K 的高斯重叠/哈密顿矩阵，
现有 `hamiltonian.cpp` 的机器就能组装），可以用 K×K 广义本征分解**精确积分**
——对 u 无条件稳定，dt 不再受相位旋转限制。剩下的 (B,R) 是慢的非线性运动，
用 RK4 粗步长即可。Strang 分裂：

```
半步精确演化 u (B,R 冻结) → 整步 RK4 演化 (B,R) (投影掉 u 方向) → 半步精确演化 u
```

- dt 的新约束 = (B,R) 的慢运动 + 分裂误差 O(dt³)，期望能放大 10–100×。
- K=16 的本征分解每步微秒级，成本可忽略；rhs 组装次数反而减少。
- **项目自己想过这条路**：summary.txt 里一直躺着 `u_split_trotter=0` 占位字段。
- 风险/工作量：与 P_perp 度规和 rcond 截断的一致性要重新推导；
  norm/energy 守恒判据重验；建议先在 Python 参考实现里原型验证
  （lohes_experience 的网格代码可当精确参照），再进 C++。

### 排序建议

```
(0) dt scan（进行中） → 定 fixed dt 基线，可能已经拿到 2–5×
(5a) 相位规范诊断 1 小时 → 决定是否做规范预处理
(1) RK45：如果还要继续换配置/加深晶格，作为安全网做上
(5b) u-split：想把 Lohse 全周期从"几天"变成"几小时"，做这个
```

---

## 6. 附：正在跑的 dt scan 是什么

`out/dtscan_lohse/dt0p00X/`：lohse_10_5 的前 T_test=1.5 code units
（恰好是 gap 最小、schedule 爬行、但相位旋转全速的窗口——对 §2a 最敏感），
dt = 0.002…0.010 共 9 组，其余参数逐位相同（含 OMP_NUM_THREADS=8）。
判据：以 dt=0.002 为参考，比较末态 basis_final.csv 的全参数矢量差
‖z(dt)−z(0.002)‖、E(t) 轨迹差、raw_resid，并检查差值随 dt 的 4 阶收敛标度。
结论（推荐 dt 及安全系数）在 scan 完成后另行给出。
