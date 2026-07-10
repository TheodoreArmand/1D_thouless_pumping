# R5 详细提案：并集选基（union + conditioning-guarded selection）K32 → K40/48

日期：2026-07-10（深夜）。状态：**纯提案，未改任何代码**。
上游文档：`vs3_n2_k32_success_roadmap_report.html` §6-R5 / §8.5；`vs3_n2_k32_failure_analysis_report.html` §4。

## Abstract (EN)

R5 enriches the N=2 K=32 path-pad basis by **union, never replacement** (the COMdiag swap made
both free and gauss worse), selecting 8–16 additional two-particle "parking spots" from a
structured candidate pool by a **greedy, conditioning-guarded fit against the actual grid
wavefunction snapshots** ψ_grid(s) through the first transition (from R4). New terms enter
parked (|u| floor ≪ u_on) so the R3 parking property is preserved; the initial state is re-fit
over the full K=40/48 set, a matched grid reference is regenerated, and the same Gate-1
short-window acceptance applies. A K-neutral alternative (R5c "pad recycling" via R7a
checkpoints) is specified for the case where cost matters more than span.

---

## 0. 触发条件（不要提前执行）

R5 只有在下面三件事都成立时才启动：

1. R1（frozen-A gauss）与 R2（A01-only、普通 basis）短窗都没有过 Gate-1；
2. R3（门控白化）重跑后 discarded 峰值已降但 δΔP 仍超标 —— 即"求解器已修好但还是不行"；
3. **R4 探针判定确有 span 缺口**：对 s∈[0.05,0.40] 的 grid 快照，现有 K32 基组的
   u-only 投影 infidelity 在某段 s 上 > 10⁻²（如果 R4 给出 ≤10⁻³，问题不在基组，
   R5 不会有用——回到求解器/参数化路线）。

R4 同时给出 R5 最重要的输入：**缺口在哪个 s、哪块位形空间**（残差密度 |r_m(x1,x2)|² 热图）。

### 0.1 Schedule 兼容性（2026-07-10 深夜切换后必须遵守）

深度约定已切换（见 CLAUDE.md）：所有新 ECG 与 grid run 读
`rice_mele_reference/Vs3Vl3_3_3/gap_adaptive_vs3vl3_full_depth_schedule.csv`
（生成器 `make_vs3vl3_full_depth_schedule.py`）。三条纪律：

1. R4 的 ψ 快照、R5 的短窗 run、匹配 grid 参考必须使用**同一个（新）schedule**——
   混用新旧 schedule 的对比全部作废。
2. 本提案正文引用的 s 数值（跃迁 s≈0.2–0.3、窗口 [0.18,0.32] 等）来自旧 maincpp
   schedule 的映射；落地前用新 schedule 的 g(φ) 重新标定 fine 窗口与 s_end。
3. **快照与验收判据一律改用 φ 定位**（schedule 不变量）：建议 R4/R5 快照取
   φ/π ∈ {0.30, 0.32, …, 0.60}（步长 0.02，16 点），Gate-1 判据保持"φ ≤ 0.52π"表述。

## 1. 设计原则（从 2026-07 的失败中总结）

- **P1 并集，不替换。** COMdiag 把 16 个单粒子 pads 换成 16 个 COM 对，free/gauss 双双变差、
  0.01-越界点提前到 s≈0.085。单粒子 pads 是第一次跃迁的承重结构，必须保留。
- **P2 停车哲学不变。** 新增项一律以"停车 reserve"身份进入（|u| 在 floor 值），
  波包到达前不移动（依赖 R3a/R3b 的门控；u 通道常开负责接收振幅）。
- **P3 选基必须带条件数守卫。** 盲目加 Gaussian 会制造近线性相关 → 复刻 rank-collapse。
  每个候选入选前检查 novelty（见 §4）。
- **P4 一次只动一个旋钮。** R5 短窗运行的其它一切（演化模式、rcond、fine 窗口、σ=1）
  锁定为 R1/R2/R3 阶段的胜者配置；R5 唯一的变化是 basis 文件与 K。
- **P5 静态 span ≠ 动态可达。** 选基通过后仍必须过短窗动力学验收；
  静态 fidelity 表只用来筛选候选，不用来宣布成功（COMdiag 的教训：静态 0.9994，动态更差）。

## 2. 位形空间几何（候选池的物理依据）

代码单位：长晶格周期 a=8，短晶格周期 d_s=4。初态两玻色子占相邻长胞：
(x₁,x₂) ≈ (0, 8)。一个周期泵浦 −8：终态 (−8, 0)。两次跃迁各移动一个短格点（−4）：

```
s≈0.20–0.30（第一次跃迁）：(0,8) → (−4,4)
s≈0.30–0.70（平台）：       停在 (−4,4)
s≈0.70–0.80（第二次跃迁）：(−4,4) → (−8,0)
```

跃迁中每个粒子在旧井/新井间离域，二粒子密度呈 4 瓣：
(0,8), (0,4), (−4,8), (−4,4)，外加势垒中点（−2、6 等）的过渡振幅。

现有 K32（ordinary path-pad）的覆盖：

| 块 | 项数 | 位形空间覆盖 |
|---|---|---|
| core（fit 乘积） | 16 | 初态 blob (≈0, ≈8) |
| 单粒子 pads（左动） | 8 | 横线 (−2/−4/−6/−8, ≈8)，pad 宽度 b/n |
| 单粒子 pads（右动） | 8 | 竖线 (≈0, 6/4/2/0)，pad 宽度 b/n |

**缺口（若 R4 证实）**：对角线上的"双动"格点 (−2,6),(−4,4),(−6,2),(−8,0)（COMdiag 已生成，
但当时是替换式引入）；以及"一格 + 半格"的交错组合 (−4,6),(−2,4)（第一跃迁）与
(−6,4),(−4,2),(−8,2),(−6,0)（第二跃迁）；以及既有一动配置缺失的联合宽度组合。

## 3. 候选池 P（结构化、分层）

宽度记号沿用 `make_vs3vl3_initial_basis_n2_comdiag.py`：b = N1 pad 宽（term 8），
n = N1 窄宽（term 9）；联合宽度四组合 {bb, bn, nb, nn}。除 P-E 外全部 A=0。
所有候选按未排序中心对去重（对称化后 (c₁,c₂) ≡ (c₂,c₁)）。

| 层 | 内容 | 数量 | 说明 |
|---|---|---|---|
| P-A | 现有 K32 全部 32 项 | 32 | **强制保留，不参与淘汰**（P1） |
| P-B | 双动对角 (−2,6),(−4,4),(−6,2),(−8,0) × {bb,bn,nb,nn} | 16 | 直接复用 COMdiag 生成器的 reserve 段 |
| P-C1 | 第一跃迁交错 (−4,6),(−2,4) × {bb,bn,nb,nn} | 8 | "一个到新井、一个在势垒"构型 |
| P-C2 | 第二跃迁交错 (−6,4),(−4,2),(−8,2),(−6,0) × {bb,bn,nb,nn} | 16 | 全周期时才需要；Gate-1 选基可先排除 |
| P-C3 | 既有一动配置的缺失联合宽度：(−4,8),(0,4),(−2,8),(0,6) × 补齐 {bb,bn,nb,nn} 中现无者 | ≤8–16 | 现有 pads 的 anchor 宽度是固定的，只覆盖 2/4 组合 |
| P-E（可选，gauss 专用） | P-B 各中心 × A=0.05·[[1,−1],[−1,1]] | 16 | 相对坐标关联 seed；**必须修好 u 补偿**（§5），且 grid loader 需支持 A（§6） |

池规模：Gate-1 版（A+B+C1+C3）≈ 64–72；全周期版（+C2）≈ 80–88；含 P-E ≈ 96–104。

## 4. 选基算法（greedy + novelty 守卫；~150 行 Python）

**目标**：从池中选出 K_target ∈ {40, 48} 项（含强制的 32），最小化
"最差快照的 u-only 投影 infidelity"，同时不引入近线性相关。

输入：
- R4 的 ψ_grid(s_m) 快照，s_m ∈ {0, 0.05, 0.075, …, 0.40}（512² 网格即可）；
  全周期版再加 s_m ∈ {0.55, 0.60, …, 0.95}；
- 候选池 CSV（basis 格式：每项 u; A(N²); B(N²); R(N); name）。

步骤（全部数值都在快照网格上做，复用 `n2_grid_reference.py` 第 91 行的
对称化 primitive 定义 (g(x₁)h(x₂)+h(x₁)g(x₂))/√2）：

```
1  对每个候选 i：在网格上装配对称化 primitive，归一化 ĝ_i。
2  I ← P-A 的 32 项。
3  对每个快照 m：解 (S_II + 1e-10·I) c = b_{I,m}，b_{i,m}=⟨ĝ_i|ψ_m⟩；
   残差 r_m = ψ_m − Σ c_i ĝ_i；记录 infid_m = 1 − |⟨ψ̂_fit|ψ̂_m⟩|²。
4  while |I| < K_target 且 max_m infid_m > 1e-3：
4a   对每个 j ∉ I：novelty_j = 1 − ‖P_I ĝ_j‖²（P_I 用同一正则化解）
4b   若 novelty_j < 0.005 → 跳过（条件数守卫，P3）
4c   score_j = max_m |⟨ĝ_j|r_m⟩|² / (novelty_j + 1e-6)
4d   选 score 最大的 j*，I ← I ∪ {j*}，重算步骤 3
5  输出选基 CSV + 选基报告
```

备选：把 4a–4d 换成对快照加权矩阵 [⟨ĝ_i|ψ_m⟩] 的 pivoted-QR（失败报告 §4 的建议）——
结果通常一致；greedy 版胜在能对每一步给出可读的"为什么选它"。

**选基报告（必须产出，进 git）** `rice_mele_reference/Vs3Vl3_3_3/k48_selection_report.md`：
- 入选项表：中心、联合宽度、层号、入选轮次、novelty、贡献最大的快照 s；
- 每快照 infidelity：K32 原值 vs K_target 终值（这一列同时是 R4 结论的复核）；
- 归一化重叠矩阵 S 的谱：λ_min、λ_max（λ_min < 10⁻⁸ 时给出警告并建议收紧 novelty 阈值）；
- 每快照残差密度热图 |r_m(x₁,x₂)|² PNG（人工检查缺口是否真的被补上）。

停止条件二选一：max_m infid ≤ 10⁻³（提前停，K 可能 < 48——更好），或 K=48 打满
（此时若 infid 仍 > 10⁻²，宣布"池不够"，走 §9 matching-pursuit 分支）。

## 5. 新项的 u 初始化与初态重拟合（Aseed 事故的防重演清单）

失败报告 §4 记录的事故：Aseed 生成器改 R 后没有补偿指数常数项，远端 reserve 的
norm contribution 从 ~10⁻⁶ 掉到 10⁻¹⁶–10⁻¹⁹，seed 实际不存在。R5 的规程：

1. **全集重拟合**：定下 K_target 项后，对初态目标（与 K32 相同的对称化相邻胞乘积目标）
   重新做一次带正则化的最小二乘，得到全部 u（生成器现成流程，扩到 K48）。
2. **u 下限**：远端新项拟合值会是 ~0。设 floor：|u_k| ← max(|u_fit|, 1e-4·max_j|u_j|)
   ≈ 5×10⁻⁵（相位取 0）。这在 norm² 里只贡献 ~10⁻⁸（不伤初态保真度），
   又远高于数值尘埃、远低于 R3 的唤醒阈 u_on≈0.02–0.05（停车逻辑自洽）。
3. **P-E 项的 u 补偿**：给中心 m、宽度 B 的项加 A 时，保持"函数峰位与峰值"不变：
   重解 R 使包络均值不动之外，还必须把被改变的指数常数吸收进 u
   （逐项验收：加 A 前后每个 primitive 的 L² norm 之比 ∈ [0.5, 2]）。
4. **三条静态验收**（smoke 前置）：
   - 初态保真度 |⟨ψ_{K48}(0)|ψ_{K32}(0)⟩|² ≥ 0.99999；
   - 初能量差 |E_{K48}(0) − E_{K32}(0)| < 10⁻⁸（free 与 gauss 各测一次）；
   - rcond=1e-4 下初始 solve rank 与 sv 谱与 K32 相比无病态跳变（smoke run 自动给出）。

## 6. 匹配 grid 参考（对照必须 apples-to-apples）

初态重拟合后 ψ(0) 有微小变化 ⇒ 照 COMdiag 的先例生成匹配参考：

- `out/vs3_n2_dt0p01_T160pi_K32_sweep/grid_ref_k48union/{free,gauss_sigma1p000}`，
  每 case ~2 h @32 cores（复用 `ecg1d_vs3_n2_grid_comdiag_20260710.sbatch` 模板）；
- live report 变体照 COMdiag 的改法：主参考指向 k48union，普通 grid 留作 audit 列；
- **前置修复（若 P-E 入选则必须，否则也建议顺手做）**：当前 grid loader 只读
  u、B 对角、R，忽略 A 与非对角 B（失败报告 §3 的 A-seed 例外正源于此）。
  改为直接按完整 primitive 公式 u·exp(−xᵀAx−(x−R)ᵀB(x−R)) 在网格上装配（~20 行 numpy），
  用解析 norm 做单元测试。修完后 A-seed 类比较从此不再有例外条款。

## 7. 运行矩阵

模式锁定 = R1/R2/R3 的胜者（例：frozen-A + rcond 1e-5 + R3a 门控；以实际结果为准）。

| 步骤 | 内容 | 成本 | 通过才继续 |
|---|---|---|---|
| R5.0 | 池生成 + 选基 + 报告 + K48 CSV | ~0.5 天（人）+ 分钟级（机） | infid ≤ 1e-3（或明确宣布池不够） |
| R5.1 | smoke ×4（free/gauss × rcond 两档；`ECG_SMOKE_TOTAL_TIME` 钩子） | 分钟级 | §5 三条静态验收 |
| R5.1b | 匹配 grid 重跑 ×2 | 2 h ×2（并行） | — |
| R5.2 | Gate-1 短窗 ×4：{gauss, free} × {rcond 胜者, 次者}，s_end=0.35，fine 0.1×[0.18,0.32] | 每发 8–18 h（K48 ≈ 2.25× pair 表 + param_dim 240/288 的装配；估 0.4–0.7 s/step × 4.6–8.1 万步） | Gate-1 全套判据（§8） |
| R5.3 | 胜者全周期（含 R7 工程件） | 1.5–3 天 wall | Gate-2 判据 |

driver 改动：给 Gate-1 短窗 driver 加一个 basis 位置参数
`ordinary | comdiag | k48union`，映射到 CSV 路径并同步设 `cfg.K`（basis 行数必须等于 K，
loader 会核对）；`config_appendix` 记录 `basis_variant=K48_union selection_report=<path>
pool_tiers=A+B+C1+C3 u_floor=1e-4max`。sbatch 照 Gate-1 模板，`--time=24:00:00`。

## 8. 验收判据（在 Gate-1 全套之上新增三条）

Gate-1 原判据：|ΔP_ECG−ΔP_grid|(φ≤0.52π) < 0.01；discarded 峰值 < 0.1；
λ_min[Re(A+B)] > 10⁻²；对 rcond 的轨迹敏感性小。新增：

1. **停车验收**（同 R3）：所有 |u_k| < u_on 的 reserve 在被占据前 ‖R_k(t)−R_k(0)‖ < 1。
2. **唤醒次序**：从 snapshots.csv 离线画每项 |u_k|(s)——pads 应按路径顺序渐次占据
   （核 → 一动 pads → 对角 pads）。乱序唤醒 = 动力学在走别的路径，要人工看。
3. **free 对照不得变差**：K48 free 短窗的 δΔP 不劣于 K32 free 同 rcond
   （COMdiag 的一票否决教训——伤害 free 的基组直接淘汰）。

## 9. 失败分支

- **池不够**（K=48 打满仍 infid > 10⁻²）：matching pursuit 发现缺失形状——
  对最差快照的残差 r_m 拟合一个自由 Gaussian（中心、宽度、A 全放开，scipy 局部优化，
  多起点），把发现的形状加回池里重选。若连自由 Gaussian 都压不下残差
  （对这个光滑问题不太可能），才谈"根本性表示极限"。
- **静态过、动态仍败**：可达性问题而非 span 问题 → 回 R3 深挖
  （截断/白化细节、保留集抖动），或走 R6（schedule 减速降低对流形的需求），
  或上 R5c。
- **条件数守卫频繁触发 / 初始谱病态**：收紧 novelty 阈值到 0.01，降 K_target 到 40。

## R5c 备选：pad 循环（K 不变，"adiabatic basis transport"）

若 K48 成本不可接受或条件数不可控，用 R7a 的 checkpoint/restart 做 K=32 的
"停车位循环"：在 checkpoint（例 s=0.40、0.55）上，离线 Python 识别波包已经过、
|u_k| < 10⁻⁴ 且持续 Δs > 0.05 的**尾部**项，把它们重新停到波前**前方**的路径位置
（第二跃迁的 P-B/P-C2 格点），u 置回 floor，重拟合当前快照后 restart。
纯离线脚本（读 snapshots.csv → 写新 basis CSV → restart），引擎零改动，
每次循环成本 ≈ 一次重启。哲学与 path-pad 一致：**停车位跟着泵浦走，而不是无限加车位。**

## 10. 依赖与工时汇总

依赖：R4 的 ψ 快照（硬依赖）；R1/R2 的胜者模式（硬依赖）；R3 门控（强烈建议，
否则新 reserve 的停车只由 rcond 副作用保证）；R7a restart（仅 R5c 需要）。

| 项 | 工时 |
|---|---|
| 池生成器 + 选基脚本 + 报告 | 0.5 天 |
| grid loader A-fix + 单测 | 0.5 天（可选，P-E 入选则必做） |
| driver/sbatch 改动 | 1–2 h |
| smoke + 匹配 grid | ~0.5 天（多为等待） |
| 短窗 ×4 | 1 天 wall（并行） |
| **合计（到 Gate-1 判定）** | **≈ 2–3 天** |

## 11. 文件清单（拟新增/修改）

- 新增 `rice_mele_reference/Vs3Vl3_3_3/make_vs3vl3_basis_n2_union_pool.py`（池生成）
- 新增 `rice_mele_reference/Vs3Vl3_3_3/select_basis_greedy.py`（选基 + 报告）
- 新增 `initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K48_union.csv`（+K40 变体）
- 新增 `rice_mele_reference/Vs3Vl3_3_3/k48_selection_report.md`（选基报告，进 git）
- 修改 Gate-1 短窗 driver（basis 参数 + cfg.K）与对应 sbatch
- 修改 grid loader（支持 A 与非对角 B，见 §6）
- （R5c 时）新增 `rice_mele_reference/Vs3Vl3_3_3/recycle_pads_from_snapshot.py`
