# Why the two "Vs3/Vl3" gap-adaptive figures differ

*Recorded 2026-07-05.*

> **Update (later 2026-07-05):** the 2-panel `figs/gap_adaptive_Vs3Vl3_maincpp.png`
> was removed as redundant. The canonical half-depth figure is the 3-panel
> `figs/gap_adaptive_grid_Vs3Vl3_maincpp.png` (displacement + band-0 population
> P₀ + schedule), now generated in-code by
> `make_gapadaptive_progress_report.py` (`fig_schedule_p0`) so it stays
> reproducible and in sync with the T=160π progress report. Mentions of the
> 2-panel file below are historical.

**Question.** `figs/gap_adaptive_Vs3Vl3_maincpp.png` and
`figs/gap_adaptive_grid_Vs3Vl3.png` are both labelled "Vs3/Vl3" gap-adaptive
pumps, yet the curves look different. Aren't they the same configuration?

**Short answer: no — they are different lattices.** The label is the same
($V_s=V_l=3\,E_R$) but the two scripts evolve potentials of **different depth**,
and they also **plot against different x-axes**. Both effects change the lines.

---

## 1. Different potential depth (the main reason)

| | `gap_adaptive_grid.py` → `gap_adaptive_grid_Vs3Vl3.png` | `make_vs3vl3_maincpp_schedule.py` → `gap_adaptive_Vs3Vl3_maincpp.png` |
|---|---|---|
| potential | `V = -Vs·cos(4πx/a) - Vl·cos(2πx/a + τ)` | `V = -0.5·Vs·cos(4πx/a) - 0.5·Vl·cos(2πx/a + φ)` |
| cosine amplitude | `Vs = 3 E_R` (**full depth**) | `0.5·Vs = 1.5 E_R` (**half depth**) |
| convention | grid benchmark (`grid_smoothness_diagnosis`, `compute_tight_parameters.py`) | **main.cpp effective** potential (matched point-for-point to the real ECG binary trace) |
| `g_min` (verified) | **0.474 E_R** | **0.562 E_R** |

Same numbers `Vs=Vl=3`, but the main.cpp version uses **half the cosine
amplitude** (the "sin²-depth" convention). This is the units gotcha the project
notes flag: the effective potential the ECG binary integrates is half-depth.

> Counterintuitive but correct: the *shallower* (half-depth) lattice has the
> *larger* gap. At the dimerized crossing (φ = π/2, 3π/2, where Δ = 0) the gap is
> `g_min = 2|J1 − J2|`, and the hoppings grow as the lattice shallows, so
> `|J1 − J2|` — and hence `g_min` — is bigger for the half-depth lattice.

### Consequence — why the red (uniform) curves differ most

At the **same** cycle time `T = 160π`:

| convention | `g_min` | uniform-drive band-0 pop `P0` | uniform curve |
|---|---|---|---|
| full depth (grid) | 0.47 E_R | **0.72** (strongly non-adiabatic) | big Stückelberg oscillations, under-pumps |
| half depth (maincpp) | 0.56 E_R | **0.999** (nearly adiabatic) | reaches −a with only a tiny wiggle |

The larger gap of the half-depth lattice makes the uniform pump nearly adiabatic
at `T = 160π`, whereas the full-depth lattice is far from adiabatic at the same
`T`. (The gap-adaptive curve reaches `P0 ≈ 0.999` in **both**.)

---

## 2. Different x-axis and layout (why even the green curves differ)

- **Top-panel x-axis.** The maincpp figure plots the centroid vs **time** `t/T`
  (0 → 1); the grid figure plots it vs **pump phase** `τ/π` (0 → 2). For the
  gap-adaptive schedule the transport is spread ~uniformly *in time* but happens
  in *steps in phase*, so the **same** dynamics look **linear** (maincpp, vs
  `t/T`) versus a **staircase** (grid, vs `τ/π`).
- **Panel layout.** maincpp = 2 panels (centroid; `g_min(φ)` + schedule rate).
  grid = 3 panels (centroid; instantaneous band-0 population `P0(τ)`; schedule).

---

## 3. Which one should you trust?

For comparing against the **ECG run, use the half-depth main.cpp figure**
(`gap_adaptive_Vs3Vl3_maincpp.png`) — that is the potential the code actually
integrates. The full-depth grid figure is the benchmark-convention version and
will systematically differ (smaller gap, less adiabatic at a given `T`).

## 4. Apples-to-apples reproduction

- Half-depth (main.cpp) for any Vs/Vl:
  `python3 make_vs3vl3_maincpp_schedule.py <Vs> <Vl>`
- Full-depth (grid) is hard-wired to Vs3/Vl3 in `gap_adaptive_grid.py`; to compare
  the *same physics*, either halve the depths there or run the maincpp script at
  double depth. Do **not** compare the two figures directly as-is — the depth
  convention differs by a factor of 2.

## 5. Update (2026-07-05): grid figure regenerated in main.cpp convention

`gap_adaptive_grid.py` now takes a `--maincpp` flag (half depth) and an optional
`--Tpi=<N>`; the default stays full-depth (for the report). Running
`python3 gap_adaptive_grid.py --maincpp` writes
`figs/gap_adaptive_grid_Vs3Vl3_maincpp.png` (3-panel, half-depth, g_min 0.56).

**Key finding of the corrected figure:** at `T = 160π` in the real half-depth
convention, the uniform pump is **already adiabatic** (`P0 = 0.999`), so uniform
and gap-adaptive essentially overlap — the dramatic "gap-adaptive rescues the
pump" contrast in the old *full-depth* figure was a **convention artifact** of the
smaller (0.47 E_R) gap. In the correct convention gap-adaptive only visibly helps
at shorter cycle times, e.g. `T = 80π`: uniform `P0 = 0.96` vs gap-adaptive
`P0 = 0.997` (run `gap_adaptive_grid.py --maincpp --Tpi=80`).

> Caveat: the report `thouless_pumping_rice_mele_reference.html` §8.5 still shows
> the FULL-depth `gap_adaptive_grid_Vs3Vl3.png` and mislabels it "main.cpp
> convention". That figure/narrative is valid *in the full-depth grid convention*;
> the label should read "full-depth grid convention" to be precise.

## Files
- `gap_adaptive_grid.py` → `figs/gap_adaptive_grid_Vs3Vl3.png` (full depth, g_min 0.47)
- `gap_adaptive_grid.py --maincpp` → `figs/gap_adaptive_grid_Vs3Vl3_maincpp.png` (half depth, g_min 0.56)
- `make_vs3vl3_maincpp_schedule.py` → `figs/gap_adaptive_Vs3Vl3_maincpp.png` (half depth, 2-panel, g_min 0.56)
- Related: `../rice_mele_grid/README.md` (grid-convention Rice–Mele parameter scans)
