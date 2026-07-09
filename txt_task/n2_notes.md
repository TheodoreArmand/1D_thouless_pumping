# N=2 Lohse Pump Notes

Implemented targets:

- `ecg1d_2ninteraction`: two bosonic particles in the Lohse lattice, no pair interaction.
- `ecg1d_2gaussian`: same initial state with Gaussian pair repulsion.

Source layout:

- Shared copied pump loop: `nointeraction_src/pump_common.cpp`.
- Gaussian-specific setup: `2gaussian_src/gaussian_terms.cpp`.
- Configs: `pumpconfig/lohse_n2_free.cpp`, `pumpconfig/lohse_n2_gauss.cpp`.
- Initial state generator: `rice_mele_reference/lohes_experience/make_lohse_initial_basis_n2.py`.
- Grid reference: `rice_mele_reference/lohes_experience/n2_grid_reference.py`.
- SLURM launcher: `slurm/n2_validation.sbatch`.

Initial state:

- Adjacent-cell bosonic product of band-0 Wannier packets centered near x=2 and x=10.
- Current Vs3/Vl3 full-cycle test uses K=32: 16 dominant product Gaussians plus
  16 path-pad products. Each packet has four downstream centers with both the
  broad and narrow N=1 path-pad widths.
- Generator check: fidelity vs the full N=1 path-pad symmetrized product is
  about 0.9999928581 for the K=32 Vs3/Vl3 state.
- Generator energy prediction on the 2D grid:
  - free initial energy: about -1.43891916
  - Gaussian interaction energy for g=0.3 E_rs, sigma=1: about 2.21e-6

Cost warning:

The engine is N-general but the current assembly rebuilds a PairCache inside
every matrix-element/permutation/stage path. N=2, K=32 has a much larger
constant factor than the N=1 K=16 run, so the committed configs are short
validation windows (`total_time=12`) only. Full-cycle N=2 pumping should wait
for PairCache/permutation memoization.

Follow-up memoization scope:

- Cache `PairCache::build(i,j,p)` per basis snapshot and reuse it across
  overlap, Hamiltonian, trace, metric, and gradient kernels during one RHS.
- Reuse generated `PermutationSet` for fixed N.
- Keep the cache local to one basis state/RHS assembly to avoid invalidation
  problems during RK stages.
