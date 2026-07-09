# Time Summary

Generated: 2026-07-10T00:01:35

Rows scanned: 84

Full per-run data: `tmp/time_summary_all.csv`.

## Family Aggregate

| family | runs | done | progress | mean % | mean s/step | sum wall h | max ETA h |
| --- | --- | --- | --- | --- | --- | --- | --- |
| bench_64core_10step | 4 | 4 | 0 | 100.0 | 8.434 | 0.094 |  |
| bench_legacy_prb_3_3_300step_8thread | 1 | 1 | 0 | 100.0 | 0.045 | 0.004 |  |
| bench_vs3_n2_64core_10step | 2 | 2 | 0 | 100.0 | 1.893 | 0.011 |  |
| bench_vs3_n2_parallel_k32 | 10 | 10 | 0 | 100.0 | 0.595 | 0.022 |  |
| bench_vs3_n2_parallel_smoke | 2 | 2 | 0 | 100.0 | 0.226 | 0.001 |  |
| dtscan_lohse | 9 | 9 | 0 | 100.0 | 1.971 | 1.581 |  |
| vs3_n2_dt0p01_T160pi_K32_sweep | 8 | 0 | 8 | 37.45 | 0.466 | 49.64 | 19.85 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20 | 8 | 0 | 8 | 24.06 | 0.489 | 56.65 | 36.18 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em5 | 2 | 0 | 2 | 0.311 | 0.462 | 0.183 | 35.40 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em6 | 2 | 0 | 2 | 0.363 | 0.472 | 0.189 | 44.55 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5 | 8 | 0 | 8 | 27.82 | 0.421 | 22.41 | 8.827 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20 | 8 | 0 | 8 | 2.101 | 0.571 | 24.65 | 226.2 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat | 8 | 0 | 8 | 0.693 | 0.583 | 8.196 | 207.3 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged | 8 | 0 | 8 | 38.61 | 0.200 | 44.10 | 11.44 |

## Completed Benchmarks

| run | steps | dim | dt | wall h | s/step |
| --- | --- | --- | --- | --- | --- |
| bench_64core_10step/n1_3_3 | 10 | 48 | 0.002 | 0.001 | 0.315 |
| bench_64core_10step/n1_lohse | 10 | 48 | 0.002 | 0.001 | 0.311 |
| bench_64core_10step/n2_lohse_free | 10 | 120 | 0.002 | 0.024 | 8.539 |
| bench_64core_10step/n2_lohse_gauss | 10 | 192 | 0.002 | 0.068 | 24.57 |
| bench_legacy_prb_3_3_300step_8thread | 300 | 48 | 0.002 | 0.004 | 0.045 |
| bench_vs3_n2_64core_10step/free | 10 | 120 | 0.002 | 0.002 | 0.591 |
| bench_vs3_n2_64core_10step/gauss | 10 | 192 | 0.002 | 0.009 | 3.194 |
| bench_vs3_n2_parallel_k32/local2/free | 2 | 160 | 0.002 | 0.000 | 0.646 |
| bench_vs3_n2_parallel_k32/local2/gauss | 2 | 256 | 0.002 | 0.001 | 1.486 |
| bench_vs3_n2_parallel_k32/schemeC_omp16/free | 20 | 160 | 0.002 | 0.001 | 0.248 |
| bench_vs3_n2_parallel_k32/schemeC_omp16/gauss | 20 | 256 | 0.002 | 0.003 | 0.525 |
| bench_vs3_n2_parallel_k32/schemeC_omp32/free | 20 | 160 | 0.002 | 0.002 | 0.299 |
| bench_vs3_n2_parallel_k32/schemeC_omp32/gauss | 20 | 256 | 0.002 | 0.003 | 0.570 |
| bench_vs3_n2_parallel_k32/schemeC_omp64/free | 20 | 160 | 0.002 | 0.002 | 0.381 |
| bench_vs3_n2_parallel_k32/schemeC_omp64/gauss | 20 | 256 | 0.002 | 0.004 | 0.691 |
| bench_vs3_n2_parallel_k32/schemeC_omp96/free | 20 | 160 | 0.002 | 0.002 | 0.447 |
| bench_vs3_n2_parallel_k32/schemeC_omp96/gauss | 20 | 256 | 0.002 | 0.004 | 0.659 |
| bench_vs3_n2_parallel_smoke/free | 10 | 120 | 0.002 | 0.000 | 0.142 |
| bench_vs3_n2_parallel_smoke/gauss | 10 | 192 | 0.002 | 0.001 | 0.310 |
| dtscan_lohse/dt0p002 | 750 | 48 | 0.002 | 0.415 | 1.993 |
| dtscan_lohse/dt0p003 | 501 | 48 | 0.0030000000000000001 | 0.270 | 1.942 |
| dtscan_lohse/dt0p004 | 375 | 48 | 0.0040000000000000001 | 0.200 | 1.922 |
| dtscan_lohse/dt0p005 | 301 | 48 | 0.0050000000000000001 | 0.159 | 1.907 |
| dtscan_lohse/dt0p006 | 250 | 48 | 0.0060000000000000001 | 0.137 | 1.971 |
| dtscan_lohse/dt0p007 | 215 | 48 | 0.0070000000000000001 | 0.121 | 2.033 |
| dtscan_lohse/dt0p008 | 188 | 48 | 0.0080000000000000002 | 0.098 | 1.885 |
| dtscan_lohse/dt0p009 | 167 | 48 | 0.0089999999999999993 | 0.090 | 1.935 |
| dtscan_lohse/dt0p010 | 150 | 48 | 0.01 | 0.090 | 2.150 |

## K32 Long-Run Progress

| run | step | est | % | elapsed h | s/step | ETA h | dt |
| --- | --- | --- | --- | --- | --- | --- | --- |
| vs3_n2_dt0p01_T160pi_K32_sweep/free | 105975 | 140747 | 75.29 | 6.207 | 0.211 | 2.037 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma0p500 | 47775 | 140747 | 33.94 | 6.206 | 0.468 | 12.08 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma0p750 | 48175 | 140747 | 34.23 | 6.208 | 0.464 | 11.93 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma1p000 | 45900 | 140747 | 32.61 | 6.208 | 0.487 | 12.83 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma1p250 | 47175 | 140747 | 33.52 | 6.205 | 0.474 | 12.31 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma1p500 | 45675 | 140747 | 32.45 | 6.205 | 0.489 | 12.92 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma2p000 | 33500 | 140747 | 23.80 | 6.199 | 0.666 | 19.85 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep/gauss_sigma3p000 | 47525 | 140747 | 33.77 | 6.200 | 0.470 | 12.16 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/free | 121150 | 241277 | 50.21 | 7.091 | 0.211 | 7.031 | 0.010 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma0p500 | 50250 | 241277 | 20.83 | 7.091 | 0.508 | 26.96 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma0p750 | 54700 | 241277 | 22.67 | 7.085 | 0.466 | 24.17 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma1p000 | 56300 | 241277 | 23.33 | 7.086 | 0.453 | 23.28 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma1p250 | 40050 | 241277 | 16.60 | 7.074 | 0.636 | 35.54 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma1p500 | 39450 | 241277 | 16.35 | 7.072 | 0.645 | 36.18 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma2p000 | 50625 | 241277 | 20.98 | 7.073 | 0.503 | 26.64 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20/gauss_sigma3p000 | 51900 | 241277 | 21.51 | 7.073 | 0.491 | 25.81 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em5/free | 1000 | 241277 | 0.414 | 0.110 | 0.395 | 26.33 | 0.010 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em5/gauss_sigma1p000 | 500 | 241277 | 0.207 | 0.074 | 0.529 | 35.40 | 0.010 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em6/free | 1250 | 241277 | 0.518 | 0.096 | 0.278 | 18.51 | 0.010 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine20_rcond1em6/gauss_sigma1p000 | 500 | 241277 | 0.207 | 0.093 | 0.666 | 44.55 | 0.010 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/free | 41900 | 90481 | 46.31 | 2.444 | 0.210 | 2.833 | 0.010 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma0p500 | 22400 | 90481 | 24.76 | 2.857 | 0.459 | 8.682 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma0p750 | 22425 | 90481 | 24.78 | 2.856 | 0.458 | 8.667 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma1p000 | 23575 | 90481 | 26.06 | 2.855 | 0.436 | 8.103 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma1p250 | 23375 | 90481 | 25.83 | 2.855 | 0.440 | 8.196 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma1p500 | 23475 | 90481 | 25.94 | 2.857 | 0.438 | 8.155 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma2p000 | 22050 | 90481 | 24.37 | 2.844 | 0.464 | 8.827 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_fine5/gauss_sigma3p000 | 22175 | 90481 | 24.51 | 2.845 | 0.462 | 8.764 | 0.002 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/free | 29250 | 1005310 | 2.910 | 3.093 | 0.381 | 103.2 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma0p500 | 19000 | 1005310 | 1.890 | 3.090 | 0.585 | 160.4 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma0p750 | 28500 | 1005310 | 2.835 | 3.062 | 0.387 | 105.0 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma1p000 | 27250 | 1005310 | 2.711 | 3.066 | 0.405 | 110.0 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma1p250 | 20750 | 1005310 | 2.064 | 3.090 | 0.536 | 146.6 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma1p500 | 14250 | 1005310 | 1.417 | 3.086 | 0.780 | 214.6 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma2p000 | 16500 | 1005310 | 1.641 | 3.090 | 0.674 | 185.2 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20/gauss_sigma3p000 | 13500 | 1005310 | 1.343 | 3.079 | 0.821 | 226.2 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/free | 14250 | 1005310 | 1.417 | 1.072 | 0.271 | 74.54 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma0p500 | 7250 | 1005310 | 0.721 | 1.069 | 0.531 | 147.1 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma0p750 | 5000 | 1005310 | 0.497 | 1.036 | 0.746 | 207.3 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma1p000 | 6750 | 1005310 | 0.671 | 1.042 | 0.556 | 154.1 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma1p250 | 6500 | 1005310 | 0.647 | 1.019 | 0.565 | 156.6 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma1p500 | 5250 | 1005310 | 0.522 | 1.006 | 0.690 | 191.7 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma2p000 | 5500 | 1005310 | 0.547 | 0.977 | 0.640 | 177.6 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_full20_repeat/gauss_sigma3p000 | 5250 | 1005310 | 0.522 | 0.974 | 0.668 | 185.6 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/free | 212500 | 277470 | 76.58 | 5.520 | 0.094 | 1.688 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma0p500 | 92250 | 277470 | 33.25 | 5.520 | 0.215 | 11.08 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma0p750 | 90500 | 277470 | 32.62 | 5.511 | 0.219 | 11.39 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma1p000 | 91750 | 277470 | 33.07 | 5.506 | 0.216 | 11.15 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma1p250 | 95250 | 277470 | 34.33 | 5.513 | 0.208 | 10.55 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma1p500 | 92000 | 277470 | 33.16 | 5.513 | 0.216 | 11.11 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma2p000 | 90250 | 277470 | 32.53 | 5.515 | 0.220 | 11.44 | 0.001 |
| vs3_n2_dt0p01_T160pi_K32_sweep_staged/gauss_sigma3p000 | 92500 | 277470 | 33.34 | 5.506 | 0.214 | 11.01 | 0.001 |
