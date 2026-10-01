# Plan 2 comparison: case33_18_25

Commit `b619864cfa81d54d9a1c637535565583e0a540b8`; case33bw load nodes [18, 25], mode 1, budget 7; partitions pp, pn, np, nn.

Settings: total 300 s wall clock per run; tau 0.005; eps 0.01; H: eps_A 0.15, share_A 0.25, eps_B 0.005, criterion auto, coverage global; Gurobi threads 1, seed 0, MISOCP cap 60 s, MIPGap 0.001; Gurobi 13.0.2.

Primary reference: AC scan (FR/MR in %, undecided AC cells excluded). SOCP columns are diagnostics: the OBBT-tightened model removes SOCP-only area, so outer-vs-SOCP MR > 0 is expected.

RCUT / RCUT2 (R + mainline cutting-plane construction instead of phase B): per network, SP scores the vertices of N_x (fixed x with its OBBT box and envelope rows) and the max-eta vertex yields a joint cut; a network stops after `patience` consecutive cuts each removing < `threshold` of vol(N_x). RCUT takes the stagnated/exact N_x as the network region (not a certified inner set); RCUT2 also runs the mainline inner certification (rays) and uses the certified hull. Inner I = (I_R ∪ sets) ∩ O_R; certified when vol(O_R∩box) - vol(I) <= eps·vol(I) (RB criterion). Runs: RCUT (threshold 1%, patience 3), RCUT-t0.5 (threshold 0.5%, patience 3), RCUT-t2 (threshold 2%, patience 3), RCUT2 (threshold 1%, patience 3), RCUT2-t0.5 (threshold 0.5%, patience 3), RCUT2-t2 (threshold 2%, patience 3).

## Runs

| method | workers | run | status | certified partitions | t_cert (s) | final gap | t gap<=10% | t gap<=5% | t gap<=2% | valid |
|---|---|---|---|---|---|---|---|---|---|---|
| H | 1 | 1 | certified | 4/4 | 75.8 | 0.0043 | 71.0 | 71.8 | 75.7 | yes |
| H | 1 | 2 | certified | 4/4 | 78.0 | 0.0043 | 73.2 | 74.1 | 77.9 | yes |
| H | 1 | 3 | certified | 4/4 | 77.9 | 0.0043 | 73.0 | 73.9 | 77.7 | yes |
| H | 16 | 1 | certified | 4/4 | 47.9 | 0.0043 | 21.0 | 21.9 | 25.7 | yes |
| H | 16 | 2 | certified | 4/4 | 73.5 | 0.0043 | 21.4 | 22.3 | 25.2 | yes |
| H | 16 | 3 | certified | 4/4 | 50.8 | 0.0042 | 21.3 | 22.1 | 26.6 | yes |
| R | 1 | 1 | certified | 4/4 | 74.0 | 0.0075 | 65.5 | 66.3 | 67.5 | yes |
| R | 1 | 2 | certified | 4/4 | 75.5 | 0.0075 | 66.8 | 67.6 | 68.9 | yes |
| R | 1 | 3 | certified | 4/4 | 76.6 | 0.0075 | 67.9 | 68.7 | 70.0 | yes |
| R | 16 | 1 | certified | 4/4 | 28.4 | 0.0075 | 19.9 | 20.7 | 21.9 | yes |
| R | 16 | 2 | certified | 4/4 | 29.0 | 0.0075 | 20.2 | 21.1 | 22.3 | yes |
| R | 16 | 3 | certified | 4/4 | 29.2 | 0.0075 | 20.2 | 21.1 | 22.4 | yes |
| RB | 16 | 1 | certified | 4/4 | 27.8 | 0.0076 | 21.6 | 22.5 | 24.8 | yes |
| RCUT | 16 | 1 | certified | 4/4 | 29.4 | 0.0070 | 21.8 | 22.7 | 26.4 | yes |
| RCUT-t0.5 | 16 | 1 | certified | 4/4 | 29.9 | 0.0057 | 21.8 | 22.7 | 26.8 | yes |
| RCUT-t2 | 16 | 1 | certified | 4/4 | 28.8 | 0.0056 | 21.6 | 22.5 | 24.6 | yes |
| RCUT2 | 16 | 1 | certified | 4/4 | 32.4 | 0.0074 | 21.5 | 22.3 | 29.4 | yes |
| RCUT2-t0.5 | 16 | 1 | certified | 4/4 | 33.0 | 0.0073 | 21.7 | 22.5 | 30.0 | yes |
| RCUT2-t2 | 16 | 1 | certified | 4/4 | 31.9 | 0.0078 | 21.4 | 22.3 | 28.9 | yes |

## Medians [min, max] over repeats

| method | workers | runs certified | t_cert (s) | final gap | t gap<=10% | t gap<=5% | t gap<=2% | inner-AC FR% | inner-AC MR% | outer-AC FR% | outer-AC MR% | inner-SOCP FR% | inner-SOCP MR% | outer-SOCP FR% | outer-SOCP MR% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H | 1 | 3/3 | 77.9 [75.8, 78] | 0.00425 | 73 [71, 73.2] | 73.9 [71.8, 74.1] | 77.7 [75.7, 77.9] | 0.809 | 0.0674 | 1.1 | 0 | 0 | 12.5 | 0.0727 | 12.2 |
| H | 16 | 3/3 | 50.8 [47.9, 73.5] | 0.00428 [0.00424, 0.00432] | 21.3 [21, 21.4] | 22.1 [21.9, 22.3] | 25.7 [25.2, 26.6] | 0.827 [0.809, 0.827] | 0.0796 [0.0735, 0.0796] | 1.12 [1.1, 1.12] | 0 | 0 | 12.5 [12.5, 12.5] | 0.0848 [0.0788, 0.0848] | 12.2 [12.2, 12.2] |
| R | 1 | 3/3 | 75.5 [74, 76.6] | 0.0075 | 66.8 [65.5, 67.9] | 67.6 [66.3, 68.7] | 68.9 [67.5, 70] | 0.75 | 0.276 | 1.08 | 0 | 0 | 12.7 | 0.0909 | 12.3 |
| R | 16 | 3/3 | 29 [28.4, 29.2] | 0.0075 | 20.2 [19.9, 20.2] | 21.1 [20.7, 21.1] | 22.3 [21.9, 22.4] | 0.75 | 0.276 | 1.08 | 0 | 0 | 12.7 | 0.0909 | 12.3 |
| RB | 16 | 1/1 | 27.8 | 0.00763 | 21.6 | 22.5 | 24.8 | 0.827 | 0.141 | 1.36 | 0 | 0 | 12.5 | 0.181 | 12.1 |
| RCUT | 16 | 1/1 | 29.4 | 0.00703 | 21.8 | 22.7 | 26.4 | 0.989 | 0.0245 | 1.43 | 0 | 0.109 | 12.4 | 0.254 | 12.1 |
| RCUT-t0.5 | 16 | 1/1 | 29.9 | 0.0057 | 21.8 | 22.7 | 26.8 | 0.935 | 0.0429 | 1.36 | 0 | 0.0546 | 12.4 | 0.181 | 12.1 |
| RCUT-t2 | 16 | 1/1 | 28.8 | 0.00559 | 21.6 | 22.5 | 24.6 | 1.11 | 0.0184 | 1.43 | 0 | 0.224 | 12.4 | 0.254 | 12.1 |
| RCUT2 | 16 | 1/1 | 32.4 | 0.00744 | 21.5 | 22.3 | 29.4 | 0.827 | 0.104 | 1.36 | 0 | 0 | 12.5 | 0.181 | 12.1 |
| RCUT2-t0.5 | 16 | 1/1 | 33 | 0.00735 | 21.7 | 22.5 | 30 | 0.827 | 0.104 | 1.36 | 0 | 0 | 12.5 | 0.181 | 12.1 |
| RCUT2-t2 | 16 | 1/1 | 31.9 | 0.00784 | 21.4 | 22.3 | 28.9 | 0.815 | 0.11 | 1.36 | 0 | 0 | 12.5 | 0.181 | 12.1 |

## Validity (final state, cells)

| method | workers | run | inner ∧ SOCP-infeasible | AC-feasible ∖ outer | SOCP-feasible ∖ outer (diagnostic) | undecided AC |
|---|---|---|---|---|---|---|
| H | 1 | 1 | 0 | 0 | 2300 | 0 |
| H | 1 | 2 | 0 | 0 | 2300 | 0 |
| H | 1 | 3 | 0 | 0 | 2300 | 0 |
| H | 16 | 1 | 0 | 0 | 2298 | 0 |
| H | 16 | 2 | 0 | 0 | 2298 | 0 |
| H | 16 | 3 | 0 | 0 | 2300 | 0 |
| R | 1 | 1 | 0 | 0 | 2306 | 0 |
| R | 1 | 2 | 0 | 0 | 2306 | 0 |
| R | 1 | 3 | 0 | 0 | 2306 | 0 |
| R | 16 | 1 | 0 | 0 | 2306 | 0 |
| R | 16 | 2 | 0 | 0 | 2306 | 0 |
| R | 16 | 3 | 0 | 0 | 2306 | 0 |
| RB | 16 | 1 | 0 | 0 | 2274 | 0 |
| RCUT | 16 | 1 | 18 (approximate inner: diagnostic) | 0 | 2274 | 0 |
| RCUT-t0.5 | 16 | 1 | 9 (approximate inner: diagnostic) | 0 | 2274 | 0 |
| RCUT-t2 | 16 | 1 | 37 (approximate inner: diagnostic) | 0 | 2274 | 0 |
| RCUT2 | 16 | 1 | 0 | 0 | 2274 | 0 |
| RCUT2-t0.5 | 16 | 1 | 0 | 0 | 2274 | 0 |
| RCUT2-t2 | 16 | 1 | 0 | 0 | 2274 | 0 |

## Time breakdown (medians; solver seconds summed over partitions and processes)

| method | workers | MISOCP s | SOCP s | LP s | other s | discovery s | cone_outer s | support s | coverage s | obbt s | ray s | zero_sp s | cut s | pool solve s | phase A s | phase B s | phase C s | phase CUT s | phase A+ s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H | 1 | 37.8 [36.8, 38.1] | 37.9 [36.8, 38] | 0.134 [0.133, 0.138] | 1.84 [1.82, 1.96] | 12.5 [12.1, 12.6] | 8.06 [7.8, 8.14] | 2.07 [2.06, 2.08] | 17.3 [16.9, 17.3] | 33.9 [32.9, 34] | 1.77 [1.71, 1.79] | 0.282 [0.266, 0.289] | 0 | 0 | 55.1 [53.3, 55.2] | 3.28 [3.27, 3.38] | 19 [18.5, 19] | - | - |
| H | 16 | 63.3 [58.7, 84] | 19.3 [18.6, 19.4] | 0.171 [0.151, 0.176] | 3.15 [3.06, 3.21] | 12.9 [12.7, 12.9] | 8.19 [8.03, 8.2] | 2.2 [1.99, 2.24] | 42.3 [38, 62.8] | 15.2 [14.5, 15.5] | 1.8 [1.76, 1.82] | 0.276 [0.271, 0.28] | 0 | 2.03 [1.84, 2.06] | 37.1 [36.8, 37.4] | 2.13 [2, 2.15] | 43.2 [40, 64.3] | - | - |
| R | 1 | 37.6 [36.8, 37.9] | 37.8 [37, 38.5] | 0 | 0.19 [0.167, 0.191] | 12.5 [12.3, 12.6] | 25.1 [24.6, 25.3] | 0 | 0 | 33.7 [33, 34.3] | 3.78 [3.68, 3.88] | 0.322 [0.315, 0.362] | 0 | 0 | - | - | - | - | - |
| R | 16 | 38.3 [37.4, 38.3] | 14.3 [14.2, 14.4] | 0 | 0.235 [0.23, 0.241] | 12.7 [12.5, 12.7] | 25.5 [24.9, 25.6] | 0 | 0 | 10.1 [10.1, 10.1] | 3.87 [3.78, 3.91] | 0.326 [0.319, 0.328] | 0 | 0 | - | - | - | - | - |
| RB | 16 | 29.5 | 17.7 | 0.143 | 2.71 | 13.9 | 15.6 | 1.81 | 0 | 12.8 | 2.88 | 0.308 | 0 | 1.67 | 37.7 | 1.8 | - | - | 8.52 |
| RCUT | 16 | 27.3 | 26 | 3.45 | 5.77 | 13.1 | 14.1 | 0 | 0 | 13.3 | 2.74 | 0.311 | 13.2 | 13.2 | 37.4 | - | - | 4.77 | 6.82 |
| RCUT-t0.5 | 16 | 28.4 | 27 | 3.81 | 6.55 | 13.1 | 15.3 | 0 | 0 | 13.2 | 2.88 | 0.314 | 14.5 | 14.5 | 37.4 | - | - | 5.51 | 8.11 |
| RCUT-t2 | 16 | 28.2 | 23.4 | 2.9 | 5.03 | 14.3 | 13.9 | 0 | 0 | 12.2 | 2.66 | 0.313 | 11.1 | 11.1 | 37.3 | - | - | 4.1 | 6.71 |
| RCUT2 | 16 | 28 | 41.7 | 3.42 | 12.1 | 12.9 | 15 | 0 | 0 | 13.3 | 18.4 | 0.445 | 12.9 | 28.6 | 37.2 | - | - | 11.1 | 7.98 |
| RCUT2-t0.5 | 16 | 28.2 | 44 | 3.77 | 13.3 | 13 | 15.1 | 0 | 0 | 13.4 | 19.8 | 0.449 | 14.2 | 31.2 | 37.4 | - | - | 12.2 | 8.02 |
| RCUT2-t2 | 16 | 28 | 38.7 | 2.93 | 11.2 | 13 | 15 | 0 | 0 | 13.3 | 16.9 | 0.457 | 11 | 25.2 | 37.1 | - | - | 10.2 | 7.97 |

## Counts (medians)

| method | workers | MISOCP | SOCP | LP | |X*| | cones | faces | unresolved faces |
|---|---|---|---|---|---|---|---|---|
| H | 1 | 54 | 233 | 126 | 19 | 10 | 167 | 0 |
| H | 16 | 56 [55, 57] | 310 [292, 315] | 158 [139, 161] | 20 [19, 20] | 10 | 197 [179, 200] | 0 |
| R | 1 | 75 | 122 | 0 | 0 | 23 | 0 | 0 |
| R | 16 | 75 | 122 | 0 | 0 | 23 | 0 | 0 |
| RB | 16 | 60 | 287 | 129 | 18 | 16 | 169 | 0 |
| RCUT | 16 | 56 | 508 | 189 | 18 | 14 | 0 | 0 |
| RCUT-t0.5 | 16 | 58 | 546 | 207 | 18 | 15 | 0 | 0 |
| RCUT-t2 | 16 | 56 | 450 | 160 | 18 | 14 | 0 | 0 |
| RCUT2 | 16 | 58 | 855 | 189 | 18 | 15 | 0 | 0 |
| RCUT2-t0.5 | 16 | 58 | 919 | 207 | 18 | 15 | 0 | 0 |
| RCUT2-t2 | 16 | 58 | 764 | 160 | 18 | 15 | 0 | 0 |

## Partitions (status, t_cert s)

- H w1 run1: pp certified 2.5 (support); pn certified 7.3 (radial); np certified 38.8 (support); nn certified 75.8 (support)
- H w1 run2: pp certified 2.5 (support); pn certified 7.5 (radial); np certified 39.8 (support); nn certified 78.0 (support)
- H w1 run3: pp certified 2.5 (support); pn certified 7.6 (radial); np certified 40.1 (support); nn certified 77.9 (support)
- H w16 run1: pp certified 3.1 (support); pn certified 4.5 (radial); np certified 47.9 (support); nn certified 25.9 (support)
- H w16 run2: pp certified 3.1 (support); pn certified 4.6 (radial); np certified 73.5 (support); nn certified 25.3 (support)
- H w16 run3: pp certified 3.1 (support); pn certified 4.5 (radial); np certified 50.8 (support); nn certified 26.8 (support)
- R w1 run1: pp certified 2.7; pn certified 7.6; np certified 33.6; nn certified 74.0
- R w1 run2: pp certified 2.7; pn certified 7.6; np certified 34.1; nn certified 75.5
- R w1 run3: pp certified 2.7; pn certified 7.7; np certified 34.6; nn certified 76.6
- R w16 run1: pp certified 2.1; pn certified 3.3; np certified 20.3; nn certified 28.4
- R w16 run2: pp certified 2.1; pn certified 3.4; np certified 20.6; nn certified 29.0
- R w16 run3: pp certified 2.2; pn certified 3.4; np certified 20.5; nn certified 29.2
- RB w16 run1: pp certified 3.4 (radial); pn certified 4.6 (radial); np certified 14.7 (inner); nn certified 27.8 (inner)
- RCUT w16 run1: pp certified 3.3 (cut); pn certified 4.4 (radial); np certified 14.4 (cut); nn certified 29.4 (cut)
- RCUT-t0.5 w16 run1: pp certified 3.3 (cut); pn certified 4.4 (radial); np certified 15.9 (cut); nn certified 29.9 (cut)
- RCUT-t2 w16 run1: pp certified 3.3 (cut); pn certified 4.4 (radial); np certified 14.1 (cut); nn certified 28.8 (cut)
- RCUT2 w16 run1: pp certified 3.8 (cut); pn certified 4.5 (radial); np certified 18.0 (cut); nn certified 32.4 (cut)
- RCUT2-t0.5 w16 run1: pp certified 4.0 (cut); pn certified 4.5 (radial); np certified 18.8 (cut); nn certified 33.0 (cut)
- RCUT2-t2 w16 run1: pp certified 3.8 (cut); pn certified 4.5 (radial); np certified 17.6 (cut); nn certified 31.9 (cut)

## RCUT / RCUT2 threshold sensitivity (final state)

| method | threshold | patience | workers | certified partitions | t_cert (s) | final gap | inner-AC FR% | inner-AC MR% | outer-AC FR% | inner ∧ SOCP-infeasible cells | phase CUT s | networks in inner / X* | shared cuts | SP | rays |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RCUT-t0.5 | 0.5% | 3 | 16 | 4/4 | 29.9 | 0.0057 | 0.935 | 0.043 | 1.360 | 9 | 5.5 | 18/18 | 207 | 451 | 0 |
| RCUT | 1% | 3 | 16 | 4/4 | 29.4 | 0.0070 | 0.989 | 0.025 | 1.431 | 18 | 4.8 | 18/18 | 189 | 416 | 0 |
| RCUT-t2 | 2% | 3 | 16 | 4/4 | 28.8 | 0.0056 | 1.115 | 0.018 | 1.431 | 37 | 4.1 | 18/18 | 160 | 358 | 0 |
| RCUT2-t0.5 | 0.5% | 3 | 16 | 4/4 | 33.0 | 0.0073 | 0.827 | 0.104 | 1.360 | 0 | 12.2 | 18/18 | 207 | 447 | 376 |
| RCUT2 | 1% | 3 | 16 | 4/4 | 32.4 | 0.0074 | 0.827 | 0.104 | 1.360 | 0 | 11.1 | 18/18 | 189 | 411 | 348 |
| RCUT2-t2 | 2% | 3 | 16 | 4/4 | 31.9 | 0.0078 | 0.815 | 0.110 | 1.360 | 0 | 10.2 | 18/18 | 160 | 354 | 314 |

## Per-network sets vs RB (same partition and scheme, final state)

RB's support phase gives P_x ⊆ R_x ⊆ O_x for the same tightened network region. For RCUT the excess vol(N_x)/vol(R_x)-1 therefore lies in [vol(N_x)/vol(O_x)-1, vol(N_x)/vol(P_x)-1]; for RCUT2, vol(N'_x)/vol(P_x)-1 compares two certified inner sets (negative: smaller than RB's P_x).

| method | quantity | networks | median | min | max |
|---|---|---|---|---|---|
| RCUT | N_x excess, lower bound | 18 | 0.22% | -0.39% | 2.13% |
| RCUT | N_x excess, upper bound | 18 | 0.71% | 0.03% | 2.35% |
| RCUT-t0.5 | N_x excess, lower bound | 18 | 0.16% | -0.39% | 1.79% |
| RCUT-t0.5 | N_x excess, upper bound | 18 | 0.45% | 0.03% | 2.01% |
| RCUT-t2 | N_x excess, lower bound | 18 | 0.53% | -0.21% | 6.15% |
| RCUT-t2 | N_x excess, upper bound | 18 | 0.94% | 0.04% | 6.45% |
| RCUT2 | N'_x vs RB P_x | 18 | -0.06% | -0.51% | 0.35% |
| RCUT2-t0.5 | N'_x vs RB P_x | 18 | 0.01% | -0.45% | 0.35% |
| RCUT2-t2 | N'_x vs RB P_x | 18 | -0.07% | -2.27% | 0.35% |

## Conclusions

- workers 1: R certified 3/3 runs, H 3/3; median final gap R 0.0075 vs H 0.0043; median t_cert R 75.5 s vs H 77.9 s; inner-AC MR R 0.276% vs H 0.067%, inner-AC FR R 0.750% vs H 0.809%. H phase medians (s): A 55.1, B 3.3, C 19.0.
- workers 16: R certified 3/3 runs, H 3/3; median final gap R 0.0075 vs H 0.0043; median t_cert R 29.0 s vs H 50.8 s; inner-AC MR R 0.276% vs H 0.080%, inner-AC FR R 0.750% vs H 0.827%. H phase medians (s): A 37.1, B 2.1, C 43.2.
- workers 16: R certified 3/3 runs, RB 1/1; median final gap R 0.0075 vs RB 0.0076; median t_cert R 29.0 s vs RB 27.8 s; inner-AC MR R 0.276% vs RB 0.141%, inner-AC FR R 0.750% vs RB 0.827%. RB phase medians (s): A 37.7, B 1.8, A+ 8.5.
- workers 16: R certified 3/3 runs, RCUT 1/1; median final gap R 0.0075 vs RCUT 0.0070; median t_cert R 29.0 s vs RCUT 29.4 s; inner-AC MR R 0.276% vs RCUT 0.025%, inner-AC FR R 0.750% vs RCUT 0.989%. RCUT phase medians (s): A 37.4, CUT 4.8, A+ 6.8.
- workers 16: R certified 3/3 runs, RCUT-t0.5 1/1; median final gap R 0.0075 vs RCUT-t0.5 0.0057; median t_cert R 29.0 s vs RCUT-t0.5 29.9 s; inner-AC MR R 0.276% vs RCUT-t0.5 0.043%, inner-AC FR R 0.750% vs RCUT-t0.5 0.935%. RCUT-t0.5 phase medians (s): A 37.4, CUT 5.5, A+ 8.1.
- workers 16: R certified 3/3 runs, RCUT-t2 1/1; median final gap R 0.0075 vs RCUT-t2 0.0056; median t_cert R 29.0 s vs RCUT-t2 28.8 s; inner-AC MR R 0.276% vs RCUT-t2 0.018%, inner-AC FR R 0.750% vs RCUT-t2 1.115%. RCUT-t2 phase medians (s): A 37.3, CUT 4.1, A+ 6.7.
- workers 16: R certified 3/3 runs, RCUT2 1/1; median final gap R 0.0075 vs RCUT2 0.0074; median t_cert R 29.0 s vs RCUT2 32.4 s; inner-AC MR R 0.276% vs RCUT2 0.104%, inner-AC FR R 0.750% vs RCUT2 0.827%. RCUT2 phase medians (s): A 37.2, CUT 11.1, A+ 8.0.
- workers 16: R certified 3/3 runs, RCUT2-t0.5 1/1; median final gap R 0.0075 vs RCUT2-t0.5 0.0073; median t_cert R 29.0 s vs RCUT2-t0.5 33.0 s; inner-AC MR R 0.276% vs RCUT2-t0.5 0.104%, inner-AC FR R 0.750% vs RCUT2-t0.5 0.827%. RCUT2-t0.5 phase medians (s): A 37.4, CUT 12.2, A+ 8.0.
- workers 16: R certified 3/3 runs, RCUT2-t2 1/1; median final gap R 0.0075 vs RCUT2-t2 0.0078; median t_cert R 29.0 s vs RCUT2-t2 31.9 s; inner-AC MR R 0.276% vs RCUT2-t2 0.110%, inner-AC FR R 0.750% vs RCUT2-t2 0.815%. RCUT2-t2 phase medians (s): A 37.1, CUT 10.2, A+ 8.0.
