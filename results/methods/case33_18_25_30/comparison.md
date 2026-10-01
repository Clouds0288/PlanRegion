# Plan 2 comparison: case33_18_25_30

Commit `4c5b2a407f4b708d3c5adc1460b206e702b75c4c`; case33bw load nodes [18, 25, 30], mode 1, budget 7; partitions ppp, ppn, pnp, pnn, npp, npn, nnp, nnn.

Settings: total 300 s wall clock per run; tau 0.005; eps 0.015; H: eps_A 0.15, share_A 0.25, eps_B 0.0075, criterion auto, coverage global; Gurobi threads 1, seed 0, MISOCP cap 60 s, MIPGap 0.001; Gurobi 13.0.2.

Primary reference: AC scan (FR/MR in %, undecided AC cells excluded). SOCP columns are diagnostics: the OBBT-tightened model removes SOCP-only area, so outer-vs-SOCP MR > 0 is expected.

RCUT / RCUT2 (R + mainline cutting-plane construction instead of phase B): per network, SP scores the vertices of N_x (fixed x with its OBBT box and envelope rows) and the max-eta vertex yields a joint cut; a network stops after `patience` consecutive cuts each removing < `threshold` of vol(N_x). RCUT takes the stagnated/exact N_x as the network region (not a certified inner set); RCUT2 also runs the mainline inner certification (rays) and uses the certified hull. Inner I = (I_R ∪ sets) ∩ O_R; certified when vol(O_R∩box) - vol(I) <= eps·vol(I) (RB criterion). Runs: RCUT (threshold 1%, patience 3), RCUT-t0.5 (threshold 0.5%, patience 3), RCUT-t2 (threshold 2%, patience 3), RCUT2 (threshold 1%, patience 3), RCUT2-t0.5 (threshold 0.5%, patience 3), RCUT2-t2 (threshold 2%, patience 3).

## Runs

| method | workers | run | status | certified partitions | t_cert (s) | final gap | t gap<=10% | t gap<=5% | t gap<=2% | valid |
|---|---|---|---|---|---|---|---|---|---|---|
| R | 16 | 1 | time_limit | 3/8 | - | 0.0369 | 116.1 | 190.5 | - | yes |
| RB | 16 | 1 | time_limit | 3/8 | - | 0.0235 | 107.1 | 169.0 | - | yes |
| RCUT | 16 | 1 | time_limit | 7/8 | - | 0.0170 | 88.4 | 126.7 | 219.9 | yes |
| RCUT-t0.5 | 16 | 1 | time_limit | 7/8 | - | 0.0165 | 105.6 | 138.2 | 252.0 | yes |
| RCUT-t2 | 16 | 1 | certified | 8/8 | 272.2 | 0.0145 | 95.3 | 116.1 | 177.8 | yes |
| RCUT2 | 16 | 1 | time_limit | 3/8 | - | 0.0364 | 130.4 | 192.8 | - | yes |
| RCUT2-t0.5 | 16 | 1 | time_limit | 3/8 | - | 0.0322 | 140.0 | 193.3 | - | yes |
| RCUT2-t2 | 16 | 1 | time_limit | 3/8 | - | 0.0326 | 126.8 | 183.6 | - | yes |

## Medians [min, max] over repeats

| method | workers | runs certified | t_cert (s) | final gap | t gap<=10% | t gap<=5% | t gap<=2% | inner-AC FR% | inner-AC MR% | outer-AC FR% | outer-AC MR% | inner-SOCP FR% | inner-SOCP MR% | outer-SOCP FR% | outer-SOCP MR% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R | 16 | 0/1 | - | 0.0369 | 116 | 191 | - | 2.42 | 0.768 | 4.2 | 0 | 0 | 21 | 0.833 | 19.6 |
| RB | 16 | 0/1 | - | 0.0235 | 107 | 169 | - | 2.58 | 0.164 | 4.35 | 0 | 0 | 20.4 | 0.919 | 19.6 |
| RCUT | 16 | 0/1 | - | 0.017 | 88.4 | 127 | 220 | 4.28 | 0 | 5.51 | 0 | 0.709 | 19.5 | 1.49 | 19.1 |
| RCUT-t0.5 | 16 | 0/1 | - | 0.0165 | 106 | 138 | 252 | 3.6 | 0 | 4.83 | 0 | 0.53 | 19.9 | 1.15 | 19.4 |
| RCUT-t2 | 16 | 1/1 | 272 | 0.0145 | 95.3 | 116 | 178 | 4.54 | 0 | 5.55 | 0 | 0.919 | 19.4 | 1.41 | 19 |
| RCUT2 | 16 | 0/1 | - | 0.0364 | 130 | 193 | - | 2.55 | 0.286 | 4.75 | 0 | 0 | 20.6 | 1.24 | 19.5 |
| RCUT2-t0.5 | 16 | 0/1 | - | 0.0322 | 140 | 193 | - | 2.57 | 0.261 | 4.79 | 0 | 0 | 20.5 | 1.12 | 19.4 |
| RCUT2-t2 | 16 | 0/1 | - | 0.0326 | 127 | 184 | - | 2.55 | 0.296 | 4.57 | 0 | 0 | 20.6 | 1.01 | 19.5 |

## Validity (final state, cells)

| method | workers | run | inner ∧ SOCP-infeasible | AC-feasible ∖ outer | SOCP-feasible ∖ outer (diagnostic) | undecided AC |
|---|---|---|---|---|---|---|
| R | 16 | 1 | 0 | 0 | 36922 | 0 |
| RB | 16 | 1 | 0 | 0 | 36814 | 0 |
| RCUT | 16 | 1 | 1081 (approximate inner: diagnostic) | 0 | 35842 | 0 |
| RCUT-t0.5 | 16 | 1 | 802 (approximate inner: diagnostic) | 0 | 36413 | 0 |
| RCUT-t2 | 16 | 1 | 1405 (approximate inner: diagnostic) | 0 | 35668 | 0 |
| RCUT2 | 16 | 1 | 0 | 0 | 36681 | 0 |
| RCUT2-t0.5 | 16 | 1 | 0 | 0 | 36419 | 0 |
| RCUT2-t2 | 16 | 1 | 0 | 0 | 36607 | 0 |

## Time breakdown (medians; solver seconds summed over partitions and processes)

| method | workers | MISOCP s | SOCP s | LP s | other s | discovery s | cone_outer s | support s | coverage s | obbt s | ray s | zero_sp s | cut s | pool solve s | phase A s | phase B s | phase C s | phase CUT s | phase A+ s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R | 16 | 1592 | 185 | 0 | 12 | 15.3 | 1576 | 0 | 0 | 62.1 | 121 | 1.23 | 0 | 0 | - | - | - | - | - |
| RB | 16 | 1393 | 198 | 1.8 | 187 | 18.3 | 1375 | 17.5 | 0 | 70.7 | 110 | 1.4 | 0 | 15.7 | 438 | 20.4 | - | - | 1283 |
| RCUT | 16 | 871 | 232 | 15.8 | 151 | 16.7 | 855 | 0 | 0 | 68.6 | 67 | 1.05 | 111 | 111 | 416 | - | - | 67 | 660 |
| RCUT-t0.5 | 16 | 1010 | 270 | 19.6 | 187 | 15.8 | 994 | 0 | 0 | 70 | 76.9 | 1.1 | 142 | 142 | 416 | - | - | 84.9 | 828 |
| RCUT-t2 | 16 | 836 | 213 | 13 | 132 | 16 | 820 | 0 | 0 | 69.2 | 64.6 | 1.08 | 91.5 | 91.5 | 417 | - | - | 53.2 | 616 |
| RCUT2 | 16 | 1284 | 400 | 16.8 | 352 | 16.1 | 1267 | 0 | 0 | 70.9 | 226 | 1.93 | 118 | 244 | 416 | - | - | 176 | 1191 |
| RCUT2-t0.5 | 16 | 1253 | 444 | 21.3 | 386 | 16.6 | 1236 | 0 | 0 | 70.1 | 242 | 1.9 | 151 | 295 | 417 | - | - | 205 | 1161 |
| RCUT2-t2 | 16 | 1324 | 372 | 14.1 | 315 | 16.5 | 1307 | 0 | 0 | 70.5 | 217 | 1.93 | 97.1 | 211 | 417 | - | - | 160 | 1225 |

## Counts (medians)

| method | workers | MISOCP | SOCP | LP | |X*| | cones | faces | unresolved faces |
|---|---|---|---|---|---|---|---|---|
| R | 16 | 2236 | 2757 | 0 | 0 | 1124 | 0 | 0 |
| RB | 16 | 1946 | 4218 | 1355 | 63 | 974 | 3702 | 0 |
| RCUT | 16 | 1240 | 5486 | 804 | 63 | 610 | 0 | 0 |
| RCUT-t0.5 | 16 | 1437 | 6819 | 1023 | 64 | 711 | 0 | 0 |
| RCUT-t2 | 16 | 1196 | 4828 | 677 | 63 | 587 | 0 | 0 |
| RCUT2 | 16 | 1873 | 9031 | 834 | 65 | 939 | 0 | 0 |
| RCUT2-t0.5 | 16 | 1830 | 10380 | 1035 | 65 | 917 | 0 | 0 |
| RCUT2-t2 | 16 | 1930 | 8209 | 698 | 65 | 969 | 0 | 0 |

## Partitions (status, t_cert s)

- R w16 run1: ppp certified 14.4; ppn time_limit; pnp certified 147.6; pnn time_limit; npp certified 130.9; npn time_limit; nnp time_limit; nnn time_limit
- RB w16 run1: ppp certified 16.4 (inner); ppn time_limit; pnp certified 153.2 (inner); pnn time_limit; npp certified 98.5 (inner); npn time_limit; nnp time_limit; nnn time_limit
- RCUT w16 run1: ppp certified 8.9 (cut); ppn certified 149.8 (cut); pnp certified 107.3 (cut); pnn certified 211.3 (cut); npp certified 46.6 (cut); npn certified 190.6 (cut); nnp time_limit; nnn certified 149.5 (cut)
- RCUT-t0.5 w16 run1: ppp certified 8.9 (cut); ppn certified 182.5 (cut); pnp certified 120.8 (cut); pnn certified 214.6 (cut); npp certified 62.5 (cut); npn certified 194.3 (cut); nnp time_limit; nnn certified 266.0 (cut)
- RCUT-t2 w16 run1: ppp certified 9.0 (cut); ppn certified 145.9 (cut); pnp certified 107.1 (cut); pnn certified 209.0 (cut); npp certified 36.7 (cut); npn certified 179.5 (cut); nnp certified 272.2 (cut); nnn certified 148.4 (cut)
- RCUT2 w16 run1: ppp certified 18.0 (cut); ppn time_limit; pnp certified 159.3 (cut); pnn time_limit; npp certified 136.2 (cut); npn time_limit; nnp time_limit; nnn time_limit
- RCUT2-t0.5 w16 run1: ppp certified 17.9 (cut); ppn time_limit; pnp certified 164.1 (cut); pnn time_limit; npp certified 131.2 (cut); npn time_limit; nnp time_limit; nnn time_limit
- RCUT2-t2 w16 run1: ppp certified 18.1 (cut); ppn time_limit; pnp certified 158.4 (cut); pnn time_limit; npp certified 141.4 (cut); npn time_limit; nnp time_limit; nnn time_limit

## RCUT / RCUT2 threshold sensitivity (final state)

| method | threshold | patience | workers | certified partitions | t_cert (s) | final gap | inner-AC FR% | inner-AC MR% | outer-AC FR% | inner ∧ SOCP-infeasible cells | phase CUT s | networks in inner / X* | shared cuts | SP | rays |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RCUT-t0.5 | 0.5% | 3 | 16 | 7/8 | - | 0.0165 | 3.597 | 0.000 | 4.830 | 802 | 84.9 | 64/64 | 1023 | 4987 | 0 |
| RCUT | 1% | 3 | 16 | 7/8 | - | 0.0170 | 4.282 | 0.000 | 5.512 | 1081 | 67.0 | 63/63 | 804 | 3886 | 0 |
| RCUT-t2 | 2% | 3 | 16 | 8/8 | 272.2 | 0.0145 | 4.542 | 0.000 | 5.546 | 1405 | 53.2 | 63/63 | 677 | 3279 | 0 |
| RCUT2-t0.5 | 0.5% | 3 | 16 | 3/8 | - | 0.0322 | 2.567 | 0.261 | 4.791 | 0 | 205.0 | 65/65 | 1035 | 5019 | 3039 |
| RCUT2 | 1% | 3 | 16 | 3/8 | - | 0.0364 | 2.550 | 0.286 | 4.749 | 0 | 176.3 | 65/65 | 834 | 4006 | 2656 |
| RCUT2-t2 | 2% | 3 | 16 | 3/8 | - | 0.0326 | 2.545 | 0.296 | 4.571 | 0 | 160.3 | 65/65 | 698 | 3358 | 2412 |

## Per-network sets vs RB (same partition and scheme, final state)

RB's support phase gives P_x ⊆ R_x ⊆ O_x for the same tightened network region. For RCUT the excess vol(N_x)/vol(R_x)-1 therefore lies in [vol(N_x)/vol(O_x)-1, vol(N_x)/vol(P_x)-1]; for RCUT2, vol(N'_x)/vol(P_x)-1 compares two certified inner sets (negative: smaller than RB's P_x).

| method | quantity | networks | median | min | max |
|---|---|---|---|---|---|
| RCUT | N_x excess, lower bound | 63 | 4.69% | 0.05% | 14.44% |
| RCUT | N_x excess, upper bound | 63 | 5.40% | 0.34% | 15.28% |
| RCUT-t0.5 | N_x excess, lower bound | 63 | 2.57% | -0.06% | 12.54% |
| RCUT-t0.5 | N_x excess, upper bound | 63 | 3.28% | 0.22% | 13.30% |
| RCUT-t2 | N_x excess, lower bound | 63 | 7.07% | 0.05% | 19.86% |
| RCUT-t2 | N_x excess, upper bound | 63 | 7.71% | 0.71% | 20.69% |
| RCUT2 | N'_x vs RB P_x | 63 | -1.88% | -13.44% | -0.13% |
| RCUT2-t0.5 | N'_x vs RB P_x | 63 | -0.75% | -13.20% | 0.03% |
| RCUT2-t2 | N'_x vs RB P_x | 63 | -2.99% | -13.75% | -0.11% |

## Conclusions

- workers 16: R certified 0/1 runs, RB 0/1; median final gap R 0.0369 vs RB 0.0235; median t_cert R - s vs RB - s; inner-AC MR R 0.768% vs RB 0.164%, inner-AC FR R 2.421% vs RB 2.579%. RB phase medians (s): A 437.7, B 20.4, A+ 1282.6.
- workers 16: R certified 0/1 runs, RCUT 0/1; median final gap R 0.0369 vs RCUT 0.0170; median t_cert R - s vs RCUT - s; inner-AC MR R 0.768% vs RCUT 0.000%, inner-AC FR R 2.421% vs RCUT 4.282%. RCUT phase medians (s): A 416.5, CUT 67.0, A+ 659.6.
- workers 16: R certified 0/1 runs, RCUT-t0.5 0/1; median final gap R 0.0369 vs RCUT-t0.5 0.0165; median t_cert R - s vs RCUT-t0.5 - s; inner-AC MR R 0.768% vs RCUT-t0.5 0.000%, inner-AC FR R 2.421% vs RCUT-t0.5 3.597%. RCUT-t0.5 phase medians (s): A 416.2, CUT 84.9, A+ 827.9.
- workers 16: R certified 0/1 runs, RCUT-t2 1/1; median final gap R 0.0369 vs RCUT-t2 0.0145; median t_cert R - s vs RCUT-t2 272.2 s; inner-AC MR R 0.768% vs RCUT-t2 0.000%, inner-AC FR R 2.421% vs RCUT-t2 4.542%. RCUT-t2 phase medians (s): A 417.1, CUT 53.2, A+ 616.3.
- workers 16: R certified 0/1 runs, RCUT2 0/1; median final gap R 0.0369 vs RCUT2 0.0364; median t_cert R - s vs RCUT2 - s; inner-AC MR R 0.768% vs RCUT2 0.286%, inner-AC FR R 2.421% vs RCUT2 2.550%. RCUT2 phase medians (s): A 416.4, CUT 176.3, A+ 1190.8.
- workers 16: R certified 0/1 runs, RCUT2-t0.5 0/1; median final gap R 0.0369 vs RCUT2-t0.5 0.0322; median t_cert R - s vs RCUT2-t0.5 - s; inner-AC MR R 0.768% vs RCUT2-t0.5 0.261%, inner-AC FR R 2.421% vs RCUT2-t0.5 2.567%. RCUT2-t0.5 phase medians (s): A 417.2, CUT 205.0, A+ 1161.0.
- workers 16: R certified 0/1 runs, RCUT2-t2 0/1; median final gap R 0.0369 vs RCUT2-t2 0.0326; median t_cert R - s vs RCUT2-t2 - s; inner-AC MR R 0.768% vs RCUT2-t2 0.296%, inner-AC FR R 2.421% vs RCUT2-t2 2.545%. RCUT2-t2 phase medians (s): A 416.5, CUT 160.3, A+ 1224.7.
