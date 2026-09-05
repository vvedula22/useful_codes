# Courtemanche–Ramirez–Nattel (CRN) human atrial myocyte

Courtemanche M, Ramirez RJ, Nattel S. *Am J Physiol* 275:H301–H321, 1998.
Model keyword: `crn` (or `courtemanche`). nX = 6, nG = 15.

Run:

```
../../bin/cep.exe in_CRN.dat        # writes log_CRN.txt (26 columns)
```

`in_CRN.dat` reproduces the paper's stimulus (2-ms, −2 nA pulses at
1000 ms) for a single 600-ms beat from rest; set
`Number of time steps: 2400000` for the paper's full 12-s protocol.
`plot_CRN.m` documents the log column map.

## Units

Clock in **ms**; currents are **absolute pA** (the paper's pA/pF
densities multiplied by `Cm = 100 pF` at the point of use), so
`dV/dt = -ΣI/Cm` is in pA/pF = mV/ms and the stimulus amplitude is
absolute pA — **not** the pA/pF density a TTP deck uses.

These are the paper's own units, kept internally consistent with the
original formulation rather than rescaled to a common system. Because
the CEP models each do this, stimulus blocks and time steps are not
portable between decks — see [`../README_CEP_UNITS.md`](../README_CEP_UNITS.md).

CRN must stay in absolute pA: Eq. 68's `Fn` thresholds (the SR release
trigger) only fire with absolute currents — a pA/pF port silently
disables calcium release. The one deliberate deviation from the
published equations is `Fn_rel_scale`, needed to reproduce the paper's
Fig. 15 (default `0.01`; set `1.0` in `params_crn.in` for Eq. 68 as
published). Both are documented in `include/PARAMS_CRN.f`.

## Figures

Generated from this Fortran implementation (steady state, 12 s of
pacing), overlaid with the validated Python reference (black dashed)
and, where available, points digitised from the published figures:

| file | contents |
|---|---|
| `fig14_ap_currents.png` / `fig14_overlay.png` | AP + membrane currents (paper Fig. 14) |
| `fig15_ca_handling.png` / `fig15_overlay.png` | Ca dynamics (paper Fig. 15; overlay flux panel: model ×1/100 onto the published axis — the published panel is mis-scaled by 100×) |
| `fig16_rate_dependence.png` | BCL 1 / 0.6 / 0.3 s (paper Fig. 16) |
| `plot_CRN_12s_table3.png` | 12-s protocol + biomarkers vs paper Table 3 |

Fortran vs Python agreement: all 26 logged quantities identical to the
log's full 10-significant-digit precision at every step (max |ΔV| =
5e-9 mV), for FE, RK4 and both `Fn_rel_scale` settings.
