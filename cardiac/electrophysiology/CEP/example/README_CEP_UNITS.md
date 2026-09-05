# CEP model unit systems

**The cellular activation models in this solver do not share a unit
system.** They differ in the clock, in whether currents are absolute or
densities, and in what role `Cm` plays in `dV/dt`. Stimulus blocks, time
steps and parameter values are therefore **not portable between input
decks**.

This file is the authoritative reference; the module headers and example
decks point here rather than repeat it.

## The table

| model | keyword | clock | currents | `Cm` | `dV/dt` | stimulus amplitude | example deck |
|---|---|---|---|---|---|---|---|
| Aliev–Panfilov | `ap` | ms | — dimensionless | 1, unused | dimensionless, rescaled by `Vscale` / `Tscale` / `Voffset` | mV/ms | `52.0`, dt `0.2` |
| Bueno-Orovio | `bo` | ms | — dimensionless | 1, unused | dimensionless, rescaled by `Vscale` / `Tscale` / `Voffset` | mV/ms | `-38.0`, dt `0.1` |
| FitzHugh–Nagumo | `fn` | deck's unit | — dimensionless | none | fully dimensionless | model units | `0.0`, dt `0.01` |
| ten Tusscher–Panfilov | `ttp` | ms | densities **pA/pF** | 0.185 | −Σ I (no `Cm`) | **pA/pF** | `-38.0`, dt `0.2` |
| Stewart Purkinje | `pfib` | ms | densities **pA/pF** | 0.185 | −Σ I (no `Cm`) | **pA/pF** | `-38.0`, dt `0.2` |
| Tong uterine | `tong` | ms | densities **pA/pF** | 1 µF/cm² | −Σ I (no `Cm`) | **pA/pF** | `-0.2` … `-0.5`, dt `1.0` / `0.2` |
| Nygren atrial | `nyg` | **s** | absolute **pA** | 0.05 **nF** | −Σ I / `Cm` | **pA** | `-280.0`, dt `1e-4` |
| Courtemanche–Ramirez–Nattel | `crn` | ms | absolute **pA** | 100 **pF** | −Σ I / `Cm` | **pA** | `-2000.0`, dt `0.005` |

Long-form keywords also work: `aliev_panfilov`, `bueno_orovio`,
`fitzhugh_nagumo`, `tentusscher_panfilov`, `purkinje`, `nygren`,
`courtemanche`. (`decoupled` selects excitation–contraction only, with
no cell model.)

## Why the clocks differ

The time base follows from the units of `Cm`:

- `pA / pF = mV/ms` → **CRN** runs in milliseconds
- `pA / nF = mV/s` → **NYG** runs in seconds
- TTP, PFIB and TONG carry current *densities*, which are already mV/ms,
  so their `dV/dt` has no `Cm` in it at all

Each model is kept internally consistent with the formulation in its own
paper, which is why the values are not rescaled to a common system.

## `Cm` means three different things

- **NYG, CRN** — divides `dV/dt`; this is what fixes the clock. CRN also
  multiplies it into the Table 1 conductances to turn densities into
  absolute currents.
- **TTP, PFIB, TONG** — appears only in the concentration updates, never
  in `dV/dt`. (`PARAMS_TTP.f` annotates its `Cm = 0.185` as µF/cm²
  while noting the CellML source gives 0.185 µF; the value is used as
  CellML defines it.)
- **AP, BO** — read from the parameters file and then never used by the
  single-cell equations. It is there for the tissue solver, which needs
  `Cm`, `sV` and `rho` to form the monodomain conductivity.

## What goes wrong in practice

A TTP amplitude of `-38` is a *density* — roughly −3800 pA absolute for a
100 pF cell. CRN's `-2000` is absolute pA. Paste one into the other's
deck and you are off by about two orders of magnitude.

Copying NYG's `Time step size: 0.0001` (seconds) into a CRN deck is worse:
NYG means 0.1 ms, but CRN reads the bare number as 0.0001 ms, so each
step is 1000× too small and the 0.006 stimulus lasts 6 µs instead of
6 ms. It will not error — it produces a flat trace.

For tissue simulation the conductivity follows the clock too: mm²/ms for
the millisecond models, mm²/s for NYG.
