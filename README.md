# Supplemental Materials for 'Transient Melt in the Shock Front: Reinterpreting SiO₂ Shock Temperatures and the Melt Curve'


This GitHub repository contains supplemental material and replication data for:<br>

Amodeo, K. M., Davies, E. J., Spaulding, D. K., Millot, M., Stewart, S. T., Transient Melt in the Shock Front: Reinterpreting SiO₂ Shock Temperatures and the Melt Curve, <i>J. Geophysical Research: Solid Earth</i>, 2026 (revised). (Zenodo DOI tbd)<p>

The supplementary materials are collected in a Jupyter book that is deployed at:<br>
https://ststewart.github.io/silica-supplement/

## Repository layout
```
data_silica_db/      literature and new shock/static data (csv/txt), read by silica-book/scripts/load_silica_data.py
py-rusbank/          RUSbank Shock Database
silica-book/         Jupyter book: EOS model, Hugoniot and melt-curve fits, manuscript figures
  scripts/           shared modules (set_paths, define_eos, load_*, PT/PVphasediagram.py, melt-curve.txt)
  model-fits/        plain-text fit records (coefficients, covariance, data) for Tables 2, 4, F1
  paper-figs/        PDFs of every python manuscript figure
  paper-schematics/  PDF schematics
  imp-match-UCD/     Monte Carlo impedance match, UC Davis gas-gun shots (Table 1), Cu and SS304 fits (Table F1)
  imp-match-Lyzenga/ recalibrated Lyzenga et al. (1983)/Boslough (1988) data, Ta fit (Table F1, Fig. F1)
  imp-match-JA79/    Jackson & Ahrens (1979) fused silica data (excluded from fits; Appendix A)
  shot50/, shot60/   gas-gun decaying-shock data
attenuation-length/  Appendix G: attenuation-length model, data and figures G1–G5 (python3 run_all.py);
                     stishovite-model Hugoniot used in Fig. 8A
```

## Figure index

N = `silica-book/additional-plots.ipynb` (the revised notebook). The `pind` values select the figure inside
`scripts/PTphasediagram.py` or `scripts/PVphasediagram.py`.

| Fig. | File in `paper-figs/` | Made by |
|---|---|---|
| 1A | `A23-paper-fig1a.pdf` | N, "Figure 1A" (`pind=23`, PTphasediagram.py) |
| 1B | `B1-SiO2-P-rho-phase-diagram-Fig1b.pdf` | N, "Figure 1B" (`pind=1`, PVphasediagram.py) |
| 2A, 2B | `A25-paper-fig2a.pdf`, `A24-paper-fig2b.pdf` | N, "Figure 2AB" (`pind` 24, 25) |
| 3 | `sample_schematic.pdf` | drawing (no code) in paper-schematics/ |
| 4 | `Fig4-all_fits-usuponly-small.pdf` | N, "Figure 4" |
| 5 | `Fig5-combined-small.pdf` | N, "Figure 5 combined" |
| 6 | `penetration depth cartoon.pdf` | drawing (no code) in paper-schematics/ |
| 7A, 7B | `A27-paper-fig7a.pdf`, `A26-paper-fig7b.pdf` | N, "Figure 7ABC" (`pind` 26, 27) |
| 8A | `Fig8A-newcombo.pdf` | N, "Figure 8A" (`pind=910`, PVphasediagram.py) |
| 8B | `Fig8B-newcombo.pdf` | N, "Figure 8B" (cell saving `Fig8B-newcombo`) |
| B1 | `A30-paper-fig9.pdf` | N, "Supplemental Figures" (`pind=30`) |
| B2 | `FIG-DFT-PBE-offset.pdf` | `dft-offset.ipynb`, "Figure for paper" |
| B3 | `MC-dV-comp.pdf` | `melt-curve-fits.ipynb`, "Volume change on melting" |
| C1 | `AP-A-wiggles.pdf` | N, "Supplemental Figures" |
| C2 | `AP-AB-omega-kink.pdf` | N, "Inflection in the quartz shock temperatures" |
| C3 | `shot50-sop-48-comp.pdf`, `shot50-4.8-model.pdf` | `shot50-decay.ipynb` |
| D1 | `cst-fits-all.pdf` | N, "Coesite supplemental figures" |
| D2 | `A28-paper-fig7c.pdf` | N, "Figure 7ABC" (`pind=28`) |
| D3 | `B13-SiO2-P-rho-phase-diagram-cst-Fig8C.pdf` | N, "Alternate version of Figure 8ABC" (`pind=13`) |
| D4 | `Stv-fits.pdf` | N, "Stishovite supplemental figures" |
| D5 | `stv-melt-crossing2.pdf`, `stv-melt-crossing3.pdf` | `hugoniot-fits-stishovite.ipynb`, "Using the linear fit to up<18 km/s" |
| E1 | `fs-cs.pdf` | N, "Stishovite supplemental figures" (`fs-cs` cell) |
| F1 | `Boslough-redo.pdf` | N, "Supplemental Figures" |
| G1 | `shot50-vis-radiance-model.pdf`, `shot50-vis-T-time.pdf` | `attenuation-length/shot50_visible_model.py` → `attenuation-length/output/` |
| G2 | `shot50-vis-corrected-TP.pdf` | same |
| G3 | `omega-qtz-T-Us-model.pdf` | `attenuation-length/omega_quartz_attenuation.py` |
| G4 | `omega-Us-T-time.pdf` | `attenuation-length/omega_extract_time.py` → `output/omega-extracted/` |
| G5 | `omega-qtz-vs-fs-attenuation.pdf` | `attenuation-length/omega_quartz_attenuation.py` |
| H1 | `fs-cs-high.pdf` | rerun notebooks with higher temperature melt curve |
| H2 | `B9-…-FS-Fig8A-high.pdf`, `B10-…-qtz-Fig8B-high.pdf`, `B13-…-cst-Fig8C-high.pdf` | re-run notebooks with higher temperature melt curve |



## Fit and table index

| Manuscript | Fit record (`silica-book/model-fits/`) | Made by |
|---|---|---|
| Table 1 (gas-gun data) | input `imp-match-UCD/quartz-all-may2025.csv`; MC output `MCIM-Silica-output-*.txt` | Monte Carlo impedance match in `imp-match-UCD/`; table printed in `tables.ipynb` |
| Table 2, fused silica | `fs_fit_stv_linear`, `fs_fitTUs_stv_quadratic`, `fs_fit_liq_quadratic`, `fs_fitTus_liq_cubic` | `hugoniot-fits-fused-silica.ipynb` |
| Table 2, α-quartz | `qtz_fit_solid_linear`, `qtz_fitTus_solid_quadratic`, `qtz_fit_liq_cubic`, `qtz_fitTus_liq_cubic` | `hugoniot-fits-quartz.ipynb` |
| Table 2, coesite | `cst_fit_cacl2_linear`, `cst_fitTus_cacl2_quadratic` (model-EOS temperatures; zero covariance), `cst_fit_liq_linear`, `cst_fitTUs_liq_quadratic` | `hugoniot-fits-coesite.ipynb` |
| Table 2, stishovite | `stv_fit_liq_linear_up18`, `stv_fitTus_liq_fourthpoly` | `hugoniot-fits-stishovite.ipynb` |
| Tables 3, D1 | (computed) | `melt-curve-fits.ipynb`, "Hugoniot crossing the melt curve"; `tables.ipynb` |
| Table 4, β-cristobalite | constant 1999 K (Bourova & Richet 1998) | `melt-curve-fits.ipynb` |
| Table 4, β-quartz | `SiO2_fithighqtzstatic_simons` (cubic) | `melt-curve-fits.ipynb`, "Stishovite low-pressure region melt curve" |
| Table 4, coesite | `SiO2_fitcoesitestatic_simons` (cubic) | same |
| Table 4, stishovite–CaCl₂ | `SiO2_fitstvstatic_simons` | same |
| Table 4, GM24 | none; published fit (Geng & Mohn 2024), coded in `melt-curve-fits.ipynb`, "Print the fiducial melt curve" | — |
| Table 4, GDG16* | `SiO2_fitG16mod_simons` (lowest point raised to meet GM24); unmodified fit `SiO2_fitG16_simons` | `melt-curve-fits.ipynb`, "Gonzales et al. 2016 - tweak to intersect GM24" |
| Model melt curve (tabulated) | `silica-book/scripts/melt-curve.txt` | `melt-curve-fits.ipynb`, "Print the fiducial melt curve" |
| Liquidus density fit (Fig. B2) | `rho_l_fit_quadratic` | `melt-curve-fits.ipynb` |
| Table F1 | `ss304_fit_linear`; Cu and Ta coefficients in notebook output | `recalibrate-data.ipynb` writes `ss304_fit_linear` (so does `imp-match-UCD/Al2024-fit.ipynb`; keep one, and check that `steel304-fit.ipynb` agrees); `imp-match-UCD/copper101-fit.ipynb`; `imp-match-Lyzenga/fit-ta.ipynb` |
| Table H1 | (computed) | rerun notebooks with higher temperature melt curve |
| Appendix G results | `attenuation-length/output/*.csv` | `attenuation-length/run_all.py` |
| Appendix I covariance matrices | the `cov:` field of each `model-fits/*-readable.txt` file above | — |

Each `*-readable.txt` file lists the equation, coefficients, covariance, fit range, N, the data used, and their
references.


## Regenerating the Jupyterbook

To run the Jupyter notebooks, build the conda environment:<p>

conda env create -f environment-silicapaper.yml<br>
conda activate silicapaper<p>

To rebuild the html Jupyterbook:<p>

jupyter-book build silica-book<p>

The output is in:<p>

silica-book/_build/html<p>

S. T. Stewart<br>
sstewa56@asu.edu<br>
10/23/2026<br>
