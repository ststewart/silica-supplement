# Attenuation-length correction of visible thermal emission from decaying shocks in silica

Scripts and data to re-run the attenuation-length fits and make the figures in the appendix
"Attenuation-Length Correction of Visible Thermal Emission from Decaying Shocks".

The model follows Ye et al. (2024, J. Appl. Phys. 136, 105904): the shocked sample is a stack
of Lagrangian slices that release isentropically behind the front, each with an attenuation
length set by its peak shock pressure. The visible radiance is the sum of the slice emission
transmitted to the free surface. The attenuation length is fitted empirically:

- **Fused silica, gas-gun shot 50:** shock-front temperature and attenuation length fitted
  jointly to the 620–660 nm radiance and the 4.8 µm temperature.
- **α-quartz, OMEGA decaying shocks 70813, 70817, 70820:** shock-front temperature fixed to the
  steady-shock Hugoniot; attenuation length fitted below the temperature inversion.

The pressure behind the front is taken equal to the front pressure (isobaric), which gives
lower limits on the attenuation length. For shot 50 the fit is repeated with the pressure
behind the front from a pyKO hydrocode simulation.

## Requirements

Python 3 with numpy, scipy and matplotlib.

## Running

```
python3 run_all.py
```

runs, in order (about 1 minute):

1. `omega_extract_time.py` – shock velocity and temperature vs time for the OMEGA shots, and the
   breakout times.
2. `shot50_visible_model.py` – fused silica shot 50, isobaric release.
3. `shot50_visible_model.py pyko` – fused silica shot 50, pressure behind the front from pyKO.
4. `omega_quartz_attenuation.py` – α-quartz OMEGA shots, isobaric release, with the sensitivity
   to the pressure behind the front. An optional argument moves the opacity onset, e.g.
   `python3 omega_quartz_attenuation.py 8.99` for the Hugoniot–solidus intersection.

All paths are relative to this folder; the scripts can be run from any directory.
Outputs are written to `output/`.

## Scripts

| File | Contents |
|---|---|
| `slab_emission.py` | Planck functions, stishovite Mie–Grüneisen release isentropes, slab radiative transfer |
| `omega_extract_time.py` | OMEGA Us(t), T(t) and breakout times |
| `shot50_visible_model.py` | Fused silica shot 50 fit |
| `omega_quartz_attenuation.py` | α-quartz OMEGA fit and comparison with fused silica |
| `run_all.py` | Runs the steps above in order |

## Appendix figures

| Appendix figure | Output file |
|---|---|
| Shot 50 radiance and temperature vs time | `output/shot50-vis-radiance-model.pdf`, `output/shot50-vis-T-time.pdf` |
| Shot 50 corrected temperature vs pressure | `output/shot50-vis-corrected-TP.pdf` |
| OMEGA apparent temperature vs shock velocity | `output/omega-qtz-T-Us-model.pdf` |
| OMEGA shock velocity and temperature vs time | `output/omega-extracted/omega-Us-T-time.pdf` |
| Attenuation length vs pressure | `output/omega-qtz-vs-fs-attenuation.pdf` |

The LaTeX appendix reads these PDFs from a `figs/` folder (same file names).

Numerical results are in the `.csv` files in `output/` and printed to the terminal. The
`shot50-pyko-*` outputs are the pyKO test; the fit results quoted in the appendix are those of
the isobaric runs.

## Data

`data/shot50/` – fused silica gas-gun shot 50 (copper flyer and driver, 5.298 mm sample per the shot record)
- `S050_Lineout_timeprofile_620-660nm.txt` – calibrated spectral radiance, 620–660 nm, vs time (µs).
- `S050_TemperatureProfile_-1000-1000ns_20nsWindow.txt` – streak-camera temperature profile
  (20 ns window), used for comparison only.
- `S050_IR_4.8um_radiance.csv` – 4.8 µm spectral radiance from the IR pyrometer analysis
  (time in s, radiance in W sr⁻¹ m⁻³), −1 to 2 µs.

`data/pyko/` – pyKO simulation of shot 50
- `shot50-cu-fs-decay-aneos.yml` – pyKO configuration (polycarbonate sabot, Cu flyer at
  5.353 km/s, Cu driver, fused silica; ANEOS equations of state).
- `shot50-pyko-aneos-fs.npz` – fused silica cells extracted from the pyKO output: arrays
  [time, cell] of position, pressure, temperature, density, particle velocity, entropy and
  ANEOS phase flag (SI units), with the time and initial cell positions.

`data/omega/<shot>/` – OMEGA α-quartz decaying shocks (shots 70813, 70817, 70820)
- `<shot>_Us1.txt`, `<shot>_Us2.txt` – shock velocity (km/s) vs time (ns) from the two VISAR
  channels, from the calibrated analysis of M. Millot. The time axes include the streak-camera
  sweep calibration. **These define the time base.**
- `<shot>_T.txt` – SOP temperature (eV) vs time (ns) on the same clock (M. Millot). Plotted
  only; not fitted.
- `<shot>_TandUs.txt` – the original decaying-shock temperature (kK) and its uncertainty vs
  shock velocity (km/s). **These are the data fitted in the analysis.**
- `omega-thicknesses.png` – target dimensions; the α-quartz thicknesses are 253.7, 278.0 and
  418.1 µm for shots 70813, 70817 and 70820.

Data sources for the OMEGA analysis:

| Quantity | Source | Use |
|---|---|---|
| Time base, slices: U_s(t) | VISAR channel 1, M. Millot (`<shot>_Us1.txt`) | Slice thicknesses and peak states |
| Breakout time | VISAR channel 2, M. Millot (`<shot>_Us2.txt`) | Table of transit times, figure |
| Fitted temperatures and uncertainties | Original T(U_s) records (`<shot>_TandUs.txt`) | Attenuation-length fit |
| SOP temperature vs time | M. Millot (`<shot>_T.txt`) | Figure only |

`omega_extract_time.py` takes the shock velocity from VISAR channel 1, from shock entry into
the quartz to the end of the reliable record. The VISAR record ends before breakout, when
the shock front loses its reflectivity; breakout is taken from VISAR channel 2, which records
the jump to the free-surface velocity. The distance travelled by the shock at the measured
shock velocity is consistent with the quartz thickness at breakout. The attenuation-length
fit compares the model with the `<shot>_TandUs.txt` temperatures and uncertainties at the
same shock velocity, using only the data within the range of the measured shock velocity
(this drops the 70813 points at 7.00–7.46 km/s). Millot's SOP T(t), evaluated along his
U_s(t), agrees with these records to a median of 13–25 K.

## Model inputs taken from the main text

Hugoniot fits for fused silica and α-quartz (Table 2), the stishovite Mie–Grüneisen model
(Appendix B; γ₀ = 1.22, q = 2.22, c_v = 3.65 × 3R/M, 4000 K isotherm of Driver et al. 2010)
and the shot 50 shock-velocity history are written into the scripts.
