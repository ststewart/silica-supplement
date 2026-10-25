# fs_stv_hugoniot.py
#
# Model Hugoniot of fused silica (2.204 g/cm3) transformed to stishovite, from 34 GPa (single
# shock wave; Wackerle 1962) to the melt curve, using the stishovite model in
# stv_model_hugoniot.py and the silica glass -> stishovite enthalpy of Akaogi & Navrotsky (1984),
# 10600 +- 680 cal/mol.
#
# Usage: python3 fs_stv_hugoniot.py [output csv]
#
# STS + Claude, Oct 2026

import os
import sys
here = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, here)
import numpy as np
from stv_model_hugoniot import (hugoniot, melt_crossing, load_melt, summary, linear_fit_hugoniot,
                                A_cases, A_base, dH_glass_stv, dH_glass_stv_err)

rho0 = 2.204   # g/cm3, reference density of the Table 2 fused silica fits
mc = load_melt()
out_csv = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '../data_silica_db/fs-stv-model-hugoniot.csv')

# fused silica solid-branch Hugoniot (Table 2): Us = a0 + a1 up, T(Us) quadratic; to the start of
# the temperature inversion (Us = 7.5 km/s, 65.8 GPa)
fs_a = [1.10707011, 1.6070127]
fs_T = [7684.97862763, -2052.46473422, 226.34570845]
fs_us_inversion = 7.5
fit = linear_fit_hugoniot(rho0, fs_a, fs_T, us_max=8.0, up_min=2.0)

tab = summary('fused silica', rho0, dH_glass_stv, dH_glass_stv_err, 34., [40., 60.], out_csv, fit, mc,
              rho_lo=4.0)
