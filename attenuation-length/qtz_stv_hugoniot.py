# qtz_stv_hugoniot.py
#
# Model Hugoniot of alpha-quartz (2.65 g/cm3) transformed to stishovite, from 70 GPa (single
# shock wave) to the melt curve, using the stishovite model in stv_model_hugoniot.py and the
# alpha-quartz -> stishovite enthalpy of Akaogi & Navrotsky (1984), 12400 +- 400 cal/mol.
#
# Usage: python3 qtz_stv_hugoniot.py [supplement folder] [melt curve] [output csv]
#   (the first two arguments are kept for compatibility; defaults are used if omitted)
#
# STS + Claude, Oct 2026

import os
import sys
here = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, here)
import numpy as np
from stv_model_hugoniot import (hugoniot as _hugoniot, melt_crossing as _melt_crossing, load_melt,
                                summary, linear_fit_hugoniot, A_cases, A_base, M,
                                dH_qtz_stv as dH_an84, dH_qtz_stv_err as dH_an84_err)

rho0 = 2.65
mc = load_melt(sys.argv[2]) if len(sys.argv) > 2 else load_melt()
out_csv = sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, '../data_silica_db/qtz-stv-model-hugoniot.csv')
dH_akaogi95 = (3.40 + 33.62)*1.e3/M   # J/kg: Akaogi et al. 1995 (qtz->coe + coe->stv), sensitivity only

# measured alpha-quartz solid-branch Hugoniot (Table 2): Us = a0 + a1 up, T(Us) quadratic
qa = [1.96905488, 1.66507746]
qT = [13817.14362212, -3554.25351558, 289.9608017]
fit = linear_fit_hugoniot(rho0, qa, qT, us_max=10., up_min=2.5)


def hugoniot(A, dE_tr):
    return _hugoniot(rho0, dE_tr, A)


def melt_crossing(P, T, rho):
    return _melt_crossing(P, T, rho, mc)


def measured(P):
    ''' Us, rho, T of the Table 2 fit at pressure P. '''
    Pf, rhof, Tf, usf, upf = fit
    return np.interp(P, Pf, usf), np.interp(P, Pf, rhof), np.interp(P, Pf, Tf)


tab = summary('alpha-quartz', rho0, dH_an84, dH_an84_err, 70., [70., 90.], out_csv, fit, mc)
P, rho, T, us, up = hugoniot(A_cases[A_base], dH_akaogi95)
print('With the Akaogi et al. 1995 enthalpy (37.0 kJ/mol): melt at %.1f GPa, %.0f K' % melt_crossing(P, T, rho)[:2])
