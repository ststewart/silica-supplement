# stv_model_hugoniot.py
#
# Model Hugoniot of a silica polymorph or glass transformed to stishovite, from the
# stishovite Mie-Grueneisen model of the main text (Driver et al. 2010 4000 K isotherm,
# gamma0 = 1.22 at rho_ref = 4.234 g/cm3 (Driver V0 = 23.558 A^3/SiO2), q = 2.22,
# cv = 3.65 x 3R/M), including the energy difference between the starting material
# and stishovite. Used by qtz_stv_hugoniot.py and fs_stv_hugoniot.py.
#
# Energy of the model stishovite relative to ambient stishovite (298 K, V_amb = Driver V0):
#   E(V,T) - E_stv,amb = A + [E4000(V) - E4000(V_amb)] + cv (T - 4000)
#   A = int_298^4000 cv_real dT at V_amb   (thermal energy to reach the 4000 K reference isotherm)
#   dE4000/dV = 4000 gamma cv / V - P4000  (energy along the 4000 K isotherm)
# Hugoniot from a starting material at rho0, 298 K, 0 GPa:
#   E - E_stv,amb = 0.5 P (V0 - V) - dE_tr,  dE_tr = E_stv,amb - E_start,amb (= dH at 0 GPa)
#
# STS + Claude, Oct 2026

import os
import sys
import numpy as np
from scipy.integrate import quad, cumulative_trapezoid

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from slab_emission import stv_pressure, stv_gamma, stv_cv, stv_rhoref, sio2mol

R = 8.314462618
M = sio2mol*1.e-3                 # kg/mol
cv_dp = 9.*R/M                    # Dulong-Petit, 3 atoms per SiO2 (J/kg/K)
rho_amb = stv_rhoref*1000.        # kg/m3, Driver V0 (300 K) = 23.558 A^3/SiO2
V_amb = 1./rho_amb
theta_D = 1109.                   # K, stishovite Debye temperature at 298 K (Akaogi et al. 2011)
cal = 4.184                       # J/cal

# transformation enthalpies to stishovite at 298 K, 0 GPa (Akaogi & Navrotsky 1984, Table IV,
# preferred values); dE = dH - P dV = dH at P = 0. J/kg.
dH_qtz_stv = 12400.*cal/M
dH_qtz_stv_err = 400.*cal/M
dH_glass_stv = 10600.*cal/M
dH_glass_stv_err = 680.*cal/M


def debye_frac(T, theta):
    ''' Debye heat capacity as a fraction of the Dulong-Petit value. '''
    x = theta/T
    return 3./x**3*quad(lambda y: y**4*np.exp(y)/(np.exp(y) - 1.)**2, 0., x)[0]


def A_debye(scale):
    ''' int_298^4000 cv dT with a Debye cv, scaled so that the high-T limit is scale x DP. '''
    return quad(lambda T: scale*cv_dp*debye_frac(T, theta_D), 298., 4000., limit=200)[0]


A_cases = {
    'Debye, high-T limit = model cv (1.22 DP)': A_debye(stv_cv/cv_dp),
    'Debye, high-T limit = Dulong-Petit':        A_debye(1.),
    'model cv constant from 298 K':              stv_cv*(4000. - 298.),
}
A_base = 'Debye, high-T limit = model cv (1.22 DP)'


def hugoniot(rho0, dE_tr, A, rho_lo=None):
    ''' Model Hugoniot from a material at rho0 (g/cm3). Returns P (GPa), rho (g/cm3), T (K),
        Us, up (km/s) on a density grid from rho_lo (default: Driver V0 density) to 6.5 g/cm3. '''
    V0 = 1./(rho0*1000.)
    rho = np.linspace(rho_amb if rho_lo is None else rho_lo*1000., 6500., 4000)  # kg/m3
    V = 1./rho
    g = stv_gamma(rho/1000.)
    P4 = stv_pressure(rho/1000., 4000.)*1.e9       # Pa
    dEdV = 4000.*g*stv_cv/V - P4
    E4 = cumulative_trapezoid(dEdV, V, initial=0.)
    E4 -= np.interp(rho_amb, rho, E4)              # E4000(V) - E4000(V_amb)
    rhs = 0.5*P4*(V0 - V) - dE_tr - A - E4
    tau = rhs/(stv_cv*(1. - 0.5*g*(V0 - V)/V))
    T = 4000. + tau
    P = P4 + g*stv_cv*tau/V
    up = np.sqrt(np.maximum(P*(V0 - V), 0.))
    us = V0*np.sqrt(np.maximum(P/(V0 - V), 0.))
    return P/1.e9, rho/1000., T, us/1000., up/1000.


def load_melt(path=os.path.join(here, '../silica-book/scripts/melt-curve.txt')):
    return np.loadtxt(path, delimiter=',', skiprows=1)


def melt_crossing(P, T, rho, mc):
    ''' First crossing of the model Hugoniot with the P-T melt curve. '''
    Tm = np.interp(P, mc[:, 0], mc[:, 1])
    k = np.where((T - Tm)[:-1]*(T - Tm)[1:] <= 0)[0][0]
    Ps = np.interp(0., (T - Tm)[k:k+2], P[k:k+2])
    return Ps, np.interp(Ps, P, T), np.interp(Ps, P, rho)


def linear_fit_hugoniot(rho0, a, Tc, us_max, up_min=1.5):
    ''' Hugoniot from Us = a0 + a1 up and T = sum Tc[i] Us^i, up to us_max. '''
    up = np.linspace(up_min, (us_max - a[0])/a[1], 2000)
    us = a[0] + a[1]*up
    T = sum(c*us**i for i, c in enumerate(Tc))
    return rho0*us*up, rho0*us/(us - up), T, us, up


def summary(label, rho0, dE, dE_err, P_from, P_report, out_csv, fit, mc, rho_lo=None):
    ''' Print the sensitivity table and write the baseline Hugoniot from P_from to the melt curve. '''
    print('%s: rho0 %.3f g/cm3; dE_tr %.4f +- %.4f MJ/kg (%.1f kJ/mol)' % (label, rho0, dE/1e6, dE_err/1e6, dE*M/1e3))
    print('rho_amb (Driver V0) = %.4f g/cm3; P4000(rho_amb) = %.1f GPa; cv = %.0f J/kg/K (DP %.0f)' %
          (rho_amb/1000., stv_pressure(np.array([rho_amb/1000.]), 4000.)[0], stv_cv, cv_dp))
    print('\nSensitivity (melt-curve crossing, and T at %s GPa):' % '/'.join('%.0f' % p for p in P_report))
    for name, A in A_cases.items():
        for lab, d in [('dE', dE), ('dE+1s', dE + dE_err), ('dE-1s', dE - dE_err)]:
            P, rho, T, us, up = hugoniot(rho0, d, A, rho_lo=rho_lo)
            Ps, Ts, rs = melt_crossing(P, T, rho, mc)
            print('  %-42s %-6s A=%.2f MJ/kg: melt at %5.1f GPa, %4.0f K, %.3f g/cm3; T = %s K' %
                  (name, lab, A/1e6, Ps, Ts, rs, ', '.join('%4.0f' % np.interp(p, P, T) for p in P_report)))
    P, rho, T, us, up = hugoniot(rho0, dE, A_cases[A_base], rho_lo=rho_lo)
    Ps, Ts, rs = melt_crossing(P, T, rho, mc)
    Pg = np.concatenate([np.arange(P_from, Ps, 1.), [Ps]])
    Pf, rhof, Tf, usf, upf = fit
    tab = np.column_stack([Pg, np.interp(Pg, P, rho), 1.e3/np.interp(Pg, P, rho), np.interp(Pg, P, T),
                           np.interp(Pg, P, us), np.interp(Pg, P, up),
                           np.interp(Pg, Pf, rhof), np.interp(Pg, Pf, Tf), np.interp(Pg, Pf, usf)])
    hdr = ('P_GPa,rho_gcc,V_cm3_per_kg,T_K,Us_kms,up_kms,rho_meas_gcc,T_meas_K,Us_meas_kms\n'
           '# %s (%.3f g/cm3, 298 K) -> stishovite model Hugoniot: Driver et al. 2010 4000 K isotherm,\n'
           '# gamma0 1.22 at rho_ref 4.234 g/cm3, q 2.22, cv 3.65x3R/M; dE(start->stv) %.1f kJ/mol (Akaogi & Navrotsky 1984);\n'
           '# thermal energy 298->4000 K at the Driver V0 from a Debye cv (theta 1109 K, Akaogi et al. 2011) scaled to the\n'
           '# model cv. Last row: crossing of the main-text melt curve. *_meas: Table 2 solid-branch fits at the same P\n'
           '# (extrapolated beyond the fit range).' % (label, rho0, dE*M/1e3))
    np.savetxt(out_csv, tab, delimiter=',', header=hdr, comments='', fmt='%.5g')
    print('\nBaseline model Hugoniot (Debye cv scaled to model cv):')
    print('   P(GPa)  rho(g/cc)  V(cm3/kg)   T(K)   Us(km/s)  up(km/s) | fit rho, T, Us')
    rows = list(range(0, len(tab) - 1, 5)) + [len(tab) - 1]
    for r in tab[rows].tolist():
        print('  %6.1f   %6.3f   %7.2f   %6.0f   %6.2f   %6.2f   | %6.3f %6.0f %6.2f' % tuple(r))
    print('Written to', out_csv)
    return tab
