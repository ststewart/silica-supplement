# omega_quartz_attenuation.py
#
# Empirical visible attenuation length of shocked alpha-quartz below the temperature
# inversion from the OMEGA decaying shock SOP data (Millot & Davies), using the same
# slab model as shot50_visible_model.py (Ye et al. 2024 approach), and comparison
# with the fused silica attenuation length from shot 50.
#
# Below the temperature inversion, the decaying-shock apparent temperatures are
# higher than the steady-shock Hugoniot. Here the shock-front temperature is fixed to
# the steady-shock Hugoniot fits and the attenuation length d(P) is fitted so that the
# slab model reproduces the decaying-shock apparent temperature.
#
# The slices are defined by the measured VISAR Us(t) extracted by omega_extract_time.py
# (one slice per VISAR sample, on the calibrated time axis). The fit uses the archived
# T vs Us records within the range of the measured Us(t). Run omega_extract_time.py and
# shot50_visible_model.py first.
#
# Assumptions:
#   - material shocked above the end of the temperature inversion (Us = 9.29 km/s;
#     opaque melt layer and liquid) is opaque, d = d_opaque
#   - below that, d depends on peak shock pressure (log-linear between two nodes)
#   - pressure behind the front is uniform and equal to the front pressure (isobaric);
#     this gives lower limits on the attenuation length. As a sensitivity test the fits
#     are repeated with the material behind the front released to a fraction
#     release_factors of the front pressure
#   - release temperatures from the stishovite Mie-Grueneisen model for all slices
#     (approximate for the liquid slices behind the transparent front)
#   - SOP is treated as monochromatic at T0 = 1.909 eV (649.5 nm; Gregor et al. 2016)
#   - shock front reflectivity below the opacity onset is negligible
#
# Usage: python3 omega_quartz_attenuation.py [opaque_us]
#   optional opaque_us (km/s) moves the opacity onset (e.g. 8.99 for the Hugoniot-solidus
#   intersection); outputs then get a separate prefix
#
# Outputs are written to output/
#
# STS + Claude, Oct 2026

import os
import sys
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from slab_emission import brightnessT, ReleaseTable, slab_radiance

##########################################################################
# user parameters

# all paths are relative to the folder containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
qtzdir = 'data/omega/'                              # archived T vs Us records
usfile = 'output/omega-extracted/omega-%s-Us-time.csv'  # from omega_extract_time.py
fsfile = 'output/shot50-visible-model.pkl'          # from shot50_visible_model.py
outdir = 'output/'
prefix = 'omega-qtz-'  # output file prefix

# shot: alpha-quartz thickness (micron) from data/omega/omega-thicknesses.png
shots = {'70813': 253.7, '70817': 278.0, '70820': 418.1}

wvln_sop = 1239.84/1.909*1.e-9  # m, T0 = 1.909 eV (Gregor et al. 2016)

d_opaque = 0.5        # micron, attenuation length of the opaque (melt layer and liquid) material
opaque_us = 9.29      # km/s, opaque above the end of the temperature inversion; None = solidus
fit_usmax = 9.6       # km/s, fit data below this shock velocity
dus_sub = 0.05        # km/s, subsample data for the fit (records are on a 0.01 km/s grid)
d_guess = 50.         # micron
release_factors = [0.95, 0.9, 0.85, 0.8, 0.75]  # sensitivity test: P behind front = factor x front pressure
release_factor = 1.0  # 1 = isobaric (baseline)

##########################################################################
# alpha-quartz Hugoniot fits (Table 2 of the main text)

rho0 = 2.65  # g/cm3
qtz_solid_usup = [1.96905488, 1.66507746]                  # solid branch, Us = a0 + a1 up
qtz_solid_TUs = [13817.14362212, -3554.25351558, 289.9608017]  # solid branch, T(Us) quadratic
qtz_liq_usup = [1.75601795e+00, 1.84956245e+00, -3.17458138e-02, 4.97145666e-04]  # liquid branch, Us(up) cubic
qtz_liq_TUs = [-4391.27584577, 1263.49808289, -81.26334251, 5.26247505]  # liquid branch, T(Us) cubic
qtz_wiggle_us1 = 9.29  # km/s, end of the solid branch (temperature inversion)
qtz_wiggle_us2 = 9.5   # km/s, start of the liquid branch
qtz_startsol_P = 100.50773755669458  # GPa, Hugoniot intersects the solidus

up_s = np.arange(2.5, 5., .001)
us_s = qtz_solid_usup[0] + qtz_solid_usup[1]*up_s
up_l = np.arange(3., 35., .005)
us_l = qtz_liq_usup[0] + qtz_liq_usup[1]*up_l + qtz_liq_usup[2]*up_l**2 + qtz_liq_usup[3]*up_l**3


def qtz_hugoniot(us):
    ''' up, P, rho, T on the alpha-quartz Hugoniot vs Us.
        Solid branch below qtz_wiggle_us1, liquid branch above qtz_wiggle_us2,
        linear in Us between them. '''
    us = np.asarray(us, dtype=float)
    ups = np.interp(us, us_s, up_s)
    Ts = qtz_solid_TUs[0] + qtz_solid_TUs[1]*us + qtz_solid_TUs[2]*us**2
    upl = np.interp(us, us_l, up_l)
    Tl = qtz_liq_TUs[0] + qtz_liq_TUs[1]*us + qtz_liq_TUs[2]*us**2 + qtz_liq_TUs[3]*us**3
    w = np.clip((us - qtz_wiggle_us1)/(qtz_wiggle_us2 - qtz_wiggle_us1), 0., 1.)
    up = (1 - w)*ups + w*upl
    T = (1 - w)*Ts + w*Tl
    P = rho0*us*up
    rho = rho0/(1. - up/us)
    return up, P, rho, T


usgrid = np.arange(6.5, 20., 0.001)
qtz_startsol_us = np.interp(qtz_startsol_P, qtz_hugoniot(usgrid)[1], usgrid)
print('alpha-quartz solidus intersection: %.1f GPa, Us = %.3f km/s' % (qtz_startsol_P, qtz_startsol_us))
if len(sys.argv) > 1:
    opaque_us = float(sys.argv[1])
    prefix = 'omega-qtz-opaque%.2f-' % opaque_us
opaque_P = qtz_startsol_P if opaque_us is None else qtz_hugoniot(opaque_us)[1]
qtz_hugoniot_us_opaque = np.interp(opaque_P, qtz_hugoniot(usgrid)[1], usgrid)
print('Opaque above %.1f GPa (Us = %.2f km/s)' % (opaque_P, qtz_hugoniot_us_opaque))

##########################################################################
# load data

data = {}
for s in shots:
    a = np.loadtxt(qtzdir + s + '/' + s + '_TandUs.txt', delimiter=',', skiprows=2)
    isort = np.argsort(a[:, 2])
    us = a[isort, 2]
    T = a[isort, 0]*1000.
    Terr = a[isort, 1]*1000.
    # subsample for the fit
    edges = np.arange(us.min(), fit_usmax + 1.e-6, dus_sub)
    usb, Tb, Teb = [], [], []
    for i in range(len(edges) - 1):
        ind = np.where((us >= edges[i]) & (us < edges[i+1]))[0]
        if len(ind):
            usb.append(us[ind].mean())
            Tb.append(T[ind].mean())
            Teb.append(Terr[ind].mean())
    data[s] = {'us': us, 'T': T, 'Terr': Terr, 'usfit': np.asarray(usb),
               'Tfit': np.asarray(Tb), 'Terrfit': np.asarray(Teb)}
    print('Shot %s: quartz %.1f micron, Us %.2f - %.2f km/s, %d fit points below %.1f km/s' %
          (s, shots[s], us.min(), us.max(), len(usb), fit_usmax))

Pmin_all = np.min([qtz_hugoniot(data[s]['us'].min())[1] for s in shots])
d_Pnodes = np.array([np.floor(Pmin_all), opaque_P])
print('Attenuation length nodes at P = %.1f and %.1f GPa' % tuple(d_Pnodes))

release = ReleaseTable(Pmin=20., Pmax=200., dP=2., Tmin=2500., Tmax=15000., dT=500.)

# slices from the measured VISAR Us(t): one slice per pixel
slices = {}
for s in shots:
    a = np.loadtxt(usfile % s, delimiter=',', skiprows=2)
    t, us = a[:, 0], a[:, 3]
    dh = us*np.gradient(t)  # Lagrangian thickness of each slice (micron); calibrated time axis
    up, P, rho, T = qtz_hugoniot(us)
    nn, jj = np.tril_indices(len(us))
    slices[s] = {'t': t, 'us': us, 'dh': dh, 'dx': dh*rho0/rho, 'P': P, 'T': T,
                 'x': np.cumsum(dh), 'opaque': P >= opaque_P, 'nn': nn, 'jj': jj}

##########################################################################
# model and fit


def atten_length(logd_nodes, P, opaque):
    if len(logd_nodes) == 1:
        d = np.exp(logd_nodes[0]) + 0.*P
    else:
        d = np.exp(np.interp(P, d_Pnodes, logd_nodes))
    return np.where(opaque, d_opaque, d)


def model_Tb(logd_nodes, sl):
    # material behind the front is released to release_factor x front pressure
    Prel = sl['P'][sl['nn']]*np.where(sl['jj'] < sl['nn'], release_factor, 1.)
    L, _ = slab_radiance(sl['dx'], sl['P'], sl['T'], atten_length(logd_nodes, sl['P'], sl['opaque']),
                         release, wvln_sop, sl['nn'], sl['jj'], Prel=Prel)
    return brightnessT(L, wvln_sop)


def residuals(logd_nodes, shotlist):
    res = []
    for s in shotlist:
        sl = slices[s]
        Tb = model_Tb(logd_nodes, sl)
        # model at the front for the last time each shock velocity is reached
        # (Us(t) from VISAR is not strictly monotonic)
        usmin = np.minimum.accumulate(sl['us'])
        last = np.where(np.diff(np.concatenate([usmin, [-1.]])) < 0)[0]
        inside = (data[s]['usfit'] >= usmin[-1]) & (data[s]['usfit'] <= usmin[0])
        Tmod = np.interp(data[s]['usfit'], np.flip(usmin[last]), np.flip(Tb[last]))
        res.append(np.where(inside, (Tmod - data[s]['Tfit'])/data[s]['Terrfit'], 0.))
    return np.concatenate(res)


def fit(shotlist, nnodes=None):
    ''' Fit log attenuation length at d_Pnodes; nnodes=1 fits one constant d below the opacity onset. '''
    p0 = np.zeros(nnodes or len(d_Pnodes)) + np.log(d_guess)
    sol = least_squares(residuals, p0, args=(shotlist,),
                        bounds=(np.log(0.1), np.log(1.e5)), x_scale=1.)
    dof = max(len(sol.fun) - len(sol.x), 1)
    chi2red = np.sum(sol.fun**2)/dof
    cov = np.linalg.pinv(sol.jac.T @ sol.jac)*max(chi2red, 1.)
    return sol.x, np.sqrt(np.diag(cov)), chi2red


for s in shots:
    sl = slices[s]
    ind = np.where(~sl['opaque'])[0]
    print('%s: VISAR record %.1f ns after shock entry, transparent region %.0f micron (%.1f-%.1f GPa), mean dUs/dx there %.4f km/s/micron' %
          (s, sl['t'][-1] + 0.5*(sl['t'][1] - sl['t'][0]), sl['dh'][ind].sum(), sl['P'][ind].min(),
           sl['P'][ind].max(), (sl['us'][ind].max() - sl['us'][ind].min())/sl['dh'][ind].sum()))

joint, joint_err, chi2 = fit(list(shots))
per = {s: fit([s]) for s in shots}
const, const_err, chi2c = fit(list(shots), nnodes=1)
per_const = {s: fit([s], nnodes=1) for s in shots}
print('\njoint fit (chi2_red = %.2f): d(%.0f GPa) = %.1f micron (x/ %.2f), d(%.0f GPa) = %.1f micron (x/ %.2f)' %
      (chi2, d_Pnodes[0], np.exp(joint[0]), np.exp(joint_err[0]), d_Pnodes[1], np.exp(joint[1]), np.exp(joint_err[1])))
print('constant d below opacity onset, joint (chi2_red = %.2f): d = %.1f micron (x/ %.2f)' %
      (chi2c, np.exp(const[0]), np.exp(const_err[0])))
for s in shots:
    p, pe, c = per_const[s]
    print('constant d, %s alone (chi2_red = %.2f): d = %.1f micron (x/ %.2f)' % (s, c, np.exp(p[0]), np.exp(min(pe[0], 50.))))

# sensitivity to the pressure behind the front (the isobaric case gives lower limits)
sens = [(1.0, joint, joint_err, const, const_err)]
for release_factor in release_factors:
    j_, je_, _ = fit(list(shots))
    c_, ce_, _ = fit(list(shots), nnodes=1)
    sens.append((release_factor, j_, je_, c_, ce_))
release_factor = 1.0
print('\nSensitivity to the pressure behind the front (joint fits):')
print('  P behind / P front   d(%.0f GPa)   d(%.0f GPa)   constant d   (micron)' % tuple(d_Pnodes))
for fac, j_, je_, c_, ce_ in sens:
    print('  %4.2f               %6.1f       %6.1f       %6.1f' % (fac, np.exp(j_[0]), np.exp(j_[1]), np.exp(c_[0])))

##########################################################################
# fused silica (shot 50) attenuation length for comparison

fs = None
if os.path.exists(fsfile):
    with open(fsfile, 'rb') as f:
        fs = pickle.load(f)
    nT = len(fs['T_Pnodes'])
    fs_logd_err = fs['logd_err_total'] if 'logd_err_total' in fs else np.sqrt(np.diag(fs['cov']))[nT:]
    for P in [56., 60., 64., 68.]:
        print('d(%.0f GPa): fused silica shot 50 = %.0f micron, alpha-quartz = %.0f micron' %
              (P, np.exp(np.interp(P, fs['d_Pnodes'], np.log(fs['d_nodes']))),
               np.exp(np.interp(P, d_Pnodes, joint))))
else:
    print('Run shot50_visible_model.py first to compare with fused silica.')

##########################################################################
# save

with open(outdir + prefix + 'release-pressure-sensitivity.csv', 'w') as f:
    f.write('P_behind_over_P_front,d_%.0fGPa_um,d_%.0fGPa_um,d_constant_um,err_factor_constant\n' % tuple(d_Pnodes))
    for fac, j_, je_, c_, ce_ in sens:
        f.write('%.2f,%.4g,%.4g,%.4g,%.3f\n' % (fac, np.exp(j_[0]), np.exp(j_[1]), np.exp(c_[0]), np.exp(ce_[0])))

with open(outdir + prefix + 'attenuation-fits.csv', 'w') as f:
    f.write('shot,d_form,chi2_red,d_%.0fGPa_um,d_%.0fGPa_um,err_factor_%.0fGPa,err_factor_%.0fGPa\n'
            % tuple(list(d_Pnodes) + list(d_Pnodes)))
    fits = [('joint', 'two-node', (joint, joint_err, chi2)),
            ('joint', 'constant', (const, const_err, chi2c))]
    fits += [(s, 'two-node', per[s]) for s in shots]
    fits += [(s, 'constant', per_const[s]) for s in shots]
    for s, form, (p, pe, c) in fits:
        p = np.resize(p, 2)
        pe = np.minimum(np.resize(pe, 2), 50.)
        f.write('%s,%s,%.3f,%.4g,%.4g,%.3g,%.3g\n' %
                (s, form, c, np.exp(p[0]), np.exp(p[1]), np.exp(pe[0]), np.exp(pe[1])))

##########################################################################
# figures

plt.rcParams.update({'font.size': 12})
qtz_col = '#2ca02c'  # alpha-quartz
fs_col = '#1845fb'   # fused silica

# apparent temperature vs Us for each shot
fig, axs = plt.subplots(1, len(shots), figsize=(15, 5), sharey=True)
ush = np.linspace(6.8, 11., 300)
for ax, s in zip(axs, shots):
    sl = slices[s]
    ax.plot(data[s]['us'], data[s]['T'], '.', color='grey', ms=2, label='OMEGA SOP (decaying)')
    ax.plot(ush, qtz_hugoniot(ush)[3], 'k-', lw=2, label='Steady-shock Hugoniot fit')
    ax.plot(sl['us'], model_Tb(joint, sl), '-', color=qtz_col, lw=2.5, label='Attenuation-length model, joint fit')
    ax.axvline(qtz_hugoniot_us_opaque, color='grey', ls=':', lw=1)
    ax.text(qtz_hugoniot_us_opaque - 0.05, 3000, 'transparent', ha='right', fontsize=9)
    ax.text(qtz_hugoniot_us_opaque + 0.05, 3000, 'opaque', fontsize=9)
    ax.set_xlim(6.8, 11.)
    ax.set_ylim(2500, 8000)
    ax.set_xlabel('Shock velocity (km/s)')
    ax.set_title('Shot %s (%.0f μm)' % (s, shots[s]))
axs[0].set_ylabel('Temperature (K)')
axs[0].legend(fontsize=9, loc='upper left')
fig.tight_layout()
fig.savefig(outdir + prefix + 'T-Us-model.png', dpi=200)
fig.savefig(outdir + prefix + 'T-Us-model.pdf')

# attenuation length comparison
fig, ax = plt.subplots(figsize=(7.5, 5.5))
Pq = np.linspace(d_Pnodes[0], d_Pnodes[1], 50)
lq = np.interp(Pq, d_Pnodes, joint)
eq = np.interp(Pq, d_Pnodes, joint_err)  # node uncertainties interpolated
ax.semilogy(Pq, np.exp(lq), '-', color=qtz_col, lw=2.5, label='α-quartz, OMEGA decaying shocks')
ax.fill_between(Pq, np.exp(lq - eq), np.exp(lq + eq), color=qtz_col, alpha=0.15)
if fs is not None:
    Pf = np.linspace(fs['d_Pnodes'][0], fs['d_Pnodes'][-1], 50)
    lf = np.interp(Pf, fs['d_Pnodes'], np.log(fs['d_nodes']))
    ef = np.interp(Pf, fs['d_Pnodes'], fs_logd_err)
    ax.semilogy(Pf, np.exp(lf), '-', color=fs_col, lw=2.5, label='Fused silica, shot 50 (vis + 4.8 μm)')
    ax.fill_between(Pf, np.exp(lf - ef), np.exp(lf + ef), color=fs_col, alpha=0.15)
# temperature inversions: fused silica 7.5-8.0 km/s and alpha-quartz 9.29-9.5 km/s
fs_inversion_P = [65.76, 73.54]  # GPa
qtz_inversion_P = [opaque_P, qtz_hugoniot(qtz_wiggle_us2)[1]]
for (p1, p2), c, name in [(fs_inversion_P, fs_col, 'fused silica'), (qtz_inversion_P, qtz_col, 'α-quartz')]:
    ax.axvspan(p1, p2, color=c, alpha=0.12, lw=0)
    ax.text(p2 + 0.8, 0.04, '%s\ninversion,\nopaque →' % name, color=c, fontsize=9,
            transform=ax.get_xaxis_transform())
ax.set_ylim(10., 1000.)
ax.set_xlim(30., 130.)
ax.set_xlabel('Peak shock pressure (GPa)')
ax.set_ylabel('Attenuation length 1/α, visible (μm)')
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(outdir + prefix + 'vs-fs-attenuation.png', dpi=200)
fig.savefig(outdir + prefix + 'vs-fs-attenuation.pdf')

print('\nOutputs written to', outdir)
