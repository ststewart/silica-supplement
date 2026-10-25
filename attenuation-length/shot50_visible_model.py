# shot50_visible_model.py
#
# Attenuation-length model for the visible (620-660 nm) thermal emission from the
# shot 50 decaying shock in fused silica, following the approach of
# Ye et al. 2024 (J. Appl. Phys. 136, 105904; supplement section S6).
#
# The shocked fused silica is treated as a stack of Lagrangian slices. Each
# slice is shocked to the Hugoniot state at the shock front when it is
# overrun, then releases isentropically. Each slice keeps an absorption
# coefficient set by its peak shock pressure (Ye et al. assumption). The
# radiance at the free surface at time t is the sum over slices of
#     B(T_slice) * (1 - exp(-delta_slice)) * exp(-sum of delta for the slices in front)
# where delta = slice thickness / attenuation length
# attenuated by the remaining unshocked fused silica.
#
# There is no measured attenuation depth for shocked fused silica, so the
# attenuation length d(P) (= 1/alpha) is fitted empirically, together with a
# smooth shock-front temperature T_H(P). The 4.8 micron back-out temperature
# (opaque-front model with absorption by the unshocked sample) is used as a joint
# constraint to break the T_H - d degeneracy.
#
# Assumptions (see notes printed at the end):
#   - shock velocity history is the piecewise model of the shot 50 analysis (main text)
#   - pressure behind the front is either uniform and equal to the front pressure
#     ('isobaric', default) or follows the pyKO simulation of shot 50 ('pyko'):
#         python3 shot50_visible_model.py pyko
#     In 'pyko' mode the pressure of each slice is the pyKO pressure at the same distance
#     behind the front, taken when the pyKO front reaches the same position, scaled by
#     the ratio of the model and pyKO front pressures. Outputs then start with shot50-pyko-.
#   - release temperatures from the stishovite Mie-Grueneisen model in
#     scripts/load_solid_models.py and scripts/define_eos.py
#   - gray body with emissivity EMISSIVITY of the optically thick limit
#
# Outputs are written to output/
#
# STS + Claude, Oct 2026

import os
import pickle
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from slab_emission import planck, brightnessT, blur, ReleaseTable, slab_radiance
from slab_emission import c2 as c2_planck

##########################################################################
# user parameters

# all paths are relative to the folder containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
visfile = 'data/shot50/S050_Lineout_timeprofile_620-660nm.txt'
calfile = 'data/shot50/S050_TemperatureProfile_-1000-1000ns_20nsWindow.txt'
irfile = 'data/shot50/S050_IR_4.8um_radiance.csv'  # 4.8 micron radiance from the IR pyrometer analysis
wvln_ir = 4.8e-6  # m
outdir = 'output/'
pykofile = 'data/pyko/shot50-pyko-aneos-fs.npz'  # pyKO simulation of shot 50, fused silica cells

# pressure behind the shock front: 'isobaric' or 'pyko'
release_pressure = sys.argv[1] if len(sys.argv) > 1 else 'isobaric'
prefix = 'shot50-pyko-' if release_pressure == 'pyko' else 'shot50-'

# streak time shift (shot 50 analysis)
tshift = -0.035  # microsec

# shock velocity history (shot 50 analysis)
tdecay = 0.45  # microsec
tkink = 0.64   # microsec
tend = 0.74    # microsec
usdecayrate = 10.  # km/s per microsec
P_IM = 68.06   # GPa impedance match pressure
thickness = 5.298  # mm sample thickness (shot record)

# visible emission
wvln_vis = 640.e-9  # m, center of 620-660 nm band
EMISSIVITY = 1.0    # gray body emissivity of the optically thick limit
alpha_unshocked_vis = 8.e-6  # 1/micron, unshocked fused silica at visible (shot 50 analysis)
sigma_inst = 3.0    # ns, Gaussian time response applied to the model
dt_model = 0.002    # microsec, slice time step
dt_bin_vis = 0.005  # microsec, binning of the visible lineout
fit_tmin = 0.0      # microsec, start of visible fit window
fit_tmax = tend - 0.005

# 4.8 micron constraint (shot 50 analysis)
USE_IR = True
alpha_unshocked_ir = 8.e-4  # 1/micron, unshocked fused silica at 4.8 micron (shot 50 analysis)
ir_tmin = 0.43   # microsec, back-out temperature used for t > ~tdecay
ir_tmax = tend - 0.01
dt_bin_ir = 0.010  # microsec
sigma_ir = 292.    # K, uncertainty of the steady-shock temperature

# fit parameterization
n_T_nodes = 9         # shock-front temperature nodes in pressure
n_d_nodes = 2         # attenuation length nodes in pressure (log-linear between nodes)
smooth_sigma = 50.    # K, allowed second difference between T nodes
d_guess = 100.        # micron, initial attenuation length

##########################################################################
# fused silica Hugoniot (scripts/define_eos.py; fits in silica-book/model-fits)
# fs_fit_stv_linear:        Us = a0 + a1 up
# fs_fitTUs_stv_quadratic:  T = a0 + a1 Us + a2 Us^2

rho0 = 2.204  # g/cm3
fs_usup = [1.10707011, 1.6070127]
fs_TUs = [7684.97862763, -2052.46473422, 226.34570845]

hug_up = np.arange(2.4, 4.5, .001)
hug_us = fs_usup[0] + fs_usup[1]*hug_up
hug_P = rho0*hug_us*hug_up
hug_rho = rho0/(1. - hug_up/hug_us)
hug_T = fs_TUs[0] + fs_TUs[1]*hug_us + fs_TUs[2]*hug_us**2


##########################################################################
# shock history on the model time grid

usstart = np.interp(P_IM, hug_P, hug_us)  # km/s
tgrid = np.arange(0., tend, dt_model)     # slice j is shocked between tgrid[j] and tgrid[j]+dt
tfront = tgrid + dt_model                  # time when slice j becomes the front slice
tmid = tgrid + 0.5*dt_model


def us_history(t):
    us = np.zeros(len(t)) + usstart
    ind = np.where(t > tdecay)[0]
    us[ind] = usstart - usdecayrate*(t[ind] - tdecay)
    ind = np.where(t > tkink)[0]
    us[ind] = usstart - usdecayrate*(tkink - tdecay)
    return us


slice_us = us_history(tmid)                         # km/s
slice_P = np.interp(slice_us, hug_us, hug_P)        # GPa peak (Hugoniot) pressure of slice
slice_up = np.interp(slice_us, hug_us, hug_up)
slice_rho = np.interp(slice_us, hug_us, hug_rho)
slice_THug = np.interp(slice_us, hug_us, hug_T)     # Hugoniot model temperature
slice_dh = slice_us*dt_model*1000.                  # micron, initial (Lagrangian) thickness
slice_dx = slice_dh*rho0/slice_rho                  # micron, compressed thickness
xfront = np.cumsum(slice_dh)                        # micron, shock position in initial coordinates
remaining = np.maximum(thickness*1000. - xfront, 0.)

nslice = len(tgrid)
# all pairs (front slice n, slice j behind it) with j <= n, ordered by n then j
nn, jj = np.tril_indices(nslice)

Pnode_lo = slice_P.min()
Pnode_hi = slice_P.max()
T_Pnodes = np.linspace(Pnode_lo, Pnode_hi, n_T_nodes)
d_Pnodes = np.linspace(Pnode_lo, Pnode_hi, n_d_nodes)

print('Shock travel distance in model (mm) = %.4f; sample thickness = %.4f mm' % (xfront[-1]/1000., thickness))
print('Us start = %.3f km/s, P range %.2f - %.2f GPa' % (usstart, Pnode_lo, Pnode_hi))



def pyko_release_pressure():
    ''' Pressure of slice j when slice n is at the front, from the pyKO simulation.
        The pyKO profile behind its front is sampled at the same distance behind the
        front when the pyKO front is at the model front position, and scaled by the
        ratio of model to pyKO front pressure. Within the numerically smeared pyKO
        front (4 cells) the pressure is the front pressure. '''
    d = np.load(pykofile)
    P = d['pres']/1.e9
    hp = (d['pos0'] - d['pos0'][0])*1.e6  # micron, Lagrangian cell position in the sample
    tp = d['time']*1.e6
    kf = np.array([np.where(P[i] > 20.)[0].max() if (P[i] > 20.).any() else -1 for i in range(len(tp))])
    iend = np.argmax(kf)  # breakout
    ok = np.where((kf[:iend + 1] >= 15))[0]
    xf = hp[kf[ok]]
    pf = np.array([P[i, kf[i] - 12:kf[i] - 4].mean() for i in ok])
    xf = np.maximum.accumulate(xf) + 1.e-6*np.arange(len(xf))  # monotonic for interpolation
    Prel = np.full((nslice, nslice), np.nan)
    hmid = xfront - 0.5*slice_dh
    for n in range(nslice):
        tn = np.interp(xfront[n], xf, tp[ok])
        i0 = np.clip(np.searchsorted(tp, tn) - 1, 0, len(tp) - 2)
        w = (tn - tp[i0])/(tp[i0 + 1] - tp[i0])
        prof = (1 - w)*P[i0] + w*P[i0 + 1]
        xfn = np.interp(tn, tp[ok], xf)
        pfn = np.interp(tn, tp[ok], pf)
        xbehind = hp[max(kf[i0] - 4, 0)]  # pyKO position behind the smeared front
        hpos = xfn - (xfront[n] - hmid[:n + 1])  # same distance behind the front in pyKO
        ratio = np.where(hpos > xbehind, 1., np.interp(hpos, hp, prof)/pfn)
        Prel[n, :n + 1] = np.clip(ratio, 0., 1.)*slice_P[n]
    return Prel


if release_pressure == 'pyko':
    Prel_matrix = pyko_release_pressure()
    Prel_pairs = Prel_matrix[nn, jj]
    release = ReleaseTable(Pmin=10.)
    print('Release pressure from pyKO: min %.1f GPa behind the front' % np.nanmin(Prel_pairs))
else:
    Prel_matrix = None
    Prel_pairs = None
    release = ReleaseTable()
Ttest = np.interp(Pnode_hi, hug_P, hug_T)
print('Release check: from %.1f GPa, %.0f K -> %.1f GPa: %.0f K' %
      (Pnode_hi, Ttest, Pnode_lo, release(np.array([Pnode_hi]), np.array([Ttest]), np.array([Pnode_lo]))[0]))

##########################################################################
# forward model


def front_T(TH_nodes, P):
    return np.interp(P, T_Pnodes, TH_nodes)


def atten_length(logd_nodes, P):
    if len(logd_nodes) == 1:
        return np.exp(logd_nodes[0]) + 0.*P
    return np.exp(np.interp(P, d_Pnodes, logd_nodes))


def forward(TH_nodes, logd_nodes, full=False):
    ''' Model visible radiance at the free surface vs tfront. '''
    TH = front_T(TH_nodes, slice_P)
    L, contrib = slab_radiance(slice_dx, slice_P, TH, atten_length(logd_nodes, slice_P),
                               release, wvln_vis, nn, jj, Prel=Prel_pairs)
    L *= EMISSIVITY*np.exp(-alpha_unshocked_vis*remaining)
    if full:
        return L, contrib
    return L


##########################################################################
# load data

streak = np.loadtxt(visfile, skiprows=2, delimiter=',')
streak[:, 0] += tshift
calstreak = np.loadtxt(calfile, skiprows=2, delimiter=',')
calstreak[:, 0] += tshift

# bin the visible lineout
edges = np.arange(fit_tmin, fit_tmax + 1.e-9, dt_bin_vis)
tvis, Lvis = [], []
for i in range(len(edges) - 1):
    ind = np.where((streak[:, 0] >= edges[i]) & (streak[:, 0] < edges[i+1]))[0]
    if len(ind) > 0:
        tvis.append(np.mean(streak[ind, 0]))
        Lvis.append(np.mean(streak[ind, 1]))
tvis = np.asarray(tvis)
Lvis = np.asarray(Lvis)
Tbvis = brightnessT(Lvis, wvln_vis)
ind = np.where((tvis > 0.1) & (tvis < 0.4))[0]
sigma_vis = np.std(Tbvis[ind])
print('Visible: steady Tb (0.1-0.4 us) = %.0f +- %.0f K (scatter of %.0f ns bins)' %
      (np.mean(Tbvis[ind]), sigma_vis, dt_bin_vis*1000))

# 4.8 micron back-out temperature: radiance corrected for the remaining unshocked sample
tir, Tir = np.array([]), np.array([])
if os.path.exists(irfile):
    ir = np.loadtxt(irfile, delimiter=',')
    t4 = ir[:, 0]*1.e6
    Dabs = np.interp(t4, tfront, remaining, left=thickness*1000., right=0.)
    Lbackout = ir[:, 1]/np.exp(-alpha_unshocked_ir*Dabs)
    Tbackout = brightnessT(Lbackout, wvln_ir)
    edges = np.arange(ir_tmin, ir_tmax + 1.e-9, dt_bin_ir)
    tir, Tir = [], []
    for i in range(len(edges) - 1):
        ind = np.where((t4 >= edges[i]) & (t4 < edges[i+1]))[0]
        tir.append(np.mean(t4[ind]))
        Tir.append(np.nanmean(Tbackout[ind]))
    tir, Tir = np.asarray(tir), np.asarray(Tir)
else:
    print('No IR file found; fitting visible only.')
    USE_IR = False
Pir = np.interp(tir, tmid, slice_P) if len(tir) else tir

##########################################################################
# fit

nT = n_T_nodes


def unpack(p):
    return p[:nT], p[nT:]


def residuals(p, use_ir):
    TH_nodes, logd_nodes = unpack(p)
    L = blur(forward(TH_nodes, logd_nodes), sigma_inst, dt_model)
    Lmod = np.interp(tvis, tfront, L, left=0.)
    Tbmod = brightnessT(np.maximum(Lmod, 1.e-30), wvln_vis)
    res = [(Tbmod - Tbvis)/sigma_vis]
    if use_ir:
        res.append((front_T(TH_nodes, Pir) - Tir)/sigma_ir)
    res.append(np.diff(TH_nodes, 2)/smooth_sigma)
    return np.concatenate(res)


def run_fit(use_ir):
    p0 = np.concatenate([np.interp(T_Pnodes, hug_P, hug_T) + 300., np.zeros(n_d_nodes) + np.log(d_guess)])
    lo = np.concatenate([np.zeros(nT) + 1500., np.zeros(n_d_nodes) + np.log(1.)])
    hi = np.concatenate([np.zeros(nT) + 8000., np.zeros(n_d_nodes) + np.log(1.e5)])
    sol = least_squares(residuals, p0, args=(use_ir,), bounds=(lo, hi), x_scale='jac')
    dof = max(len(sol.fun) - len(sol.x), 1)
    chi2red = np.sum(sol.fun**2)/dof
    J = sol.jac
    cov = np.linalg.pinv(J.T @ J)*max(chi2red, 1.)
    return sol, cov, chi2red


sol, cov, chi2red = run_fit(USE_IR)
TH_nodes, logd_nodes = unpack(sol.x)
perr = np.sqrt(np.diag(cov))
print('\nJoint fit (USE_IR=%s): reduced chi2 = %.2f' % (USE_IR, chi2red))
for i in range(nT):
    print('  T_H(%5.1f GPa) = %5.0f +- %4.0f K   (Hugoniot model %5.0f K)' %
          (T_Pnodes[i], TH_nodes[i], perr[i], np.interp(T_Pnodes[i], hug_P, hug_T)))
for i in range(n_d_nodes):
    d = np.exp(logd_nodes[i])
    print('  attenuation length d(%5.1f GPa) = %7.1f micron (x/ %.2f)' % (d_Pnodes[i], d, np.exp(perr[nT+i])))

# systematic uncertainty from the absolute calibrations: refit with the visible radiance
# scaled by its 1-sigma calibration factor and with the 4.8 micron temperatures shifted by
# sigma_ir; the two calibrations are independent and are added in quadrature with the
# statistical uncertainty
ind = np.where((calstreak[:, 0] > 0.1) & (calstreak[:, 0] < 0.4))[0]
f_cal = np.mean((c2_planck/(wvln_vis*calstreak[ind, 1]**2))*calstreak[ind, 2])
print('\nVisible calibration: %.0f K at %.0f K -> %.0f%% in 640 nm radiance' %
      (np.mean(calstreak[ind, 2]), np.mean(calstreak[ind, 1]), 100*f_cal))
Tbvis0, Tir0 = Tbvis.copy(), Tir.copy()
shifts = {}
for name, fac, dT in [('vis+', 1. + f_cal, 0.), ('vis-', 1. - f_cal, 0.), ('ir+', 1., sigma_ir), ('ir-', 1., -sigma_ir)]:
    Tbvis = brightnessT(Lvis*fac, wvln_vis)
    Tir = Tir0 + dT
    shifts[name] = run_fit(USE_IR)[0].x
Tbvis, Tir = Tbvis0, Tir0
sys_vis = 0.5*(np.abs(shifts['vis+'] - sol.x) + np.abs(shifts['vis-'] - sol.x))
sys_ir = 0.5*(np.abs(shifts['ir+'] - sol.x) + np.abs(shifts['ir-'] - sol.x))
perr_tot = np.sqrt(perr**2 + sys_vis**2 + sys_ir**2)
print('Uncertainties (1 sigma): statistical, visible calibration, 4.8 micron calibration, total')
for i in range(nT):
    print('  T_H(%5.1f GPa) = %5.0f K: %4.0f %4.0f %4.0f -> %4.0f K' %
          (T_Pnodes[i], TH_nodes[i], perr[i], sys_vis[i], sys_ir[i], perr_tot[i]))
for i in range(n_d_nodes):
    k = nT + i
    print('  d(%5.1f GPa) = %7.1f micron: x/ %.2f %.2f %.2f -> x/ %.2f' %
          (d_Pnodes[i], np.exp(logd_nodes[i]), np.exp(perr[k]), np.exp(sys_vis[k]), np.exp(sys_ir[k]),
           np.exp(perr_tot[k])))

# visible-only fit for comparison
if USE_IR:
    sol_v, cov_v, chi2red_v = run_fit(False)
    TH_nodes_v, logd_nodes_v = unpack(sol_v.x)
    print('\nVisible-only fit: reduced chi2 = %.2f' % chi2red_v)
    print('  T_H nodes:', np.round(TH_nodes_v).astype(int))
    print('  d nodes (micron):', np.round(np.exp(logd_nodes_v), 1))

##########################################################################
# results on the model grid

Lfit, contrib_fit = forward(TH_nodes, logd_nodes, full=True)
Lfit_b = blur(Lfit, sigma_inst, dt_model)
Tb_fit = brightnessT(Lfit_b, wvln_vis)
Tfront_fit = front_T(TH_nodes, slice_P)

# 1-sigma band on T_H(P): statistical (covariance of the nodes) and calibration terms
Pplot = np.linspace(Pnode_lo, Pnode_hi, 200)
W = np.array([np.interp(Pplot, T_Pnodes, np.eye(nT)[i]) for i in range(nT)]).T
TH_P = W @ TH_nodes
TH_P_stat = np.sqrt(np.einsum('ij,jk,ik->i', W, cov[:nT, :nT], W))
TH_P_err = np.sqrt(TH_P_stat**2 + (W @ sys_vis[:nT])**2 + (W @ sys_ir[:nT])**2)
TH_t_err = np.interp(slice_P, Pplot, TH_P_err)

# radiance-weighted mean depth of the emission behind the front
depth_mean = np.zeros(nslice)
for n in range(nslice):
    ind = np.where(nn == n)[0]
    xb = np.cumsum(slice_dx[:n+1][::-1])[::-1] - 0.5*slice_dx[:n+1]  # depth of slice j behind front
    w = contrib_fit[ind]
    depth_mean[n] = np.sum(w*xb)/np.sum(w) if np.sum(w) > 0 else np.nan

# apparent (uncorrected) visible Tb on the model grid
Tbvis_grid = np.interp(tfront, tvis, Tbvis)

##########################################################################
# save

os.makedirs(outdir, exist_ok=True)
header = ('time_us,P_GPa,Us_kms,up_kms,Tb_vis_obs_K,Tb_vis_model_K,T_front_corrected_K,'
          'T_front_err_K,T_hugoniot_model_K,atten_length_um,mean_emission_depth_um')
np.savetxt(outdir + prefix + 'visible-corrected-T.csv',
           np.column_stack([tfront, slice_P, slice_us, slice_up, Tbvis_grid, Tb_fit, Tfront_fit,
                            TH_t_err, slice_THug, atten_length(logd_nodes, slice_P), depth_mean]),
           delimiter=',', header=header, comments='', fmt='%.6g')
np.savetxt(outdir + prefix + 'visible-corrected-TP-nodes.csv',
           np.column_stack([T_Pnodes, TH_nodes, perr_tot[:nT], perr[:nT], sys_vis[:nT], sys_ir[:nT]]),
           delimiter=',', header='P_GPa,T_front_K,T_err_total_K,T_err_stat_K,T_err_vis_cal_K,T_err_ir_cal_K',
           comments='', fmt='%.6g')

with open(outdir + prefix + 'visible-model.pkl', 'wb') as f:
    pickle.dump({'time': tfront, 'P': slice_P, 'us': slice_us, 'up': slice_up,
                 'Tb_vis_obs': Tbvis_grid, 'Tb_vis_model': Tb_fit,
                 'T_front': Tfront_fit, 'T_front_err': TH_t_err,
                 'T_Pnodes': T_Pnodes, 'TH_nodes': TH_nodes, 'cov': cov,
                 'TH_err_total': perr_tot[:nT], 'logd_err_total': perr_tot[nT:],
                 'd_Pnodes': d_Pnodes, 'd_nodes': np.exp(logd_nodes),
                 'tir': tir, 'Tir': Tir, 'Pir': Pir, 'USE_IR': USE_IR,
                 'dh': slice_dh, 'dx': slice_dx, 'x': xfront, 'remaining': remaining,
                 'alpha_unshocked_vis': alpha_unshocked_vis, 'wvln': wvln_vis,
                 'release_pressure': release_pressure, 'Prel': Prel_matrix}, f)

##########################################################################
# figures

plt.rcParams.update({'font.size': 12})
panel_size = (4.8, 3.0)  # inches; shown at half the page width, fonts print at ~8 pt
ind_t = np.where(tfront <= tend)[0]

# radiance vs time
fig, ax = plt.subplots(figsize=panel_size)
ib = streak[:, 0] <= tend  # up to shock breakout
ax.plot(streak[ib, 0], streak[ib, 1]/1.e12, color='lightgrey', lw=0.8, label='Visible lineout (raw)')
ax.plot(tvis, Lvis/1.e12, 'k.', ms=4, label='Visible, %.0f ns bins' % (dt_bin_vis*1000))
ax.plot(tfront, Lfit_b/1.e12, 'r-', lw=2, label='Attenuation-length model')
ax.plot(tfront, EMISSIVITY*planck(Tfront_fit, wvln_vis)/1.e12, 'k--', label='Opaque shock front at fitted $T_H$')
ax.set_xlim(-0.1, 0.8)
ax.set_ylim(0, 1.3*np.nanmax(Lvis)/1.e12)
ax.set_xlabel('Time (μs)')
ax.set_ylabel('Radiance at 640 nm\n(10$^{12}$ W/(sr m$^3$))')
ax.legend(loc='lower left', bbox_to_anchor=(0.12, 0.02), fontsize=8)
fig.tight_layout()
fig.text(0.01, 0.99, 'A', fontsize=14, fontweight='bold', ha='left', va='top')  # panel letter
fig.savefig(outdir + prefix + 'vis-radiance-model.png', dpi=200)
fig.savefig(outdir + prefix + 'vis-radiance-model.pdf')

# temperature vs time
fig, ax = plt.subplots(figsize=panel_size)
ax.plot(tvis, Tbvis, 'k.', ms=4, label='Visible Tb (apparent)')
ib = calstreak[:, 0] <= tend  # up to shock breakout
ax.plot(calstreak[ib, 0], calstreak[ib, 1], color='grey', lw=1, label='Streak T file (20 ns window)')
ax.plot(tfront, Tb_fit, 'r-', lw=2, label='Model apparent Tb')
ax.plot(tfront, Tfront_fit, 'b-', lw=2, label='Corrected shock-front T')
ax.fill_between(tfront, Tfront_fit - TH_t_err, Tfront_fit + TH_t_err, color='b', alpha=0.15)
if len(tir):
    ax.plot(tir, Tir, '-', color='cornflowerblue')
    ax.fill_between(tir, Tir - sigma_ir, Tir + sigma_ir, color='cornflowerblue', alpha=0.25, label='4.8 μm')
for tt in [tdecay, tkink, tend]:
    ax.axvline(tt, color='grey', ls=':', lw=1)
ax.set_xlim(-0.05, 0.8)
ax.set_ylim(2500, 6000)
ax.set_xlabel('Time (μs)')
ax.set_ylabel('Temperature (K)')
ax.legend(loc='lower left', bbox_to_anchor=(0.08, 0.02), fontsize=8)
fig.tight_layout()
fig.text(0.01, 0.99, 'B', fontsize=14, fontweight='bold', ha='left', va='top')  # panel letter
fig.savefig(outdir + prefix + 'vis-T-time.png', dpi=200)
fig.savefig(outdir + prefix + 'vis-T-time.pdf')

# temperature vs pressure, colours as in Figure 7A of the main text
col_vis = 'darkgoldenrod'   # decaying shock, transparent front (visible)
col_ir = 'cornflowerblue'   # decaying shock, 4.8 micron
col_fs = '#1845fb'          # fused silica
fig, ax = plt.subplots(figsize=(8, 6))
ind = np.where((tfront > 0.025) & (tfront < tend))[0]
isort = ind[np.argsort(slice_P[ind], kind='stable')]
Tvis_err = np.interp(tfront, calstreak[:, 0], calstreak[:, 2])  # streak temperature uncertainty
ax.plot(slice_P[ind], Tbvis_grid[ind], '-', color=col_vis)
ax.fill_between(slice_P[isort], Tbvis_grid[isort] - Tvis_err[isort], Tbvis_grid[isort] + Tvis_err[isort],
                color=col_vis, alpha=0.35, label='Visible, apparent (uncorrected)')
if len(tir):
    ax.plot(Pir, Tir, '-', color=col_ir)
    ax.fill_between(Pir, Tir - sigma_ir, Tir + sigma_ir, color=col_ir, alpha=0.25, label='4.8 μm')
ax.plot(Pplot, TH_P, '-', color=col_fs, lw=2, label='Visible, corrected with attenuation-length model')
ax.fill_between(Pplot, TH_P - TH_P_err, TH_P + TH_P_err, color=col_fs, alpha=0.2)
ax.plot(hug_P, hug_T, 'k-', label='Fused silica Hugoniot fit (solid branch)')
ax.set_xlim(30, 75)
ax.set_ylim(2500, 6000)
ax.set_xlabel('Pressure (GPa)')
ax.set_ylabel('Temperature (K)')
ax.legend(fontsize=9, loc='upper left')
fig.tight_layout()
fig.savefig(outdir + prefix + 'vis-corrected-TP.png', dpi=200)
fig.savefig(outdir + prefix + 'vis-corrected-TP.pdf')

# attenuation length and emission depth
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].semilogy(Pplot, atten_length(logd_nodes, Pplot), 'k-', lw=2, label='Joint fit')
if USE_IR:
    ax[0].semilogy(Pplot, atten_length(logd_nodes_v, Pplot), 'k:', label='Visible-only fit')
ax[0].set_xlabel('Peak shock pressure (GPa)')
ax[0].set_ylabel('Attenuation length 1/α at 640 nm (μm)')
ax[0].legend(fontsize=9)
ax[1].plot(tfront[ind_t], depth_mean[ind_t], 'k-')
ax[1].set_xlabel('Time (μs)')
ax[1].set_ylabel('Radiance-weighted emission depth (μm)')
fig.tight_layout()
fig.savefig(outdir + prefix + 'vis-attenuation.png', dpi=200)
fig.savefig(outdir + prefix + 'vis-attenuation.pdf')

print('\nCorrected T at the end of decay (%.1f GPa): %.0f +- %.0f K; apparent visible Tb: %.0f K' %
      (Pnode_lo, TH_nodes[0], perr_tot[0], np.nanmean(Tbvis[(tvis > tkink) & (tvis < tend - 0.01)])))
print('Outputs written to', outdir)
