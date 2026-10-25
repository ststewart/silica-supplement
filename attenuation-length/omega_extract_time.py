# omega_extract_time.py
#
# Shock velocity and temperature vs time for the OMEGA alpha-quartz decaying shocks
# from the time-resolved VISAR and SOP analysis of M. Millot (data/omega/<shot>/):
#   <shot>_Us1.txt   VISAR channel 1 shock velocity (km/s) vs time (ns)
#   <shot>_Us2.txt   VISAR channel 2 shock velocity (km/s) vs time (ns)
#   <shot>_T.txt     SOP temperature (eV) vs time (ns)
# The time axes include the streak-camera sweep calibration.
#
# The shock velocity used in the analysis is VISAR channel 1 from shock entry into the
# quartz (first point above 5 km/s) to the end of the reliable record. Isolated spikes are
# replaced by the running median, and the record is ended where Us turns upward at the
# end of the record. The VISAR record ends before breakout because the shock front loses
# its reflectivity. Breakout is taken from channel 2, which records the jump to the
# free-surface velocity; it is placed on the channel 1 clock using the shock-entry times
# of the two channels.
#
# Outputs (output/omega-extracted/):
#   omega-<shot>-Us-time.csv      time after shock entry (ns), Lagrangian shock position (micron),
#                                 Us (km/s), Us 9-point running mean (km/s)
#   omega-<shot>-T-time.csv       time after shock entry (ns), SOP temperature (K), from entry to breakout
#   omega-time-base-summary.csv   shock entry, end of VISAR record and breakout times for each shot
#   omega-Us-T-time.png/pdf       shock velocity and temperature vs time, with the breakout time
#
# STS + Claude, Oct 2026

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# all paths are relative to the folder containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))
rawdir = 'data/omega/'
outdir = 'output/omega-extracted/'

# shot: alpha-quartz thickness (micron) from data/omega/omega-thicknesses.png
shots = {'70813': 253.7, '70817': 278.0, '70820': 418.1}

eV = 11604.5          # K per eV
us_entry = 5.         # km/s, shock entry into the quartz
nsmooth = 9           # points, running mean for the smoothed Us
spike = 0.4           # km/s, isolated points further than this from the running median are replaced
cut_rise = 0.3        # km/s, end the record where Us rises this far above its running minimum
qtz_col = '#2ca02c'   # alpha-quartz, as in the main text figures


def load_ty(f):
    ''' Two-column tab-separated file with a header; blank values are returned as NaN. '''
    lines = open(f).read().replace('\r\n', '\n').replace('\r', '\n').split('\n')[1:]
    out = []
    for l in lines:
        p = l.split('\t')
        if len(p) < 2 or p[0].strip() == '':
            continue
        try:
            v = float(p[1])
        except ValueError:
            v = np.nan
        out.append((float(p[0]), v))
    return np.array(out)


def running_mean(y, n):
    ypad = np.concatenate([np.zeros(n//2) + y[0], y, np.zeros(n//2) + y[-1]])
    return np.convolve(ypad, np.ones(n)/n, mode='valid')


os.makedirs(outdir, exist_ok=True)
summary = []
res = {}
for s, H in shots.items():
    U1 = load_ty(rawdir + s + '/' + s + '_Us1.txt')
    U2 = load_ty(rawdir + s + '/' + s + '_Us2.txt')
    T = load_ty(rawdir + s + '/' + s + '_T.txt')

    # VISAR channel 1 from shock entry
    i0 = np.argmax(U1[:, 1] > us_entry)
    t0 = U1[i0, 0]
    t = U1[i0:, 0] - t0
    us = U1[i0:, 1].copy()
    med = np.array([np.median(us[max(0, k - 5):k + 6]) for k in range(len(us))])
    bad = np.abs(us - med) > spike
    us[bad] = med[bad]
    sm = running_mean(us, nsmooth)
    runmin = np.minimum.accumulate(sm)
    rise = np.where((sm - runmin > cut_rise) & (np.arange(len(us)) > 0.8*len(us)))[0]
    iend = np.argmin(sm[:rise[0]]) if len(rise) else len(us) - 1
    t, us, sm = t[:iend + 1], us[:iend + 1], sm[:iend + 1]
    x = np.concatenate([[0.], np.cumsum(0.5*(us[1:] + us[:-1])*np.diff(t))])

    # breakout from VISAR channel 2: the last time Us2 is away from its final (free-surface) level
    t2, u2 = U2[:, 0], U2[:, 1]
    t0_2 = t2[np.argmax(u2 > us_entry)]
    ufs = np.nanmedian(u2[-50:])
    tb = t2[np.where(np.abs(u2 - ufs) > 1.)[0][-1] + 1] - t0_2  # after shock entry

    # SOP temperature from shock entry to breakout
    ok = np.isfinite(T[:, 1]) & (T[:, 0] - t0 >= 0.) & (T[:, 0] - t0 <= tb)
    tT = T[ok, 0] - t0
    TK = T[ok, 1]*eV

    print('%s: shock entry at %.3f ns; VISAR record ends %.2f ns after entry at Us = %.2f km/s, x = %.0f um; '
          'breakout %.2f ns after entry (mean Us %.2f km/s over the remaining %.0f um)' %
          (s, t0, t[-1], sm[-1], x[-1], tb, (H - x[-1])/(tb - t[-1]), H - x[-1]))
    summary.append([s, H, t0, t[-1], sm[-1], x[-1], tb])

    np.savetxt(outdir + 'omega-%s-Us-time.csv' % s, np.column_stack([t, x, us, sm]), delimiter=',',
               header=('time_ns,x_um,Us_kms,Us_smooth_kms\n'
                       '# shot %s; VISAR channel 1 (M. Millot); t = 0 at shock entry into the quartz; '
                       'x = Lagrangian shock position' % s),
               comments='', fmt='%.6g')
    np.savetxt(outdir + 'omega-%s-T-time.csv' % s, np.column_stack([tT, TK]), delimiter=',',
               header=('time_ns,T_K\n# shot %s; SOP temperature (M. Millot); t = 0 at shock entry into '
                       'the quartz (VISAR channel 1 clock); shock entry to breakout' % s),
               comments='', fmt='%.6g')
    res[s] = {'t': t, 'us': us, 'tT': tT, 'T': TK, 'tb': tb}

with open(outdir + 'omega-time-base-summary.csv', 'w') as f:
    f.write('shot,thickness_um,shock_entry_ns,visar_end_after_entry_ns,Us_at_visar_end_kms,'
            'x_at_visar_end_um,breakout_after_entry_ns\n')
    for r in summary:
        f.write('%s,%.1f,%.3f,%.3f,%.3f,%.1f,%.3f\n' % tuple(r))

##########################################################################
# figures: data used in the analysis, with the breakout time

plt.rcParams.update({'font.size': 12})
fig, axs = plt.subplots(1, len(shots), figsize=(16, 4.8))
for ax, s in zip(axs, shots):
    r = res[s]
    ax.plot(r['t'], r['us'], '-', color='k', lw=1.5)
    ax.set_ylim(6, 20)
    ax.set_xlim(-1, r['tb'] + 2)
    ax.set_xlabel('Time after shock entry (ns)')
    ax.set_ylabel('Shock velocity (km/s)')
    ax.set_title('Shot %s (%.0f μm quartz)' % (s, shots[s]))
    axT = ax.twinx()
    axT.plot(r['tT'], r['T'], '-', color=qtz_col, lw=1)
    axT.set_ylim(2000, 25000)
    axT.set_ylabel('Temperature (K)', color=qtz_col)
    axT.tick_params(axis='y', colors=qtz_col)
    ax.axvline(r['tb'], color='grey', ls='--', lw=1)
fig.tight_layout()
fig.savefig(outdir + 'omega-Us-T-time.png', dpi=200)
fig.savefig(outdir + 'omega-Us-T-time.pdf')

print('Outputs written to', outdir)
