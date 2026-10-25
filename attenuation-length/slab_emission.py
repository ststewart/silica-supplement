# slab_emission.py
#
# Shared functions for the attenuation-length (slab) model of thermal emission from
# decaying shocks, following Ye et al. 2024 (J. Appl. Phys. 136, 105904; supplement S6).
#
# Used by shot50_visible_model.py and omega_quartz_attenuation.py
#
# STS + Claude, Oct 2026

import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Planck constants as used in shot50-decay-Copy1.ipynb
c1 = 1.19118e-16  # W m^2 / sr  (2hc^2)
c2 = 1.4388e-2    # m K


def planck(T, wvln):
    ''' Spectral radiance W/(sr m^3). '''
    return c1/((wvln**5)*(np.exp(c2/(wvln*T)) - 1.))


def brightnessT(L, wvln):
    L = np.asarray(L, dtype=float)
    out = np.full(L.shape, np.nan)
    ind = L > 0
    out[ind] = (c2/wvln)/np.log(c1/(L[ind]*wvln**5) + 1.)
    return out


def blur(y, sigma_ns, dt_us):
    ''' Gaussian time response; y is zero before the first point. '''
    if sigma_ns <= 0:
        return y
    sig = sigma_ns/1000./dt_us
    k = np.arange(-int(4*sig) - 1, int(4*sig) + 2)
    kern = np.exp(-0.5*(k/sig)**2)
    kern /= kern.sum()
    pad = len(k)
    ypad = np.concatenate([np.zeros(pad), y, np.zeros(pad) + y[-1]])
    return np.convolve(ypad, kern, mode='same')[pad:-pad]


##########################################################################
# stishovite Mie-Grueneisen model (scripts/load_solid_models.py, define_eos.py)
# P(rho,T) = P4000(rho) + gamma rho cv (T - 4000)
# gamma = g0 (rhoref/rho)^q

sio2mol = 28.085 + 2*15.999
stv_g0 = 1.22
stv_q = 2.22
stv_cvfactor = 3.65
stv_cv = stv_cvfactor*3*8.314/sio2mol/1.e-3  # J/kg/K
v0dft = 159  # bohr^3/SiO2
stv_rhoref = 1./(v0dft/sio2mol*np.power(0.529177249, 3)/1.6605391/1000.)/1.e3  # g/cm3

# Driver et al. 2010 4000 K isotherm (silica data_silica_db/stvcacl2.txt)
stvcacl2 = np.array([[150, 122.55], [143.28358208955223, 123.75], [129.40298507462686, 126.3],
                     [116.86567164179104, 128.7], [105, 131.25], [101.64179104477611, 132],
                     [77.91044776119404, 137.85], [59.99999999999999, 143.25],
                     [48.80597014925373, 147.14999999999998], [34.92537313432836, 152.85],
                     [23.731343283582085, 158.55], [20.14925373134328, 160.46977611940298],
                     [16.119402985074622, 162.75], [10.746268656716417, 166.35],
                     [7.164179104477611, 168.75], [5.8208955223880565, 169.95],
                     [0, 175], [-5, 180], [-10, 185], [-15, 190]])
iso_P = stvcacl2[:, 0]
iso_rho = sio2mol/(stvcacl2[:, 1]/np.power(1.8897161646321, 3))*1.6605391


def stv_gamma(rho):
    return stv_g0*np.power(stv_rhoref/rho, stv_q)


def stv_pressure(rho, T):
    ''' Pressure (GPa) for density rho (g/cm3) and temperature T (K).
        The 4000 K isotherm is extrapolated linearly above 150 GPa. '''
    rr = np.flip(iso_rho)
    pp = np.flip(iso_P)
    p4000 = np.interp(rho, rr, pp)
    slope = (pp[-1] - pp[-2])/(rr[-1] - rr[-2])
    ind = rho > rr[-1]
    p4000 = np.where(ind, pp[-1] + slope*(rho - rr[-1]), p4000)
    return p4000 + stv_cv*(T - 4000.)*stv_gamma(rho)*rho*1000./1.e9


def stv_isentrope_T(rho, rhostart, Tstart):
    ''' Temperature along the isentrope through (rhostart, Tstart); dlnT/dlnrho = gamma. '''
    return Tstart*np.exp(-(stv_g0/stv_q)*(np.power(stv_rhoref/rho, stv_q) - np.power(stv_rhoref/rhostart, stv_q)))


class ReleaseTable:
    ''' T_release/T_H on a (P_H, T_H, P) grid for isentropic release from P_H to P.
        Queries outside the grid are clipped to the grid bounds. '''

    def __init__(self, Pmin=30., Pmax=76., dP=1., Tmin=1500., Tmax=8500., dT=250.):
        rhoarr = np.linspace(3.6, 7.0, 3000)
        self.PHg = np.arange(Pmin, Pmax + 0.1*dP, dP)
        self.THg = np.arange(Tmin, Tmax + 0.1*dT, dT)
        self.Pg = self.PHg
        ratio = np.ones((len(self.PHg), len(self.THg), len(self.Pg)))
        for i in range(len(self.PHg)):
            for j in range(len(self.THg)):
                piso = stv_pressure(rhoarr, self.THg[j])
                rhoH = np.interp(self.PHg[i], piso, rhoarr)
                tarr = stv_isentrope_T(rhoarr, rhoH, self.THg[j])
                parr = stv_pressure(rhoarr, tarr)
                ratio[i, j, :] = np.interp(self.Pg, parr, tarr)/self.THg[j]
        # no heating on "release" to a higher pressure
        for i in range(len(self.PHg)):
            ratio[i, :, self.Pg > self.PHg[i]] = 1.
        self.interp = RegularGridInterpolator((self.PHg, self.THg, self.Pg), ratio)

    def __call__(self, PH, TH, P):
        pts = np.column_stack([np.clip(PH, self.PHg[0], self.PHg[-1]),
                               np.clip(TH, self.THg[0], self.THg[-1]),
                               np.clip(P, self.Pg[0], self.Pg[-1])])
        return TH*self.interp(pts)


##########################################################################
# slab radiative transfer


def slab_radiance(dx, P_H, T_H, d, release, wvln, nn, jj, Prel=None, Trel=None):
    ''' Radiance at the shock front for each front slice n.

        dx      compressed thickness of each slice (micron), in order of shocking
        P_H,T_H peak (Hugoniot) state of each slice
        d       attenuation length of each slice (micron), set by its peak pressure
        release ReleaseTable; slices release isentropically to the current front pressure
        nn, jj  index pairs (front slice n, slice j behind it) with j <= n
        Prel    optional pressure of slice j when slice n is at the front, for each (n, j) pair;
                default is the front pressure P_H[nn] (isobaric)
        Trel    optional temperature of slice j when slice n is at the front, for each (n, j)
                pair; overrides the isentropic release temperature

        Returns L (W/(sr m^3)) for each front slice and the contribution of each (n, j) pair.
    '''
    delta = dx/d  # slice thickness in units of the attenuation length
    cum = np.cumsum(delta)
    tau = cum[nn] - cum[jj]  # attenuation by the slices between slice j and the front
    if Trel is None:
        Trel = release(P_H[jj], T_H[jj], P_H[nn] if Prel is None else Prel)
    contrib = planck(Trel, wvln)*(1. - np.exp(-delta[jj]))*np.exp(-tau)
    L = np.bincount(nn, weights=contrib, minlength=len(dx))
    return L, contrib
