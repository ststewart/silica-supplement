# This script defines key variables used in the EOS models and plots
#
# Reference points and values are defined in
#   dft-offset.ipynb
#   melt-curve-fits.ipynb
#   Hugoniot fit notebooks for fs, qtz, cst
#
# Updated 10/20/2025 sts

# define constants
Rgas    = 8.314           # J/K/mol
sio2mol = 28.085+2*15.999 # g per formula unit

##########################################################
# Adjustments to the ab initio reference models

# In dft-offset.ipynb
# Liquid silica DFT model with LDA function has a density offset
# Shift to match fused silica and alpha-quartz Hugoniots
# This shift is applied to the Karki et al 2007 liquid PVT EOS model
# LDA is more dense than data
lda_liquid_density_shift =  -0.1460732034660399  # g/cm3
lda_liquid_density_shift_std =  0.04566345115307819  # g/cm3

# In additional-plots.ipynb:
# FS, qtz, cst Hugoniot intersections with MC
# This shift is applied to melt curve liquidus from G24 and G16
# Shift liquidus densities from Geng and Mohn 2024 shift to match
# new with MC errors and actual GM24 points
# PBE is under-dense compared to data
pbe_liquid_density_shift =  0.1348471632115231  # g/cm3
pbe_liquid_density_shift_std =  0.0271945065487611  # g/cm3

# The stv and sft densities are taken from QMC simulations
# in Driver et al. 2010
# These calculations appear to be in good agreement with the shock data
# no shift required
stv_density_shift = 0.0 # g/cm3
sft_density_shift = -0.02 # g/cm3 # to align with Park et al. 2025 no density shift on melt curve

#############################################################################
# TRIPLE POINTS AND SOLID ISOTHERM TIE POINTS
#
# Coesite to stishovite on the melt curve
# estimate from extrapolated T04 GGA boundary; same value as 
# Belonoshko & Drubovinsky 1995
#Pt1 = 14.5  # GPa: coesite to stishovite on the melt curve
#Tt1 = 3100  # K
# Andrault et al. 2022
#Pt1 = 15  # GPa: coesite to stishovite on the melt curve
#Tt1 = 2900  # K
# 
# Defined in melt-curve-fits.ipynb
# using fitted melt curve at pressure from Zhang et al. JGR 1993
# Cst-Stv TP:
Pt1 = 13.7 # GPa - assumed
Tt1 =  3220.7975619137674  # K

# 2nd order LDM to HDM - THIS IS NOT AN ABRUPT TRANSITION IN THE LIQUID
# Andrault 40 GPa, 4800 K
# Calculated in melt-curve-fits.ipynb
# LDM-HDM Transition estimated near the change in slope of entropy on the melt curve
Ptm = 50 # GPa - estimated from dS/dP change on melt curve
Ttm =  4841.313561513845  # K

# stishovite to CaCl2-type on the melt curve
# Andrault et al. 2022
#Pt2 = 85   # GPa: stishovite to CaCl2-type on the melt curve
#Tt2 = 5200  # K
# In melt-curve-fits.ipynb
# Using GGA-PBE stishovite to CaCl2-type boundary intersection with fiducial melt curve
# Stv-CaCl2 TP:
Pt2 =  88.91220821388308  # GPa
Tt2 =  5582.379237642232  # K

# melt-curve-fits.ipynb
# use extrapolation of GGA phase boundary to melt curve
# CaCl2-sft TP:
Pt3 =  141.47365401777512  # GPa
Tt3 =  6077.910180343331  # K

# upper limit for sft calculation sft->pyrite -- this is a 5% volume change!
# https://www.science.org/doi/10.1126/science.1114879
# pyrite is a 6 close Si-O and 2 far Si-O coordinated structure
# the sft-pyrite-liquid triple point is not known
# check melt-curve.ipynb for clean solidus and liquidus
Pt4 = 220 # GPa # just a limit for the sft calculation in load_solid_models.py
# do not have a pyrite model Park et al. 25 TP4 is at 165 GPa for PBEsol calcs

# solid-solid tie points from Tsuchiya et al. 2004 GGA
# stv to CaCl2
#Pt12shift = 0 # remove frm code
# in melt-curve-fits.ipynb
Pt12 = 62.8270 # GPa 
Tt12 = 4000 # K
# CaCl2 to seifertite
#Pbbshift = 0 # remove frm code
# in melt-curve-fits.ipynb
Pt23 =  113.4110155408729 # GPa
Tt23 = 4000 # K

#############################################################################
# 0K/300 K transition pressures
# Prefer value slightly lower than GGA; LDA is too low transition pressure
Pt2togo = 53.3  # GPa: stishovite to CaCl2-type from Togo et al. 2008
Pt2lda  = 44    # GPa: CaCl2-type to seifertite LDA DFT Oganov et al. 2005
Pt3lda  = 84    # GPa: CaCl2-type to seifertite LDA DFT Oganov et al. 2005
Pt4lda  = 204   # GPa: seifertite to pyrite-type LDA DFT Oganov et al. 2005
Pt2gga  = 55    # GPa: stishovite to CaCl2-type GGA DFT Oganov et al. 2005
Pt3gga  = 97    # GPa: CaCl2-type to seifertite GGA DFT Oganov et al. 2005
Pt4gga  = 217   # GPa: seifertite to pyrite-type GGA DFT Oganov et al. 2005
Pt2gga  = 56.1    # GPa: stishovite to CaCl2-type GGA DFT Tsuchiya et al. 2004
Pt3gga  = 106.3   # GPa: CaCl2-type to seifertite GGA DFT Tsuchiya et al. 2004
Pt4P25  = 198 # GPa sft-pyrite PBEsol Park et al. 2020

# preferred values
# Using GGA transition pressures for self-consistency with melt curve PBE
# using 300 K reference densities from Driver et al. 2010 
T_0k   = 300 # temperature corresponding to reference density
Pt1_0k = 6.3 # from where?
Pt2_0k = Pt2gga
Pt3_0k = Pt3gga # Togo 2008 estimates transition pressure at 96.2 GPa; similar to GGA value of 97
Pt4_0k = Pt4P25

# Define 0 K isotherms for the solid phases for pressure-volume plots
# these are for plotting purposes and not used in the melt curve solution
# [V0, K0, dK/dP] in [A^3, GPa, -]
# Oganov are slightly shifted from Driver et al. 2010; and transition pressures are too low with LDA
#stv   = [23.1559, 318.326, 4.3726] # Oganov et al. 2005 Table II
#cacl2 = [23.5751, 258.354, 4.6135] # Oganov et al. 2005 Table II transition 44 GPa
#pbo2  = [22.7805, 324.426, 4.2327] # Oganov et al. 2005 Table II transition 84 GPa
# use Driver values for 0 K Vinet equation; however density is 300 K value
stv   = [23.55848, 305, 3.7] # Driver et al. 2010 Table 1
cacl2 = [23.55848, 305, 3.7] # Driver et al. 2010 Table 1
pbo2  = [22.9409,  329, 4.1] # Driver et al. 2010 Table 1
pyrite  = [12.97,  393, 3.83] # Das et al. 2020 Table 2
#
#########################################################################
# ENTROPY PLOTS NOT COMPLETE AT THIS TIME
# Entropy tie points
# will be completed in future work

# Kraus et al. 2012 liquid reference point
krausrefpt = [1/0.2153,70,4500,3820] # (rho,p,t,s) # Kraus 2012 orig
#
# Liquid reference point
# 70 GPa and 4500 K were chosen to match Kraus et al. JGR 2012
# liquid reference point used to calculate isentropes
# match Connolly JGR 2016 vapor curve boiling point at 1 bar: 
#   T=3130$\pm$30 K and S=3560$\pm$100 J/K/kg.
# density of the reference point rho = 4.645 g/cm3
# here is close to Connolly 1 bar boiling point with Karki isotherms
# this is the Kraus 2012 liquid reference point that uses a too small dS of melting 
# at 70 GPa and 4500 K
#liquid_reference_point = [-1,70,4500,3820] # (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg) # 

# this work
# empirically adjusted to match tie points
liquid_reference_point = [-1,70,4500,3980] # (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg) # 
# This value empirically meets the 1 bar boiling point and Kraus liquidus crossing qtz Hugoniot entropy.
liquid_cvfactor = 4.5 # prefactor to nR for constant heat capacity across model PT range
#
# stishovite/cacl2 reference point
# the density of the tie point is defined by the stishovite EOS model at the P and T
# Kraus et al. 2012 Table 2 State B stv on the melt curve at 4500 K
#stv_reference_point = [-1,70,4500,3567] # (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg) 
#stv_reference_point = [-1,70,4500,3267] # (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg) 
stv_reference_point = [-1,70,4500,3267] # (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg) 
#
# seifertite reference point
# tie it to the triple point for simpler calculations
# adjust melt curve slope and latent heat across CaCl2-sft
sft_reference_point = [-1,Pt3,Tt3,-1]  #  (rho,p,t,s)
#sft_dS = -400 # J/K/kg dS chosen to align to the stv solidus
sft_dS = -625 # J/K/kg dS chosen to align to the sft melt curve Park et al. 2025

#
# Solid models constant specific heat capacity
# cv = solid_cvfactor * (n=3) * Rgas
# Kraus et al. 2012 used 3.65
# In this work, the prefactor for cv and the reference solid entropy 
# are linked to be consistent with the slopes of the melt curve
stv_cvfactor = 3.65 # value cited by Kraus et al. 2012
sft_cvfactor = 3.5 # seifertite has a lower heat capacity than stishovite (Driver et al. 2010)
#
#########################################################################
# ***                                       ***
# ***       HUGONIOT Us-Up and Us-T FITS    ***
# ***                                       ***
# For each phase XXX, fits are calculated in hugoniot-fits-XXX.ipynb 
# Hugoniot fits are stored in individual data files in outputfitpath
# When satisfied, these fit data files are moved to modelfitpath 
# for use across all notebooks
#
# Checks included in each notebook to confirm using the latest fit
# ***
##################################################################
##################################################################
#  FUSED SILICA FITS
##################################################################
##################################################################
# NEW FITS FOR AMODEO ET AL.
# -------  Fused Silica Liquid Region -----------
# read in Us-up fit record
fs_fit_liq0 = hp.HPeos_Fit()
fs_fit_liq0.outputfitpath = modelfitpath # in load_paths.py
fs_fit_liq0.filename = 'fs_fit_liq_quadratic'
fs_fit_liq0.loadjson()

# create reference Hugoniot from fit
fs_hug_liq = hp.HPeos_Shock()
fs_hug_liq.rho_ref = 2.204 # g/cm3
fs_hug_liq.up = np.arange(3,20,.001) # km/s
fs_hug_liq.us = fs_fit_liq0.calcf(fs_hug_liq.up)
fs_hug_liq.calcPrho()

# read in T-Us fit record
fs_fitTus_liq0 = hp.HPeos_Fit()
fs_fitTus_liq0.outputfitpath = modelfitpath # in load_paths.py
fs_fitTus_liq0.filename = 'fs_fitTus_liq_cubic'
fs_fitTus_liq0.loadjson()
fs_hug_liq.T = fs_fitTus_liq0.calcf(fs_hug_liq.us) # calculate temperatures on the Hugoniot


# Kraus et al. 2012 entropy on fused silica Hugoniot in liquid region - only 70 to 110 
fs_hug_liq_K12 = deepcopy(fs_hug_liq)
fs_hug_liq_K12.S = 2284+26.82*fs_hug_liq.P-0.0723*fs_hug_liq.P**2
fs_hug_liq_K12.S_std = 127.3+0.328*fs_hug_liq.P-0.000183*fs_hug_liq.P*fs_hug_liq.P


if 0:
    # -------  Fused Silica Transition Region between stv and liquid -----------
    # read in fit record
    fs_fit_tra0 = hp.HPeos_Fit()
    fs_fit_tra0.outputfitpath = modelfitpath # in load_paths.py
    fs_fit_tra0.filename = 'fs_fit_trans_linear'
    fs_fit_tra0.loadjson()
    
    # Fused Silica TRANSITION stv-liquid region
    fs_hug_tra = hp.HPeos_Shock()
    fs_hug_tra.rho_ref = 2.204 # g/cm3
    fs_hug_tra.up = np.arange(3.8,4.9,.01) # km/s
    fs_hug_tra.us = fs_fit_tra0.calcf(fs_hug_tra.up)
    fs_hug_tra.calcPrho()
    # calculate temperatures in the transition region using the liquid Hugoniot fit
    fs_hug_tra.T = fs_fitTus_liq0.calcf(fs_hug_tra.us) 

# -------  Fused Silica Stishovite Region  -----------
# Fused Silica STV region
# read in fit record
fs_fit_stv0 = hp.HPeos_Fit()
fs_fit_stv0.outputfitpath = modelfitpath # in load_paths.py
fs_fit_stv0.filename = 'fs_fit_stv_linear'
fs_fit_stv0.loadjson()

fs_hug_stv = hp.HPeos_Shock()
fs_hug_stv.rho_ref = 2.204 # g/cm3
fs_hug_stv.up = np.arange(2.4,4.5,.01) # km/s
fs_hug_stv.us = fs_fit_stv0.calcf(fs_hug_stv.up)
fs_hug_stv.calcPrho()

fs_fitTus_stv0 = hp.HPeos_Fit()
fs_fitTus_stv0.outputfitpath = modelfitpath # in load_paths.py
fs_fitTus_stv0.filename = 'fs_fitTus_stv_quadratic'
fs_fitTus_stv0.loadjson()
fs_hug_stv.T = fs_fitTus_stv0.calcf(fs_hug_stv.us) # calculate temperatures on the solid branch Hugoniot

# -------  Fused Silica Transitions  -----------
#
# FS Pressure bounds on the Omega data wiggles
fs_wiggle_us1 = 7.5 # km/s # march 2025
fs_wiggle_us2 = 8. # km/s # march 2025
fs_wiggle_up1 = np.interp(fs_wiggle_us1,fs_hug_stv.us,fs_hug_stv.up)
fs_wiggle_up2 = np.interp(fs_wiggle_us2,fs_hug_liq.us,fs_hug_liq.up)
fs_wiggle_P1 = np.interp(fs_wiggle_us1,fs_hug_stv.us,fs_hug_stv.P)
fs_wiggle_P2 = np.interp(fs_wiggle_us2,fs_hug_liq.us,fs_hug_liq.P)
fs_wiggle_T1 = np.interp(fs_wiggle_us1,fs_hug_stv.us,fs_hug_stv.T)
fs_wiggle_T2 = np.interp(fs_wiggle_us2,fs_hug_liq.us,fs_hug_liq.T)

fs_wiggle_rho1 = np.interp(fs_wiggle_us1,fs_hug_stv.us,fs_hug_stv.rho)
fs_wiggle_rho2 = np.interp(fs_wiggle_us2,fs_hug_liq.us,fs_hug_liq.rho)

# FS enters liquid field
# in melt-curve-fits.ipynb
fs_startliq_P =  90.5101601196313  # GPa
fs_startliq_us = np.interp(fs_startliq_P,fs_hug_liq.P,fs_hug_liq.us)
fs_startliq_up = np.interp(fs_startliq_P,fs_hug_liq.P,fs_hug_liq.up)
fs_startliq_T = np.interp(fs_startliq_P,fs_hug_liq.P,fs_hug_liq.T)
fs_startliq_rho = np.interp(fs_startliq_P,fs_hug_liq.P,fs_hug_liq.rho)
# amorphous to stv transition
fs_trans_stv_P = 34 # GPa - onset of stv on FS Hugoniot Tracy et al. 2018 same as onset of single shock wave Wackerle 1962
fs_trans_stv_rho = np.interp(fs_trans_stv_P,fs_hug_stv.P,fs_hug_stv.rho)
fs_trans_stv_us = np.interp(fs_trans_stv_P,fs_hug_stv.P,fs_hug_stv.us)
fs_trans_stv_up = np.interp(fs_trans_stv_P,fs_hug_stv.P,fs_hug_stv.up)
fs_trans_stv_T = np.interp(fs_trans_stv_P,fs_hug_stv.P,fs_hug_stv.T)
# intersection with solidus rho-P and stv region Hugoniot
# in melt-curve-fits.ipynb
fs_startsol_P =  72.32862722602232  # GPa
fs_startsol_T = np.interp(fs_startsol_P,fs_hug_stv.P,fs_hug_stv.T)
fs_startsol_us = np.interp(fs_startsol_P,fs_hug_stv.P,fs_hug_stv.us)
fs_startsol_up = np.interp(fs_startsol_P,fs_hug_stv.P,fs_hug_stv.up)

fs_solidus_P = fs_startsol_P
fs_liquidus_P = fs_startliq_P 

##################################################################
##################################################################
#  QUARTZ FITS
##################################################################
##################################################################
# NEW FITS FOR AMODEO ET AL.
# -------  Quartz Liquid Region including supercooled shock front -----------
# read in fit record
qtz_fit_liq0 = hp.HPeos_Fit()
qtz_fit_liq0.outputfitpath = modelfitpath # in load_paths.py
qtz_fit_liq0.filename = 'qtz_fit_liq_cubic'
qtz_fit_liq0.loadjson()

qtz_hug_liq = hp.HPeos_Shock()
qtz_hug_liq.rho_ref = 2.65 # g/cm3
qtz_hug_liq.up = np.arange(3,35,.005) # km/s
qtz_hug_liq.us = qtz_fit_liq0.calcf(qtz_hug_liq.up)
qtz_hug_liq.calcPrho()

qtz_fitTus_liq0 = hp.HPeos_Fit()
qtz_fitTus_liq0.outputfitpath = modelfitpath # in load_paths.py
qtz_fitTus_liq0.filename = 'qtz_fitTus_liq_cubic'
qtz_fitTus_liq0.loadjson()
qtz_hug_liq.T = qtz_fitTus_liq0.calcf(qtz_hug_liq.us) # calculate temperatures on the solid branch Hugoniot

# Kraus entropy recalculated for Davies et al. 2020 Table 3:
qtz_hug_liq_K12 = deepcopy(qtz_hug_liq)
qtz_hug_liq_K12.S    =(2820.*np.power(qtz_hug_liq.P,-0.5)+468.3*np.power(qtz_hug_liq.P,0.5)-6.68*qtz_hug_liq.P-593.) # J/kg/K
qtz_hug_liq_K12.S_std    = (127.3+0.328*qtz_hug_liq.P-0.000193*np.power(qtz_hug_liq.P,2)) # J/K/kg


# -------  Quartz SOLID Region (dominated by CaCl2) -----------
# read in fit record
qtz_fit_solid0 = hp.HPeos_Fit()
qtz_fit_solid0.outputfitpath = modelfitpath # in load_paths.py
qtz_fit_solid0.filename = 'qtz_fit_solid_linear'
qtz_fit_solid0.loadjson()

qtz_hug_solid = hp.HPeos_Shock()
qtz_hug_solid.rho_ref = 2.65 # g/cm3
qtz_hug_solid.up = np.arange(2.5,5,.001) # km/s
qtz_hug_solid.us = qtz_fit_solid0.calcf(qtz_hug_solid.up)
qtz_hug_solid.calcPrho()

qtz_fitTus_solid0 = hp.HPeos_Fit()
qtz_fitTus_solid0.outputfitpath = modelfitpath # in load_paths.py
qtz_fitTus_solid0.filename = 'qtz_fitTus_solid_quadratic'
qtz_fitTus_solid0.loadjson()
qtz_hug_solid.T = qtz_fitTus_solid0.calcf(qtz_hug_solid.us) # calculate temperatures on the solid branch Hugoniot

# -------  Quartz CaCl2 Region -----------
# alpha-quartz CaCl2 region
# read in fit record
qtz_fit_cacl20 = hp.HPeos_Fit()
qtz_fit_cacl20.outputfitpath = modelfitpath # in load_paths.py
qtz_fit_cacl20.filename = 'qtz_fit_cacl2_linear'
qtz_fit_cacl20.loadjson()

qtz_hug_cacl2 = hp.HPeos_Shock()
qtz_hug_cacl2.rho_ref = 2.65 # g/cm3
qtz_hug_cacl2.up = np.arange(2.5,5,.001) # km/s
qtz_hug_cacl2.us = qtz_fit_cacl20.calcf(qtz_hug_cacl2.up)
qtz_hug_cacl2.calcPrho()

qtz_fitTus_cacl20 = hp.HPeos_Fit()
qtz_fitTus_cacl20.outputfitpath = modelfitpath # in load_paths.py
qtz_fitTus_cacl20.filename = 'qtz_fitTus_cacl2_linear'
qtz_fitTus_cacl20.loadjson()
qtz_hug_cacl2.T = qtz_fitTus_cacl20.calcf(qtz_hug_cacl2.us) # calculate temperatures on the solid branch Hugoniot


# -------  Quartz Stv Region -----------

# alpha-Quartz STV region
# read in fit record
qtz_fit_stv0 = hp.HPeos_Fit()
qtz_fit_stv0.outputfitpath = modelfitpath # in load_paths.py
qtz_fit_stv0.filename = 'qtz_fit_stv_linear'
qtz_fit_stv0.loadjson()

qtz_hug_stv = hp.HPeos_Shock()
qtz_hug_stv.rho_ref = 2.65 # g/cm3
qtz_hug_stv.up = np.arange(2.5,5,.001) # km/s
qtz_hug_stv.us = qtz_fit_stv0.calcf(qtz_hug_stv.up)
qtz_hug_stv.calcPrho()

qtz_fitTus_stv0 = hp.HPeos_Fit()
qtz_fitTus_stv0.outputfitpath = modelfitpath # in load_paths.py
qtz_fitTus_stv0.filename = 'qtz_fitTus_stv_quadratic'
qtz_fitTus_stv0.loadjson()
qtz_hug_stv.T = qtz_fitTus_stv0.calcf(qtz_hug_stv.us) # calculate temperatures on the solid branch Hugoniot


# -------  Quartz Transitions -----------


# Qtz Pressure bounds on the Omega data wiggles
qtz_wiggle_us1 = 9.29 # km/s
#qtz_wiggle_us1 = 9.22 # km/s
qtz_wiggle_us2 = 9.5 # km/s
qtz_wiggle_P1 = np.interp(qtz_wiggle_us1,qtz_hug_solid.us,qtz_hug_solid.P)
qtz_wiggle_up1 = np.interp(qtz_wiggle_us1,qtz_hug_solid.us,qtz_hug_solid.up)
qtz_wiggle_rho1 = np.interp(qtz_wiggle_us1,qtz_hug_solid.us,qtz_hug_solid.rho)
qtz_wiggle_T1 = np.interp(qtz_wiggle_us1,qtz_hug_solid.us,qtz_hug_solid.T)

qtz_wiggle_P2 = np.interp(qtz_wiggle_us2,qtz_hug_liq.us,qtz_hug_liq.P)
qtz_wiggle_up2 = np.interp(qtz_wiggle_us2,qtz_hug_liq.us,qtz_hug_liq.up)
qtz_wiggle_rho2 = np.interp(qtz_wiggle_us2,qtz_hug_liq.us,qtz_hug_liq.rho)
qtz_wiggle_T2 = np.interp(qtz_wiggle_us2,qtz_hug_liq.us,qtz_hug_liq.T)

# Hugoniot enters the liquid field
# in melt-curve-fits.ipynb
qtz_startliq_P =  145.80682496155464  # GPa
qtz_startliq_rho = np.interp(qtz_startliq_P,qtz_hug_liq.P,qtz_hug_liq.rho)
qtz_startliq_us = np.interp(qtz_startliq_P,qtz_hug_liq.P,qtz_hug_liq.us)
qtz_startliq_up = np.interp(qtz_startliq_P,qtz_hug_liq.P,qtz_hug_liq.up)
qtz_startliq_T = np.interp(qtz_startliq_P,qtz_hug_liq.P,qtz_hug_liq.T)

# Estimate transition to CaCl2 region (XRD indicates metastable phases form)
# in melt-curve-fits.ipynb
qtz_trans_cacl2_P =  82.58142670722148  # GPa
qtz_trans_cacl2_us = np.interp(qtz_trans_cacl2_P,qtz_hug_stv.P,qtz_hug_stv.us)
qtz_trans_cacl2_T = np.interp(qtz_trans_cacl2_P,qtz_hug_stv.P,qtz_hug_stv.T)
qtz_trans_cacl2_rho = np.interp(qtz_trans_cacl2_P,qtz_hug_stv.P,qtz_hug_stv.rho)
qtz_trans_cacl2_up = np.interp(qtz_trans_cacl2_P,qtz_hug_stv.P,qtz_hug_stv.up)

qtz_low_cacl2_T = np.interp(qtz_trans_cacl2_P,qtz_hug_cacl2.P,qtz_hug_cacl2.T)

# in melt-curve-fits.ipynb
qtz_startsol_P =  100.50773755669458  # GPa
qtz_startsol_us = np.interp(qtz_startsol_P,qtz_hug_solid.P,qtz_hug_solid.us)
qtz_startsol_up = np.interp(qtz_startsol_P,qtz_hug_solid.P,qtz_hug_solid.up)
qtz_startsol_rho = np.interp(qtz_startsol_P,qtz_hug_solid.P,qtz_hug_solid.rho)
qtz_startsol_T = np.interp(qtz_startsol_P,qtz_hug_solid.P,qtz_hug_solid.T)

qtz_solidus_P = qtz_startsol_P
qtz_liquidus_P = qtz_startliq_P 

qtz_singleshock_P = 70 # GPa Wackerle 1962 70.3 GPa lowest single shock datum
qtz_singleshock_up = np.interp(qtz_singleshock_P,qtz_hug_solid.P,qtz_hug_solid.up)
qtz_singleshock_us = np.interp(qtz_singleshock_P,qtz_hug_solid.P,qtz_hug_solid.us)
qtz_singleshock_T = np.interp(qtz_singleshock_P,qtz_hug_solid.P,qtz_hug_solid.T)

# Estimated triple point for L-Stv-CaCl2
# from GGA phase boundaries Tsuchiya et al.
tpcacl2_P = 89.32772328

#######################################################################################
######  COESITE HUGONIOT FITS
#######################################################################################

# Thermal inversion
cst_wiggle_us1 = 9.8 # km/s
cst_wiggle_us2 = 10.3 # km/s

# -------  Coesite Hugoniot CaCl2 Region -----------

# read in fit record
cst_fit_cacl20 = hp.HPeos_Fit()
cst_fit_cacl20.outputfitpath = modelfitpath # in load_paths.py
cst_fit_cacl20.filename = 'cst_fit_cacl2_linear'
#cst_fit_cacl20.filename = 'cst_fit_stv_linear'
cst_fit_cacl20.loadjson()

cst_hug_cacl2 = hp.HPeos_Shock()
cst_hug_cacl2.rho_ref = 2.92 # g/cm3
cst_hug_cacl2.up = np.arange(2.5,5,.001) # km/s
cst_hug_cacl2.us = cst_fit_cacl20.calcf(cst_hug_cacl2.up)
cst_hug_cacl2.calcPrho()

cst_fitTus_cacl20 = hp.HPeos_Fit()
cst_fitTus_cacl20.outputfitpath = modelfitpath # in load_paths.py
cst_fitTus_cacl20.filename = 'cst_fitTus_cacl2_quadratic'
#cst_fitTus_cacl20.filename = 'cst_fitTus_stv_quadratic'
cst_fitTus_cacl20.loadjson()
cst_hug_cacl2.T = cst_fitTus_cacl20.calcf(cst_hug_cacl2.us) # calculate temperatures on the solid branch Hugoniot

# -------  Coesite Hugoniot Liquid Region -----------

# read in fit record
cst_fit_liq0 = hp.HPeos_Fit()
cst_fit_liq0.outputfitpath = modelfitpath # in load_paths.py
cst_fit_liq0.filename = 'cst_fit_liq_linear'
cst_fit_liq0.loadjson()

cst_hug_liq = hp.HPeos_Shock()
cst_hug_liq.rho_ref = 2.92 # g/cm3
cst_hug_liq.up = np.arange(4,14,.001) # km/s
cst_hug_liq.us = cst_fit_liq0.calcf(cst_hug_liq.up)
cst_hug_liq.calcPrho()

cst_fitTus_liq0 = hp.HPeos_Fit()
cst_fitTus_liq0.outputfitpath = modelfitpath # in load_paths.py
cst_fitTus_liq0.filename = 'cst_fitTus_liq_quadratic'
cst_fitTus_liq0.loadjson()
cst_hug_liq.T = cst_fitTus_liq0.calcf(cst_hug_liq.us) # calculate temperatures on the solid branch Hugoniot


# -------  Coesite intersection points -----------

cst_wiggle_up1 = np.interp(cst_wiggle_us1,cst_hug_cacl2.us,cst_hug_cacl2.up)
cst_wiggle_rho1 = np.interp(cst_wiggle_us1,cst_hug_cacl2.us,cst_hug_cacl2.rho)
cst_wiggle_P1 = np.interp(cst_wiggle_us1,cst_hug_cacl2.us,cst_hug_cacl2.P)

cst_wiggle_up2 = np.interp(cst_wiggle_us2,cst_hug_liq.us,cst_hug_liq.up)
cst_wiggle_P2 = np.interp(cst_wiggle_us2,cst_hug_liq.us,cst_hug_liq.P)
cst_wiggle_rho2 = np.interp(cst_wiggle_us2,cst_hug_liq.us,cst_hug_liq.rho)

# in melt-curve-fits.ipynb
cst_solidus_P =  125.53423196515627  # GPa
cst_solidus_T = np.interp(cst_solidus_P,cst_hug_cacl2.P,cst_hug_cacl2.T)


# from melt-curve-fits.ipynb
cst_liquidus_P =  183.61347134392335  # GPa
cst_liquidus_T = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.T)
cst_liquidus_rho = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.rho)
cst_liquidus_up = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.up)
cst_liquidus_us = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.us)

# intersection point between solid Hugoniot and solidus
cst_sol_cross = hp.HPeos_Shock()
cst_sol_cross.P = cst_solidus_P
cst_sol_cross.up = np.interp(cst_solidus_P,cst_hug_cacl2.P,cst_hug_cacl2.up)
cst_sol_cross.us = np.interp(cst_solidus_P,cst_hug_cacl2.P,cst_hug_cacl2.us)
cst_sol_cross.T = np.interp(cst_solidus_P,cst_hug_cacl2.P,cst_hug_cacl2.T)
cst_sol_cross.rho = np.interp(cst_solidus_P,cst_hug_cacl2.P,cst_hug_cacl2.rho)

# intersection point between liquid Hugoniot and liquidus
cst_liq_cross = hp.HPeos_Shock()
cst_liq_cross.P = cst_liquidus_P
cst_liq_cross.T = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.T)
cst_liq_cross.up = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.up)
cst_liq_cross.us = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.us)
cst_liq_cross.rho = np.interp(cst_liquidus_P,cst_hug_liq.P,cst_hug_liq.rho)


# -------  Stishovite Hugoniot Liquid Region -----------

# read in fit record
stv_fit_liq0 = hp.HPeos_Fit()
stv_fit_liq0.outputfitpath = modelfitpath # in load_paths.py
stv_fit_liq0.filename = 'stv_fit_liq_linear_up18'
stv_fit_liq0.loadjson()

stv_hug_liq = hp.HPeos_Shock()
stv_hug_liq.rho_ref = stv_fit_liq0.rho_ref # g/cm3
stv_hug_liq.up = np.arange(5,20,.01) # km/s
stv_hug_liq.us = stv_fit_liq0.calcf(stv_hug_liq.up)
stv_hug_liq.calcPrho()

stv_fitTus_liq0 = hp.HPeos_Fit()
stv_fitTus_liq0.outputfitpath = modelfitpath # in load_paths.py
stv_fitTus_liq0.filename = 'stv_fitTus_liq_fourthpoly'
stv_fitTus_liq0.loadjson()
stv_hug_liq.T = stv_fitTus_liq0.calcf(stv_hug_liq.us) # calculate temperatures on the solid branch Hugoniot
stv_hug_liq.T_std = np.interp(stv_hug_liq.us,shockdat[iM15stv.i].us,shockdat[iM15stv.i].T_std)

# -------  Stishovite intersection points -----------

# these are just based on the wiggles in the decaying shock from Millot et al. 2015
stv_wiggle2 = hp.HPeos_Point()
stv_wiggle2.us = 17.4 # km/s
stv_wiggle2.P = np.interp(stv_wiggle2.us,stv_hug_liq.us,stv_hug_liq.P)
stv_wiggle2.rho = np.interp(stv_wiggle2.us,stv_hug_liq.us,stv_hug_liq.rho)
stv_wiggle2.up = np.interp(stv_wiggle2.us,stv_hug_liq.us,stv_hug_liq.up)
stv_wiggle2.note = 'km/s, km/s, GPa, g/cm3'

# new estimate of solidus, liquidus crossings May 2025
# linear_up18 fit and shifted Gonzales et al. 2016 melt curve with conservative errors
stv_liquidus_P =  632.9429427843702  # GPa this is fitted to P-T which is similar for the cottunite and fitted curves
# assumes 2% increase in density for solid
# stv_solidus_P =  609.2078653823245  # GPa this is for the fitted solidus density that is lower density than the DFT points in this pressure range
stv_solidus_P =  490.63813073715124  # GPa this is fitted directly to the cottunite density points shifted by 2%

#######################################################################################
######  END OF FILE
#######################################################################################



