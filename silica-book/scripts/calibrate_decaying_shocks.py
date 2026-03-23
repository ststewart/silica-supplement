# Calibrate decaying shock data
#
from copy import deepcopy
import numpy as np
#####################
# Knudson and Desjarlais PRB 2013

qtz_hug_KD      = hp.HPeos_Shock()
qtz_hug_KD.rho_ref = 2.650 # g/cm3
qtz_hug_KD.up   = np.arange(2.5,10.,.001) # km/s
# cubic fit from Knudson and Desjarlais PRB 2013
res=[1.754,1.862,-3.364e-2,5.666e-4]
qtz_hug_KD.us = res[0]+res[1]*qtz_hug_KD.up+res[2]*(qtz_hug_KD.up**2)+res[3]*(qtz_hug_KD.up**3)
qtz_hug_KD.P    = qtz_hug_KD.rho_ref*qtz_hug_KD.up*qtz_hug_KD.us # Pa
qtz_hug_KD.rho  = 1/((1/qtz_hug_KD.rho_ref)*(1-(qtz_hug_KD.up/qtz_hug_KD.us))) # kg/m3


# Calibrate Millot & Davies unpublished decaying shocks in alpha-quartz, z cut
# high pressure uses the liquid branch of the new Hugoniot fit
shockdat_list[MDqtz70811['id']].P = np.interp(shockdat_list[MDqtz70811['id']].us,qtz_hug_liq.us,qtz_hug_liq.P)
shockdat_list[MDqtz70811['id']].rho = np.interp(shockdat_list[MDqtz70811['id']].us,qtz_hug_liq.us,qtz_hug_liq.rho)
shockdat_list[MDqtz70811['id']].up = np.interp(shockdat_list[MDqtz70811['id']].us,qtz_hug_liq.us,qtz_hug_liq.up)

shockdat_list[MDqtz70813['id']].P = np.interp(shockdat_list[MDqtz70813['id']].us,qtz_hug_liq.us,qtz_hug_liq.P)
shockdat_list[MDqtz70813['id']].rho = np.interp(shockdat_list[MDqtz70813['id']].us,qtz_hug_liq.us,qtz_hug_liq.rho)
shockdat_list[MDqtz70813['id']].up = np.interp(shockdat_list[MDqtz70813['id']].us,qtz_hug_liq.us,qtz_hug_liq.up)

shockdat_list[MDqtz70817['id']].P = np.interp(shockdat_list[MDqtz70817['id']].us,qtz_hug_liq.us,qtz_hug_liq.P)
shockdat_list[MDqtz70817['id']].rho = np.interp(shockdat_list[MDqtz70817['id']].us,qtz_hug_liq.us,qtz_hug_liq.rho)
shockdat_list[MDqtz70817['id']].up = np.interp(shockdat_list[MDqtz70817['id']].us,qtz_hug_liq.us,qtz_hug_liq.up)

shockdat_list[MDqtz70820['id']].P = np.interp(shockdat_list[MDqtz70820['id']].us,qtz_hug_liq.us,qtz_hug_liq.P)
shockdat_list[MDqtz70820['id']].rho = np.interp(shockdat_list[MDqtz70820['id']].us,qtz_hug_liq.us,qtz_hug_liq.rho)
shockdat_list[MDqtz70820['id']].up = np.interp(shockdat_list[MDqtz70820['id']].us,qtz_hug_liq.us,qtz_hug_liq.up)

#M70811liq = deepcopy(shockdat_list[MDqtz70811['id']])
#M70813liq = deepcopy(shockdat_list[MDqtz70813['id']])
#M70817liq = deepcopy(shockdat_list[MDqtz70817['id']])
#M70820liq = deepcopy(shockdat_list[MDqtz70820['id']])

ind =np.where(shockdat_list[MDqtz70811['id']].P < qtz_solidus_P)[0]
if len(ind)>0:
    shockdat_list[MDqtz70811['id']].P[ind] = np.interp(shockdat_list[MDqtz70811['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.P)
    shockdat_list[MDqtz70811['id']].rho[ind] = np.interp(shockdat_list[MDqtz70811['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.rho)
    shockdat_list[MDqtz70811['id']].up[ind] = np.interp(shockdat_list[MDqtz70811['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.up)

ind =np.where(shockdat_list[MDqtz70813['id']].P < qtz_solidus_P)[0]
if len(ind)>0:
    shockdat_list[MDqtz70813['id']].P[ind] = np.interp(shockdat_list[MDqtz70813['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.P)
    shockdat_list[MDqtz70813['id']].rho[ind] = np.interp(shockdat_list[MDqtz70813['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.rho)
    shockdat_list[MDqtz70813['id']].up[ind] = np.interp(shockdat_list[MDqtz70813['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.up)

ind =np.where(shockdat_list[MDqtz70817['id']].P < qtz_solidus_P)[0]
if len(ind)>0:
    shockdat_list[MDqtz70817['id']].P[ind] = np.interp(shockdat_list[MDqtz70817['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.P)
    shockdat_list[MDqtz70817['id']].rho[ind] = np.interp(shockdat_list[MDqtz70817['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.rho)
    shockdat_list[MDqtz70817['id']].up[ind] = np.interp(shockdat_list[MDqtz70817['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.up)

ind =np.where(shockdat_list[MDqtz70820['id']].P < qtz_solidus_P)[0]
if len(ind)>0:
    shockdat_list[MDqtz70820['id']].P[ind] = np.interp(shockdat_list[MDqtz70820['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.P)
    shockdat_list[MDqtz70820['id']].rho[ind] = np.interp(shockdat_list[MDqtz70820['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.rho)
    shockdat_list[MDqtz70820['id']].up[ind] = np.interp(shockdat_list[MDqtz70820['id']].us[ind],qtz_hug_solid.us,qtz_hug_solid.up)

M70811 = deepcopy(shockdat_list[MDqtz70811['id']])
M70813 = deepcopy(shockdat_list[MDqtz70813['id']])
M70817 = deepcopy(shockdat_list[MDqtz70817['id']])
M70820 = deepcopy(shockdat_list[MDqtz70820['id']])


# Calibrate Hicks et al. 2006 Decaying Shocks
shockdat_list[H06fs['id']].P    = np.interp(shockdat_list[H06fs['id']].us,fs_hug_liq.us,fs_hug_liq.P)
shockdat_list[H06fs['id']].rho  = np.interp(shockdat_list[H06fs['id']].us,fs_hug_liq.us,fs_hug_liq.rho)
shockdat_list[H06fs['id']].up   = np.interp(shockdat_list[H06fs['id']].us,fs_hug_liq.us,fs_hug_liq.up)

shockdat_list[H06qtz['id']].P   = np.interp(shockdat_list[H06qtz['id']].us,qtz_hug_liq.us,qtz_hug_liq.P)
shockdat_list[H06qtz['id']].rho = np.interp(shockdat_list[H06qtz['id']].us,qtz_hug_liq.us,qtz_hug_liq.rho)
shockdat_list[H06qtz['id']].up  = np.interp(shockdat_list[H06qtz['id']].us,qtz_hug_liq.us,qtz_hug_liq.up)


# Calibrate Coesite Decaying shock temperature data from Feng et al. 2024

C43939.P = np.interp(C43939.us,cst_hug_liq.us,cst_hug_liq.P)
C43941.P = np.interp(C43941.us,cst_hug_liq.us,cst_hug_liq.P)
C43943.P = np.interp(C43943.us,cst_hug_liq.us,cst_hug_liq.P)

C43939.up = np.interp(C43939.us,cst_hug_liq.us,cst_hug_liq.up)
C43941.up = np.interp(C43941.us,cst_hug_liq.us,cst_hug_liq.up)
C43943.up = np.interp(C43943.us,cst_hug_liq.us,cst_hug_liq.up)

C43939.rho = np.interp(C43939.us,cst_hug_liq.us,cst_hug_liq.rho)
C43941.rho = np.interp(C43941.us,cst_hug_liq.us,cst_hug_liq.rho)
C43943.rho = np.interp(C43943.us,cst_hug_liq.us,cst_hug_liq.rho)

# only this shot had pressures below the solidus
ind = np.where(C43939.us<cst_sol_cross.us)[0]
if len(ind)>1:
    C43939.P[ind] = np.interp(C43939.us[ind],cst_hug_cacl2.us,cst_hug_cacl2.P)
    C43939.up[ind] = np.interp(C43939.us[ind],cst_hug_cacl2.us,cst_hug_cacl2.up)
    C43939.rho[ind] = np.interp(C43939.us[ind],cst_hug_cacl2.us,cst_hug_cacl2.rho)


# Coesite fit from Feng et al. 2024
# two part fit
cst_hug_F24      = hp.HPeos_Shock()
cst_hug_F24.rho_ref = 2.920 # kg/cm3
cst_hug_F24.up   = np.arange(2.69,15.,.01) # km/s
# Us = 2.52 + 1.70Up, 7.1 km/s < Us < 10.4 km/s. 
# Us = 2.097(±0.099) + 1.815(±0.022) Up – 0.020(±0.001) U2p, Us ≥ 10.4 km/s.
cst_hug_F24.us   = 2.097 + 1.815*cst_hug_F24.up - 0.02*(cst_hug_F24.up)*(cst_hug_F24.up) # km/s
ind = np.where(cst_hug_F24.us < 10.4) # km/s
cst_hug_F24.us[ind]   = 2.52 + 1.70*cst_hug_F24.up[ind]
cst_hug_F24.P    = cst_hug_F24.rho_ref*cst_hug_F24.up*cst_hug_F24.us # GPa
cst_hug_F24.rho  = 1/((1/cst_hug_F24.rho_ref)*(1-(cst_hug_F24.up/cst_hug_F24.us))) # g/m3

# Calibrate stishovite decaying shock temperature data from Millot et al. 2015

stv_hug_M15rev = deepcopy(shockdat[iM15stv.i])
stv_hug_M15rev.P = np.interp(stv_hug_M15rev.us,stv_hug_liq.us,stv_hug_liq.P)
stv_hug_M15rev.up = np.interp(stv_hug_M15rev.us,stv_hug_liq.us,stv_hug_liq.up)
stv_hug_M15rev.rho = np.interp(stv_hug_M15rev.us,stv_hug_liq.us,stv_hug_liq.rho)


######## END OF FILE ########
