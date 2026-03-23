# hugoniot-calcs.py



# Hicks, D.G., Boehly, T.R., Eggert, J.H., Miller, J.E., Celliers, P.M., Collins, G.W., 2006. Dissociation of liquid silica at high pressures and temperatures. Phys. Rev. Lett. 97, 3–6. https://doi.org/10.1103/PhysRevLett.97.025502
# Hicks Quartz Shock Temperatures PRL 2006
fs_hug_Hicks06_file = silicadatadir+'Hicks2006-qtz-fs-hug.txt'
qtz_hug_Hicks06 = hp.HPeos_Shock()
qtz_hug_Hicks06.rho_ref = 2.65 # g/cm3
qtz_hug_Hicks06.us = np.loadtxt(fs_hug_Hicks06_file,delimiter=',',skiprows=5,usecols=[0]) # km/s
qtz_hug_Hicks06.T = np.loadtxt(fs_hug_Hicks06_file,delimiter=',',skiprows=5,usecols=[5]) # K
qtz_hug_Hicks06.T_std = np.loadtxt(fs_hug_Hicks06_file,delimiter=',',skiprows=5,usecols=[6]) # K

# this is set in define_eos.py
#qtz_hug_Hicks06.P = np.interp(qtz_hug_Hicks06.us,qtz_hug_KD.us,qtz_hug_KD.P) # GPa
#qtz_hug_Hicks06.up = np.interp(qtz_hug_Hicks06.us,qtz_hug_KD.us,qtz_hug_KD.up) # km/s
#qtz_hug_Hicks06.rho = np.interp(qtz_hug_Hicks06.us,qtz_hug_KD.us,qtz_hug_KD.rho) # g/cm3

# FUSED SILICA FS
# Root, S., Townsend, J.P., Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. J. Appl. Phys. 126. https://doi.org/10.1063/1.5126205
# Root et al. JAP 2019; entropy Kraus et al. 2012
fs_hug_RTK19  = hp.HPeos_Shock()
fs_hug_RTK19.rho_ref = 2.204 # g/cm3
fs_hug_RTK19.up   = np.arange(4.,18.,.1) # km/s
a = 4.925 # km/s
b = 1.234
c = 2.074
d=0.352 # 1/(km/s)
fs_hug_RTK19.us   = (a + b*fs_hug_RTK19.up - c*(fs_hug_RTK19.up)*np.exp(-d*(fs_hug_RTK19.up))) # km/s
fs_hug_RTK19.P    = fs_hug_RTK19.rho_ref*fs_hug_RTK19.up*fs_hug_RTK19.us # GPa
fs_hug_RTK19.rho  = 1/((1/fs_hug_RTK19.rho_ref)*(1-(fs_hug_RTK19.up/fs_hug_RTK19.us))) # kg/m3
# using Kaitlyn fit to the temperatures boslough + hicks -- only a small range in pressure; melting to XXX GPa
fs_hug_RTK19.T = -1.99167591e+04 + (6.37525746e+03*(fs_hug_RTK19.us)) + (-6.71164718e+02*(fs_hug_RTK19.us)**2) + ((3.64336344e+01)*(fs_hug_RTK19.us)**3) + ((-5.41405759e-01)*(fs_hug_RTK19.us)**4)
# Kraus et al. 2012 entropy on fused silica Hugoniot in liquid region - only 70 to 110 GPa
fs_hug_RTK19.S = 2284+26.82*fs_hug_RTK19.P-0.0723*fs_hug_RTK19.P**2
fs_hug_RTK19.S_std = 127.3+0.328*fs_hug_RTK19.P-0.000183*fs_hug_RTK19.P*fs_hug_RTK19.P


# Hicks, D.G., Boehly, T.R., Eggert, J.H., Miller, J.E., Celliers, P.M., Collins, G.W., 2006. Dissociation of liquid silica at high pressures and temperatures. Phys. Rev. Lett. 97, 3–6. https://doi.org/10.1103/PhysRevLett.97.025502
# Hicks Fused Quartz Shock Temperatures PRL 2006
fs_hug_Hicks06_file = silicadatadir+'Hicks2006-qtz-fs-hug.txt'
fs_hug_Hicks06 = hp.HPeos_Shock()
fs_hug_Hicks06.rho_ref = 2.204 # g/cm3
fs_hug_Hicks06.us = np.loadtxt(fs_hug_Hicks06_file,delimiter=',',skiprows=5,usecols=[0]) # km/s
fs_hug_Hicks06.T = np.loadtxt(fs_hug_Hicks06_file,delimiter=',',skiprows=5,usecols=[7]) # K
fs_hug_Hicks06.T_std = np.loadtxt(fs_hug_Hicks06_file,delimiter=',',skiprows=5,usecols=[8]) # K
# this is set in define_eos.py
#fs_hug_Hicks06.P = np.interp(fs_hug_Hicks06.us,fs_hug_RTK19.us,fs_hug_RTK19.P) # GPa
#fs_hug_Hicks06.up = np.interp(fs_hug_Hicks06.us,fs_hug_RTK19.us,fs_hug_RTK19.up) # km/s
#fs_hug_Hicks06.rho = np.interp(fs_hug_Hicks06.us,fs_hug_RTK19.us,fs_hug_RTK19.rho) # g/cm3


if 0:
    # update Millot decaying shocks with state info
    
    qtz_hug_KD      = hp.HPeos_Shock()
    qtz_hug_KD.rho_ref = 2.650 # g/cm3
    qtz_hug_KD.up   = np.arange(2.5,10.,.001) # km/s
    if (qtz_liq_hug == 'KD13') or  (qtz_liq_hug == 'NEW'):
        # cubic fit from Knudson and Desjarlais PRB 2013
        # NOTE FIT INCLUDES STRANGE LOW TEMP HIGH ENTROPY PHASE DATA
        res=[1.754,1.862,-3.364e-2,5.666e-4]
        qtz_hug_KD.us = res[0]+res[1]*qtz_hug_KD.up+res[2]*(qtz_hug_KD.up**2)+res[3]*(qtz_hug_KD.up**3)
    if qtz_liq_hug == 'KD09':
        # exponential fit from Knudson and Desjarlais PRL 2009
        a = 6.26 # km/s
        b = 1.2
        c = 2.56
        d=0.37 # 1/(km/s)
        qtz_hug_KD.us   = (a + b*qtz_hug_KD.up - c*(qtz_hug_KD.up)*np.exp(-d*(qtz_hug_KD.up))) # km/s
    qtz_hug_KD.P    = qtz_hug_KD.rho_ref*qtz_hug_KD.up*qtz_hug_KD.us # Pa
    qtz_hug_KD.rho  = 1/((1/qtz_hug_KD.rho_ref)*(1-(qtz_hug_KD.up/qtz_hug_KD.us))) # kg/m3
    
    if 0:
        # Knudson & Desjarlais PRL 2009; temps Brygoo et al. JAP 2015; entropy Kraus et al. 2012
        qtz_hug_KD      = hp.HPeos_Shock()
        qtz_hug_KD.rho_ref = 2.650 # g/cm3
        #qtz_hug_KD.up   = np.arange(2.5,22.,.25) # km/s
        qtz_hug_KD.up   = np.arange(2.5,10.,.001) # km/s
        a = 6.26 # km/s
        b = 1.2
        c = 2.56
        d=0.37 # 1/(km/s)
        qtz_hug_KD.us   = (a + b*qtz_hug_KD.up - c*(qtz_hug_KD.up)*np.exp(-d*(qtz_hug_KD.up))) # km/s
        qtz_hug_KD.P    = qtz_hug_KD.rho_ref*qtz_hug_KD.up*qtz_hug_KD.us # Pa
        qtz_hug_KD.rho  = 1/((1/qtz_hug_KD.rho_ref)*(1-(qtz_hug_KD.up/qtz_hug_KD.us))) # kg/m3
    qtz_hug_KD.T    = 1860+3.56*np.power(qtz_hug_KD.us,3.036) # K Brygoo et al. JAP 2015
    # Kraus entropy only good from 113 to 800 GPa -- goes negative at higher values Kraus et al. JGR 2012
    #qtz_hug_KD.S    = (1793.+22.40*qtz_hug_KD.P-0.0426*np.power(qtz_hug_KD.P,2)+4.319E-5*np.power(qtz_hug_KD.P,3)-1.742E-8*np.power(qtz_hug_KD.P,4))/1.e3 # kJ/kg/K
    tmp_s    = (1793.+22.40*qtz_hug_KD.P-0.0426*np.power(qtz_hug_KD.P,2)+4.319E-5*np.power(qtz_hug_KD.P,3)-1.742E-8*np.power(qtz_hug_KD.P,4))/1.e3 # kJ/kg/K
    # Kraus entropy recalculated for Davies et al. 2020 Table 3:
    qtz_hug_KD.S    = (2820.*np.power(qtz_hug_KD.P,-0.5)+468.3*np.power(qtz_hug_KD.P,0.5)-6.68*qtz_hug_KD.P-593.)/1.e3 # kJ/kg/K
    qtz_hug_KD.S_std    = (127.3+0.328*qtz_hug_KD.P-0.000193*np.power(qtz_hug_KD.P,2))/1.e3 # kJ/K/kg
    

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


## Calculate specific entropies on the quartz Hugoniot

if 1:
    # calculate entropies
    # refpt is 4500 K; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(qtz_hug_solid.rho)
    garr = np.interp(rarr,stvpvt.rhoarr,stvpvt.garr)
    tarr = deepcopy(qtz_hug_solid.T)
    parr = deepcopy(qtz_hug_solid.P)
    sarr = np.zeros(len(qtz_hug_solid.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        rind = np.where(stvpvt.rhoarr >= rho)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = stvpvt.st4500[rind]
        sisen = stvpvt.refpt[3]
        dS = stvpvt.cv*(np.log(tarr[iii])-np.log(tisen))
        sarr[iii]=sisen+dS
    qtz_hug_solid.S = deepcopy(sarr)


if 1:
    # calculate entropies
    # refpt is 4500 K; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(qtz_hug_cacl2.rho)
    garr = np.interp(rarr,stvpvt.rhoarr,stvpvt.garr)
    tarr = deepcopy(qtz_hug_cacl2.T)
    parr = deepcopy(qtz_hug_cacl2.P)
    sarr = np.zeros(len(qtz_hug_cacl2.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        rind = np.where(stvpvt.rhoarr >= rho)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = stvpvt.st4500[rind]
        sisen = stvpvt.refpt[3]
        dS = stvpvt.cv*(np.log(tarr[iii])-np.log(tisen))
        sarr[iii]=sisen+dS
    qtz_hug_cacl2.S = deepcopy(sarr)


if 1:
    # calculate entropies LOW BRANCH OF QTZ
    # refpt is 4500 K; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(qtz_hug_stv.rho)
    garr = np.interp(rarr,stvpvt.rhoarr,stvpvt.garr)
    tarr = deepcopy(qtz_hug_stv.T)
    parr = deepcopy(qtz_hug_stv.P)
    sarr = np.zeros(len(qtz_hug_stv.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        rind = np.where(stvpvt.rhoarr >= rho)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = stvpvt.st4500[rind]
        sisen = stvpvt.refpt[3]
        dS = stvpvt.cv*(np.log(tarr[iii])-np.log(tisen))
        sarr[iii]=sisen+dS
    qtz_hug_stv.S = deepcopy(sarr)


if 0:
    # calculate entropies SFT BRANCH OF QTZ
    # refpt is 5900 K; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(qtz_hug_sft.rho)
    garr = np.interp(rarr,sftpvt.rhoarr,sftpvt.garr)
    tarr = deepcopy(qtz_hug_sft.T)
    parr = deepcopy(qtz_hug_sft.P)
    sarr = np.zeros(len(qtz_hug_sft.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        rind = np.where(sftpvt.rhoarr >= rho)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = sftpvt.st5900[rind]
        sisen = sftpvt.refpt[3]
        dS = sftpvt.cv*(np.log(tarr[iii])-np.log(tisen))
        sarr[iii]=sisen+dS
    qtz_hug_sft.S = deepcopy(sarr)

## Calculate specific entropies on the fused silica Hugoniot

if 1:
    # refpt is 4500k; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(fs_hug_stv.rho)
    garr = np.interp(rarr,stvpvt.rhoarr,stvpvt.garr)
    tarr = deepcopy(fs_hug_stv.T)
    parr = deepcopy(fs_hug_stv.P)
    sarr = np.zeros(len(fs_hug_stv.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        rind = np.where(stvpvt.rhoarr >= rho)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = stvpvt.st4500[rind]
        sisen = stvpvt.refpt[3]
        dS = stvpvt.cv*(np.log(tarr[iii])-np.log(tisen))
        #print(iii,rho,stvpvt.rhoarr[rind],tisen,tarr[iii],sisen,dS)
        sarr[iii]=sisen+dS
    fs_hug_stv.S = deepcopy(sarr)

# =======================================================================================
# CALCULATE ENTROPY ON THE NEW FUSED SILICA LIQUID HUGONIOT FIT CHANGE TO REFPT
# dS = dE/T + PdV/T
# 
if 1:
    # NO: NEED TO USE KARKI MODEL FOR SELF-CONSISTENCY WITH MELT CURVE
    if 0:
        # set in define_eos.py fs_p_startliq = 80 # GPa
        P_start = fs_wiggle_P2#fs_p_startliq # GPa based on wiggle
        rho_start = np.interp(P_start,stsmcsmooth.P,stsmcsmooth.rho)
        S_start = np.interp(P_start,stsmcsmooth.P,stsmcsmooth.S_s)
        print(P_start,S_start,rho_start)
    if 1:
        # calculate the starting point relative to the karkipvt.refpt
        #P_start = fs_wiggle_P2#fs_p_startliq # GPa based on wiggle
        P_start = fs_wiggle_P2#fs_p_startliq # GPa based on wiggle
        rho_start = np.interp(P_start,fs_hug_liq.P,fs_hug_liq.rho)
        t_start = np.interp(P_start,fs_hug_liq.P,fs_hug_liq.T)
        # calculate S from Karki
        rind = np.where(karkipvt.rhoarr >= rho_start)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = karkipvt.st4500[rind]
        sisen = karkipvt.refpt[3]
        dS = karkipvt.cv*(np.log(t_start)-np.log(tisen))
        S_start = sisen+dS
    
        
    #ind = np.where(fs_hug_liq.rho >= rho_start)[0]
    ind = np.where(fs_hug_liq.P >= P_start)[0]
    fs_hug_liq.S = np.zeros(len(fs_hug_liq.rho))
    fs_hug_liq.S[ind[0]]=S_start # J/K/kg
    indref = ind[0]
    varr = 1/(1000.*fs_hug_liq.rho) # m3/kg
    parr = fs_hug_liq.P*1.e9 # Pa
    tarr = fs_hug_liq.T # K
    uarr = fs_hug_liq.U # J/kg
    for iii in range(1,len(ind)):
        fs_hug_liq.S[ind[iii]] = fs_hug_liq.S[ind[iii-1]] + (uarr[ind[iii]]-uarr[ind[iii-1]])/tarr[ind[iii]] + \
            parr[ind[iii]]*(varr[ind[iii]]-varr[ind[iii-1]])/tarr[ind[iii]]
    #ind = np.concatenate((np.where(fs_hug_liq.rho < rho_start)[0],[indref]))
    ind = np.concatenate((np.where(fs_hug_liq.P < P_start)[0],[indref]))
    for iii in range(len(ind)-2,-1,-1):
        fs_hug_liq.S[ind[iii]] = fs_hug_liq.S[ind[iii+1]] + (uarr[ind[iii]]-uarr[ind[iii+1]])/tarr[ind[iii]] + \
            parr[ind[iii]]*(varr[ind[iii]]-varr[ind[iii+1]])/tarr[ind[iii]]

if 1:
    fs_hug_liq_karki = deepcopy(fs_hug_liq)
    # FUSED SILICA LIQUID BRANCH REFERENCE TO KARKI LIQUID MODEL - SHOULD BE SIMILAR TO INTEGRATING THE HUGONIOT
    # calculate entropies USE LIQUID REFPT everywhere
    # refpt is 4500k; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(fs_hug_liq.rho)
    garr = np.interp(rarr,karkipvt.rhoarr,karkipvt.garr)
    tarr = deepcopy(fs_hug_liq.T)
    parr = deepcopy(fs_hug_liq.P)
    sarr = np.zeros(len(fs_hug_liq.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        #print(rho)
        tmp = np.where(karkipvt.rhoarr >= rho)[0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        if len(tmp)>0:
            rind = tmp[0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
            tisen = karkipvt.st4500[rind]
            sisen = karkipvt.refpt[3]
            dS = karkipvt.cv*(np.log(tarr[iii])-np.log(tisen))
            sarr[iii]=sisen+dS
    fs_hug_liq_karki.S = deepcopy(sarr)


# =======================================================================================

# CALCULATE ENTROPY ON THE NEW ALPHA QUARTZ LIQUID HUGONIOT FIT
# dS = dE/T + PdV/T
if 1:
    #
    if 0:
        # set in define_eos.py qtz_p_startliq = 143 # GPa IMPOSED CROSSING INTO LIQUID matches wiggle
        P_start = qtz_p_startliq # GPa based on wiggle
        rho_start = np.interp(P_start,stsmcsmooth.P,stsmcsmooth.rho)
        S_start = np.interp(P_start,stsmcsmooth.P,stsmcsmooth.S_l)
    if 1:
        # calculate the starting point relative to the karkipvt.refpt
        #P_start = qtz_p_startliq # GPa based on wiggle CHECK VALUE
        P_start = qtz_wiggle_P2#qtz_p_startliq # GPa based on wiggle CHECK VALUE
        rho_start = np.interp(P_start,qtz_hug_liq.P,qtz_hug_liq.rho)
        t_start = np.interp(P_start,qtz_hug_liq.P,qtz_hug_liq.T)
        # calculate S from Karki
        rind = np.where(karkipvt.rhoarr >= rho_start)[0][0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        tisen = karkipvt.st4500[rind]
        sisen = karkipvt.refpt[3]
        dS = karkipvt.cv*(np.log(t_start)-np.log(tisen))
        S_start = sisen+dS
    
    ind = np.where(qtz_hug_liq.P >= P_start)[0]
    qtz_hug_liq.S = np.zeros(len(qtz_hug_liq.rho))
    qtz_hug_liq.S[ind[0]]=S_start # J/K/kg
    indref = ind[0]
    varr = 1/(1000.*qtz_hug_liq.rho) # m3/kg
    parr = qtz_hug_liq.P*1.e9 # Pa
    tarr = qtz_hug_liq.T # K
    uarr = qtz_hug_liq.U # J/kg
    for iii in range(1,len(ind)):
        qtz_hug_liq.S[ind[iii]] = qtz_hug_liq.S[ind[iii-1]] + (uarr[ind[iii]]-uarr[ind[iii-1]])/tarr[ind[iii]] + \
            parr[ind[iii]]*(varr[ind[iii]]-varr[ind[iii-1]])/tarr[ind[iii]]
    #ind = np.concatenate((np.where(qtz_hug_liq.rho < rho_start)[0],[indref]))
    ind = np.concatenate((np.where(qtz_hug_liq.P < P_start)[0],[indref]))
    for iii in range(len(ind)-2,-1,-1):
        qtz_hug_liq.S[ind[iii]] = qtz_hug_liq.S[ind[iii+1]] + (uarr[ind[iii]]-uarr[ind[iii+1]])/tarr[ind[iii]] + \
            parr[ind[iii]]*(varr[ind[iii]]-varr[ind[iii+1]])/tarr[ind[iii]]
    print('Density at the crossover point for qtz Hug to liquid = ',np.interp(P_start,qtz_hug_liq.P,qtz_hug_liq.rho))

if 1:
    qtz_hug_liq_karki = deepcopy(qtz_hug_liq)
    # ALPHA QUARTZ LIQUID BRANCH REFERENCE TO KARKI LIQUID MODEL
    # calculate entropies USE LIQUID REFPT everywhere
    # refpt is 4500k; the density arrays are not the same for the Hugoniot and reference isentrope
    # for each rho-T --> calculate dE at that volume from the reference isentrope Kraus point B
    rarr = deepcopy(qtz_hug_liq.rho)
    garr = np.interp(rarr,karkipvt.rhoarr,karkipvt.garr)
    tarr = deepcopy(qtz_hug_liq.T)
    parr = deepcopy(qtz_hug_liq.P)
    sarr = np.zeros(len(qtz_hug_liq.rho))
    for iii in range(0,len(rarr)):
        rho = rarr[iii]
        #print(rho)
        tmp = np.where(karkipvt.rhoarr >= rho)[0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
        if len(tmp)>0:
            rind = tmp[0]  # THIS CAN BE IMPROVED NEED TO INTERPOLATE
            tisen = karkipvt.st4500[rind]
            sisen = karkipvt.refpt[3]
            dS = karkipvt.cv*(np.log(tarr[iii])-np.log(tisen))
            sarr[iii]=sisen+dS
    qtz_hug_liq_karki.S = deepcopy(sarr)


#######################################################################################
## CALCULATE DISORDER REGION QTZ HUGONIOT
#######################################################################################
# Entropy is pegged by the adjacent branches; assumed linear with P
    qtz_hug_glass = deepcopy(qtz_hug_liq)
if 0:
    
    # Use the liquidus to estimate the Hugoniot P-V-Us-up
    v0 = 1./2650.
    #ind = np.where((stsmcsmooth.P > qtz_wiggle_P1*.5) & (stsmcsmooth.P <qtz_p_startliq*1.5))[0]
    #ind = np.where((stsmcsmooth.P > qtz_wiggle_P1) & (stsmcsmooth.P <qtz_p_startliq))[0]
    #rarr = qtz_hug_glass.rho
    #varr = 1/(1000*qtz_hug_glass.rho) # m3/kg
    #parr = np.interp(qtz_hug_glass.rho,stsmcsmooth.rho[ind],stsmcsmooth.P[ind])*1.e9 # Pa
    #uparr = np.sqrt(parr*(v0-varr))
    #usarr = v0*np.sqrt(parr/(v0-varr))
    #rescov = np.polyfit(uparr,usarr,1,cov=True)
    #print(res[0])
    res = [1.30952823, 3.72635724] # old liquidus
    ind2 = np.where(qtz_hug_glass.P <=qtz_p_startliq)[0]
    usarr = qtz_hug_glass.up[ind2]*res[0]+res[1]
    parr = 2.650*usarr*qtz_hug_glass.up[ind2]
    rarr = 1./((1./2.650)*(1-qtz_hug_glass.up[ind2]/usarr))
    qtz_hug_glass.P[ind2] = parr
    qtz_hug_glass.us[ind2] = usarr
    qtz_hug_glass.rho[ind2] = rarr

    # update the shock temperatures with the new us; liquidus fit
    #rest = [  996.66086547, -4653.8528543 ]
    #qtz_hug_glass.T[ind2] = rest[0]*usarr+rest[1] # K and km/s
    # update the entropy on this branch for continuity
#    S_start = np.interp(qtz_usmax,qtz_hug_cacl2.us,qtz_hug_cacl2.S)
#    S_start = np.interp(9.23,qtz_hug_cacl2.us,qtz_hug_cacl2.S)
    S_start = np.interp(qtz_wiggle_P2,stsmcsmooth.P,stsmcsmooth.S_s)
    S_end = np.interp(qtz_us_startliq,qtz_hug_liq.us,qtz_hug_liq.S)
    P_start = np.interp(qtz_usmax,qtz_hug_glass.us,qtz_hug_glass.P)
    P_end = np.interp(qtz_us_startliq,qtz_hug_glass.us,qtz_hug_glass.P)
    qtz_hug_glass.S[ind2] = S_start + (S_end-S_start)/(P_end-P_start)* (qtz_hug_glass.P[ind2]-P_start)
    # do not change T; just map to new bulk S
    #qtz_hug_glass.T[ind2] = np.interp(qtz_hug_glass.P[ind2],stsmcsmooth.P,stsmcsmooth.T)
    #


#######################################################################################
## CALCULATE Coesite Hugoniot
#######################################################################################


##############################
## Coesite
##############################
if 0:
    # Thermal inversion
    cst_wiggle_us1 = 9.5 # km/s
    cst_wiggle_us2 = 10.3 # km/s
    
    # FITS THIS WORK
    # ODR linear solid Us-Up Stv region km/s 
    cst_fit_stv = [1.68055327, 2.59532182] 
    cst_fitcov_stv = [[ 0.01334595, -0.04520109], [-0.04520109,  0.15506528]]
    # ODR linear liquid Us-Up km/s
    cst_fit_liq = [1.40740132, 4.00563031] 
    cst_fitcov_liq = [[ 0.00338923, -0.03182138], [-0.03182138,  0.32098646]]
    # intersection point between solid and liquid Hugoniot segments
    cst_hug_cross = hp.HPeos_Shock()
    cst_hug_cross.up = 5.16309156 # km/s
    cst_hug_cross.us = 11.2721722 # km/s
    cst_hug_cross.P = 168.77783255 # GPa
    cst_hug_cross.rho = 5.3878389 # g/cm3
    # intersection point between solid Hugoniot and solidus
    cst_sol_cross = hp.HPeos_Shock()
    cst_sol_cross.up = 4.34974206 # km/s
    cst_sol_cross.us = 9.90529504 # km/s
    cst_sol_cross.P = 125.80960931 # GPa
    cst_sol_cross.T = 5951.9755829 # K
    cst_sol_cross.rho = 5.20622541 # g/cm3
    # intersection point between liquid Hugoniot and liquidus
    cst_liq_cross = hp.HPeos_Shock()
    cst_liq_cross.P = 179.13720416 # GPa
    cst_liq_cross.T = 6365.61720996 # K
    cst_liq_cross.up = 5.33082122 # km/s
    cst_liq_cross.us = 11.50823513 # km/s
    cst_liq_cross.rho = 5.43982415 # g/cm3
    
    
    # NEW FIT Us-T with melt curve
    # ODR linear
    # cst_tfit_liq = [  1916.72752735, -16244.72022032]
    # ODR quadratic liquid region Feng et al. 2024 Us-T data km/s, K
    # fit from cst_wiggle_us2 to 18 km/s
    cst_tfit_liq =[   56.90163361,   208.76478146, -3572.92059707]
    
    # ODR quadratic solid region this work fit to solid EOS model
    cst_tfit_stv =[  465.12905585, -7068.51146808, 30331.58709626] 
    cst_tfitcov_stv = [[ 2.72842368e+01, -4.82734776e+02,  2.12051066e+03],\
     [-4.82734776e+02,  8.55160232e+03, -3.76147499e+04],\
     [ 2.12051066e+03, -3.76147499e+04,  1.65690889e+05]]
    
    # REFIT WITH MELT CURVE CONSTRAINT
    # liquid region
    cst_hug_liq = hp.HPeos_Shock()
    cst_hug_liq.rho_ref = 2.92 # g/cm3
    cst_hug_liq.up = np.arange(cst_sol_cross.up,18,.01) # km/s
    cst_hug_liq.us = cst_fit_liq[1]+cst_fit_liq[0]*cst_hug_liq.up # km/s
    cst_hug_liq.P = cst_hug_liq.rho_ref*cst_hug_liq.up*cst_hug_liq.us
    cst_hug_liq.rho = 1./((1./cst_hug_liq.rho_ref)*(1.-cst_hug_liq.up/cst_hug_liq.us))
    cst_hug_liq.U = 0.5*cst_hug_liq.P*1.e9*((1/(1000*cst_hug_liq.rho_ref))-(1/(1000*cst_hug_liq.rho))) # J/kg
    # NEW FIT Us-T with melt curve
    # linear
    #cst_hug_liq.T = cst_tfit_liq[1] + (cst_tfit_liq[0]*(cst_hug_liq.us))
    # quadratic
    cst_hug_liq.T = cst_tfit_liq[2] + (cst_tfit_liq[1]*(cst_hug_liq.us))+ (cst_tfit_liq[0]*(cst_hug_liq.us)*(cst_hug_liq.us))
    
    
    # stv region
    cst_hug_stv = hp.HPeos_Shock()
    cst_hug_stv.rho_ref = 2.92 # g/cm3
    cst_hug_stv.up = np.arange(2.86,cst_liq_cross.us,.01) # km/s
    cst_hug_stv.us = cst_fit_stv[1]+cst_fit_stv[0]*cst_hug_stv.up # km/s
    cst_hug_stv.P = cst_hug_stv.rho_ref*cst_hug_stv.up*cst_hug_stv.us
    cst_hug_stv.rho = 1./((1./cst_hug_stv.rho_ref)*(1.-cst_hug_stv.up/cst_hug_stv.us))
    cst_hug_stv.U = 0.5*cst_hug_stv.P*1.e9*((1/(1000*cst_hug_stv.rho_ref))-(1/(1000*cst_hug_stv.rho))) # J/kg
    # NEW FIT Us-T to stv EOS model
    # quadratic
    cst_hug_stv.T = cst_tfit_stv[2] + (cst_tfit_stv[1]*(cst_hug_stv.us))+ (cst_tfit_stv[0]*(cst_hug_stv.us)*(cst_hug_stv.us))
    
    # CST model Hugoniot
    cst_hug_new = hp.HPeos_Shock()
    cst_hug_new.rho_ref = 2.92 # g/cm3
    cst_hug_new.up = np.arange(2.86,18,.01) # km/s
    cst_hug_new.us = cst_fit_liq[1]+cst_fit_liq[0]*cst_hug_new.up # km/s
    cst_hug_new.T = cst_tfit_liq[2] + (cst_tfit_liq[1]*(cst_hug_new.us))+ (cst_tfit_liq[0]*(cst_hug_new.us)*(cst_hug_new.us))
    ind = np.where(cst_hug_new.up < cst_hug_cross.up)
    cst_hug_new.us[ind] = cst_fit_stv[1]+cst_fit_stv[0]*cst_hug_new.up[ind] # km/s
    cst_hug_new.T[ind] = cst_tfit_stv[2] + (cst_tfit_stv[1]*(cst_hug_new.us[ind]))+ (cst_tfit_stv[0]*(cst_hug_new.us[ind])*(cst_hug_new.us[ind]))
    
    
    cst_hug_new.P = cst_hug_new.rho_ref*cst_hug_new.up*cst_hug_new.us
    cst_hug_new.rho = 1./((1./cst_hug_new.rho_ref)*(1.-cst_hug_new.up/cst_hug_new.us))
    cst_hug_new.U = 0.5*cst_hug_new.P*1.e9*((1/(1000*cst_hug_new.rho_ref))-(1/(1000*cst_hug_new.rho))) # J/kg
    
    # follow the melt curve region
    ind = np.where((cst_hug_new.up > cst_sol_cross.up) & (cst_hug_new.up < cst_liq_cross.up))
    cst_hug_new.T[ind] = np.interp(cst_hug_new.P[ind],stsmcsmooth.P,stsmcsmooth.T)

    # new decaying shock pressure calibration - REVISED PRESSURES
    C43939R = deepcopy(C43939)
    C43939R.P = np.interp(C43939R.us,cst_hug_new.us,cst_hug_new.P)
    C43941R = deepcopy(C43941)
    C43941R.P = np.interp(C43941R.us,cst_hug_new.us,cst_hug_new.P)
    C43943R = deepcopy(C43943)
    C43943R.P = np.interp(C43943R.us,cst_hug_new.us,cst_hug_new.P)


#######################################################################################
## CALCULATE MIXED REGIONS OF HUGONIOTS
#######################################################################################
if 0:
    # FUSED SILICA LIQUID TO CACL2
    hugcrossfs = deepcopy(stsmcsmooth)
    
    pind = np.flip(np.where((hugcrossfs.P >= 45) & (hugcrossfs.P <= fs_wiggle_P1))[0])
    # Latent heat is T*dS
    dUstart = hugcrossfs.T[pind[0]]*(hugcrossfs.S_l[pind[0]]-hugcrossfs.S_s[pind[0]])
    print(dUstart/1.e6)
    karkipvt.isenval_liquid(hugcrossfs.S_l[pind[0]])
    liq_uref = np.interp(hugcrossfs.P[pind[0]],karkipvt.isenval_p,karkipvt.isenval_u)
    stvpvt.isenval_stv(hugcrossfs.S_s[pind[0]])
    stv_uref = np.interp(hugcrossfs.P[pind[0]],stvpvt.isenval_p,stvpvt.isenval_u)
    print(liq_uref,stv_uref)
    # U liquid is the new zero = 0.5*P*dV + E0 on the Hugoniot
    hugcrossfs.U_l = np.zeros(len(hugcrossfs.P))
    hugcrossfs.U_s = np.zeros(len(hugcrossfs.P))
    hugrho = np.zeros(len(hugcrossfs.P))
    hugfrac = np.zeros(len(hugcrossfs.P))
    hugcrossfs.U_l[pind[0]] = 0
    hugcrossfs.U_s[pind[0]] = -dUstart
    hugE0 = 0.5* hugcrossfs.P[pind[0]]*1.e9 * (1./(2200)-1/(1000*hugcrossfs.rho[pind[0]]))
    print('Hug E0 = ',hugE0)
    for iii in pind:
        karkipvt.isenval_liquid(hugcrossfs.S_l[iii])
        liq_iii = np.interp(hugcrossfs.P[iii],karkipvt.isenval_p,karkipvt.isenval_u)
        stvpvt.isenval_stv(hugcrossfs.S_s[iii])
        stv_iii = np.interp(hugcrossfs.P[iii],stvpvt.isenval_p,stvpvt.isenval_u)
        hugcrossfs.U_l[iii] = liq_iii-liq_uref
        hugcrossfs.U_s[iii] = stv_iii-stv_uref-dUstart
        sfrac = np.arange(11)/10
        rhonew = np.zeros(11)
        Unew = np.zeros(11)
        Ehug = np.zeros(11)
        for isf in range(11):
            rhonew[isf] = hugcrossfs.rho[iii]*(1-sfrac[isf])+hugcrossfs.rho_s[iii]*(sfrac[isf])
            Unew[isf] = hugcrossfs.U_l[iii]*(1-sfrac[isf])+hugcrossfs.U_s[iii]*(sfrac[isf])
            Ehug[isf] = 0.5*hugcrossfs.P[iii]*1.e9*((1./2200)-1./(1000*rhonew[isf]))-hugE0
        res=hp.Intersection(rhonew,Unew,rhonew,Ehug)
        #print(res)
        if len(res[0])>0:
            hugrho[iii]=res[0][0]
            hugfrac[iii]=np.interp(hugrho[iii],rhonew,sfrac)
            #print(hugcrossfs.P[iii],hugfrac[iii])
            
    hugrhofs = deepcopy(hugrho)
    hugfracfs = deepcopy(hugfrac)
    
    fs_hug_mix = hp.HPeos_Shock()
    ind = np.where(hugrhofs > 0)[0]
    fs_hug_mix.rho = hugrhofs[ind] # g/cm3
    fs_hug_mix.T = stsmcsmooth.T[ind] # K
    fs_hug_mix.P = stsmcsmooth.P[ind] # GPa
    fs_hug_mix.rho0 = 2.2 # g/cm3
    fs_hug_mix.us = (1./(1000*fs_hug_mix.rho0))*np.sqrt(fs_hug_mix.P*1.e9/((1./(1000*fs_hug_mix.rho0))-(1./(1000.*fs_hug_mix.rho))))/1.e3
    fs_hug_mix.up = np.sqrt(fs_hug_mix.P*1.e9*((1./(1000*fs_hug_mix.rho0))-(1./(1000.*fs_hug_mix.rho))))/1.e3
    fs_hug_mix.phasefrac = hugfrac[ind]
    fs_hug_mix.S = stsmcsmooth.S_l[ind]*(1-hugfracfs[ind])+stsmcsmooth.S_s[ind]*(hugfracfs[ind])

    fs_hug_mix1 = deepcopy(fs_hug_mix)
hugrhofs=[]
hugfracfs=[]

if 0:
    # FUSED SILICA LIQUID TO CACL2
    hugcrossfs = deepcopy(stsmcsmooth)
    
    pind = np.flip(np.where((hugcrossfs.P >= 45) & (hugcrossfs.P <= fs_p_startliq))[0])
    # Latent heat is T*dS
    dUstart = hugcrossfs.T[pind[0]]*(hugcrossfs.S_l[pind[0]]-hugcrossfs.S_s[pind[0]])
    print(dUstart/1.e6)
    karkipvt.isenval_liquid(hugcrossfs.S_l[pind[0]])
    liq_uref = np.interp(hugcrossfs.P[pind[0]],karkipvt.isenval_p,karkipvt.isenval_u)
    stvpvt.isenval_stv(hugcrossfs.S_s[pind[0]])
    stv_uref = np.interp(hugcrossfs.P[pind[0]],stvpvt.isenval_p,stvpvt.isenval_u)
    print(liq_uref,stv_uref)
    # U liquid is the new zero = 0.5*P*dV + E0 on the Hugoniot
    hugcrossfs.U_l = np.zeros(len(hugcrossfs.P))
    hugcrossfs.U_s = np.zeros(len(hugcrossfs.P))
    hugrho = np.zeros(len(hugcrossfs.P))
    hugfrac = np.zeros(len(hugcrossfs.P))
    hugcrossfs.U_l[pind[0]] = 0
    hugcrossfs.U_s[pind[0]] = -dUstart
    hugE0 = 0.5* hugcrossfs.P[pind[0]]*1.e9 * (1./(2200)-1/(1000*hugcrossfs.rho[pind[0]]))
    print('Hug E0 = ',hugE0)
    for iii in pind:
        karkipvt.isenval_liquid(hugcrossfs.S_l[iii])
        liq_iii = np.interp(hugcrossfs.P[iii],karkipvt.isenval_p,karkipvt.isenval_u)
        stvpvt.isenval_stv(hugcrossfs.S_s[iii])
        stv_iii = np.interp(hugcrossfs.P[iii],stvpvt.isenval_p,stvpvt.isenval_u)
        hugcrossfs.U_l[iii] = liq_iii-liq_uref
        hugcrossfs.U_s[iii] = stv_iii-stv_uref-dUstart
        sfrac = np.arange(11)/10
        rhonew = np.zeros(11)
        Unew = np.zeros(11)
        Ehug = np.zeros(11)
        for isf in range(11):
            rhonew[isf] = hugcrossfs.rho[iii]*(1-sfrac[isf])+hugcrossfs.rho_s[iii]*(sfrac[isf])
            Unew[isf] = hugcrossfs.U_l[iii]*(1-sfrac[isf])+hugcrossfs.U_s[iii]*(sfrac[isf])
            Ehug[isf] = 0.5*hugcrossfs.P[iii]*1.e9*((1./2200)-1./(1000*rhonew[isf]))-hugE0
        res=hp.Intersection(rhonew,Unew,rhonew,Ehug)
        #print(res)
        if len(res[0])>0:
            hugrho[iii]=res[0][0]
            hugfrac[iii]=np.interp(hugrho[iii],rhonew,sfrac)
            #print(hugcrossfs.P[iii],hugfrac[iii])
            
    hugrhofs = deepcopy(hugrho)
    hugfracfs = deepcopy(hugfrac)
    
    fs_hug_mix = hp.HPeos_Shock()
    ind = np.where(hugrhofs > 0)[0]
    fs_hug_mix.rho = hugrhofs[ind] # g/cm3
    fs_hug_mix.T = stsmcsmooth.T[ind] # K
    fs_hug_mix.P = stsmcsmooth.P[ind] # GPa
    fs_hug_mix.rho0 = 2.2 # g/cm3
    fs_hug_mix.us = (1./(1000*fs_hug_mix.rho0))*np.sqrt(fs_hug_mix.P*1.e9/((1./(1000*fs_hug_mix.rho0))-(1./(1000.*fs_hug_mix.rho))))/1.e3
    fs_hug_mix.up = np.sqrt(fs_hug_mix.P*1.e9*((1./(1000*fs_hug_mix.rho0))-(1./(1000.*fs_hug_mix.rho))))/1.e3
    fs_hug_mix.phasefrac = hugfrac[ind]
    fs_hug_mix.S = stsmcsmooth.S_l[ind]*(1-hugfracfs[ind])+stsmcsmooth.S_s[ind]*(hugfracfs[ind])

    fs_hug_mix2 = deepcopy(fs_hug_mix)
hugrhofs=[]
hugfracfs=[]
fs_hug_mix = []

if 0:
    #--------------------------------------------------------------
    # FUSED SILICA STV TO LIQUID
    hugcrossfs = deepcopy(stsmcsmooth)
    
    pind = np.where((hugcrossfs.P >= fs_wiggle_P2) & (hugcrossfs.P <= fs_wiggle_P1))[0]
    # Latent heat is T*dS
    dUstart = hugcrossfs.T[pind[0]]*(hugcrossfs.S_l[pind[0]]-hugcrossfs.S_s[pind[0]])
    print('dU,Ss,Sl=',dUstart/1.e6,hugcrossfs.S_s[pind[0]],hugcrossfs.S_l[pind[0]])
    karkipvt.isenval_liquid(hugcrossfs.S_l[pind[0]])
    liq_uref = np.interp(hugcrossfs.P[pind[0]],karkipvt.isenval_p,karkipvt.isenval_u)
    stvpvt.isenval_stv(hugcrossfs.S_s[pind[0]])
    stv_uref = np.interp(hugcrossfs.P[pind[0]],stvpvt.isenval_p,stvpvt.isenval_u)
    print(liq_uref,stv_uref)
    # U liquid is the new zero = 0.5*P*dV + E0 on the Hugoniot
    hugcrossfs.U_l = np.zeros(len(hugcrossfs.P))
    hugcrossfs.U_s = np.zeros(len(hugcrossfs.P))
    hugrho = np.zeros(len(hugcrossfs.P))
    hugfrac = np.zeros(len(hugcrossfs.P))
    hugcrossfs.U_l[pind[0]] = dUstart
    hugcrossfs.U_s[pind[0]] = 0
    hugE0 = 0.5* hugcrossfs.P[pind[0]]*1.e9 * (1./(2200)-1/(1000*hugcrossfs.rho_s[pind[0]]))
    print('Hug E0 = ',hugE0)
    for iii in pind[0:2]:
        karkipvt.isenval_liquid(hugcrossfs.S_l[iii])
        liq_iii = np.interp(hugcrossfs.P[iii],karkipvt.isenval_p,karkipvt.isenval_u)
        stvpvt.isenval_stv(hugcrossfs.S_s[iii])
        stv_iii = np.interp(hugcrossfs.P[iii],stvpvt.isenval_p,stvpvt.isenval_u)
        hugcrossfs.U_l[iii] = liq_iii-liq_uref+dUstart
        hugcrossfs.U_s[iii] = stv_iii-stv_uref
        sfrac = np.arange(11)/10
        rhonew = np.zeros(11)
        Unew = np.zeros(11)
        Ehug = np.zeros(11)
        for isf in range(11):
            rhonew[isf] = hugcrossfs.rho_s[iii]*(1-sfrac[isf])+hugcrossfs.rho[iii]*(sfrac[isf])
            Unew[isf] = hugcrossfs.U_s[iii]*(1-sfrac[isf])+hugcrossfs.U_l[iii]*(sfrac[isf])
            Ehug[isf] = 0.5*hugcrossfs.P[iii]*1.e9*((1./2200)-1./(1000*rhonew[isf]))-hugE0
        res=hp.Intersection(rhonew,Unew,rhonew,Ehug)
        print('Ul,Us=',hugcrossfs.U_l[iii],hugcrossfs.U_s[iii]) 
        #print(res)
        if len(res[0])>0:
            hugrho[iii]=res[0][0]
            hugfrac[iii]=np.interp(hugrho[iii],rhonew,sfrac)
            #print('P,frac=',hugcrossfs.P[iii],hugfrac[iii])
            
    hugrhofs = deepcopy(hugrho)
    hugfracfs = deepcopy(hugfrac)
    
    fs_hug_mix2 = hp.HPeos_Shock()
    ind = np.where(hugrhofs > 0)[0]
    fs_hug_mix2.rho = hugrhofs[ind] # g/cm3
    fs_hug_mix2.T = stsmcsmooth.T[ind] # K
    fs_hug_mix2.P = stsmcsmooth.P[ind] # GPa
    fs_hug_mix2.rho0 = 2.2 # g/cm3
    fs_hug_mix2.us = (1./(1000*fs_hug_mix2.rho0))*np.sqrt(fs_hug_mix2.P*1.e9/((1./(1000*fs_hug_mix2.rho0))-(1./(1000.*fs_hug_mix2.rho))))/1.e3
    fs_hug_mix2.up = np.sqrt(fs_hug_mix2.P*1.e9*((1./(1000*fs_hug_mix2.rho0))-(1./(1000.*fs_hug_mix2.rho))))/1.e3
    fs_hug_mix2.phasefrac = hugfrac[ind]
    fs_hug_mix2.S = stsmcsmooth.S_l[ind]*(1-hugfracfs[ind])+stsmcsmooth.S_s[ind]*(hugfracfs[ind])
        
    hugrhofs=[]
    hugfracfs=[]
    fs_hug_mix = []
    
    
    ##-------------------------------------------------------------------------------
    ## QUARTZ LIQUID TO SFT
    
    hugcross = deepcopy(stsmcsmooth)
    
    #pind = np.flip(np.where((hugcross.P > 111) & (hugcross.P < 143))[0])
    pind = np.flip(np.where((hugcross.P > 110) & (hugcross.P < 135))[0])
    # Latent heat is T*dS
    dUstart = hugcross.T[pind[0]]*(hugcross.S_l[pind[0]]-hugcross.S_s[pind[0]])
    print(dUstart/1.e6)
    karkipvt.isenval_liquid(hugcross.S_l[pind[0]])
    liq_uref = np.interp(hugcross.P[pind[0]],karkipvt.isenval_p,karkipvt.isenval_u)
    sftpvt.isenval_sft(hugcross.S_s[pind[0]])
    sft_uref = np.interp(hugcross.P[pind[0]],sftpvt.isenval_p,sftpvt.isenval_u)
    print(liq_uref,sft_uref)
    # U liquid is the new zero = 0.5*P*dV + E0 on the Hugoniot
    hugcross.U_l = np.zeros(len(hugcross.P))
    hugcross.U_s = np.zeros(len(hugcross.P))
    hugrho = np.zeros(len(hugcross.P))
    hugfrac = np.zeros(len(hugcross.P))
    hugcross.U_l[pind[0]] = 0
    hugcross.U_s[pind[0]] = -dUstart
    hugE0 = 0.5* hugcross.P[pind[0]]*1.e9 * (1./(2650)-1/(1000*hugcross.rho[pind[0]]))
    print('Hug E0 = ',hugE0)
    for iii in pind:
        karkipvt.isenval_liquid(hugcross.S_l[iii])
        liq_iii = np.interp(hugcross.P[iii],karkipvt.isenval_p,karkipvt.isenval_u)
        sftpvt.isenval_sft(hugcross.S_s[iii])
        sft_iii = np.interp(hugcross.P[iii],sftpvt.isenval_p,sftpvt.isenval_u)
        hugcross.U_l[iii] = liq_iii-liq_uref
        hugcross.U_s[iii] = sft_iii-sft_uref-dUstart
        sfrac = np.arange(11)/10
        rhonew = np.zeros(11)
        Unew = np.zeros(11)
        Ehug = np.zeros(11)
        for isf in range(11):
            rhonew[isf] = hugcross.rho[iii]*(1-sfrac[isf])+hugcross.rho_s[iii]*(sfrac[isf])
            Unew[isf] = hugcross.U_l[iii]*(1-sfrac[isf])+hugcross.U_s[iii]*(sfrac[isf])
            Ehug[isf] = 0.5*hugcross.P[iii]*1.e9*((1./2650)-1./(1000*rhonew[isf]))-hugE0
        res=hp.Intersection(rhonew,Unew,rhonew,Ehug)
        #print(res)
        if len(res[0])>0:
            hugrho[iii]=res[0][0]
            hugfrac[iii]=np.interp(hugrho[iii],rhonew,sfrac)
            #print(hugcross.P[iii],hugfrac[iii])
    
    qtz_hug_mix = hp.HPeos_Shock()
    ind = np.where(hugrho > 0)[0]
    qtz_hug_mix.rho = hugrho[ind] # g/cm3
    qtz_hug_mix.T = stsmcsmooth.T[ind] # K
    qtz_hug_mix.P = stsmcsmooth.P[ind] # GPa
    qtz_hug_mix.rho0 = 2.65 # g/cm3
    qtz_hug_mix.us = (1./(1000*qtz_hug_mix.rho0))*np.sqrt(qtz_hug_mix.P*1.e9/((1./(1000*qtz_hug_mix.rho0))-(1./(1000.*qtz_hug_mix.rho))))/1.e3
    qtz_hug_mix.up = np.sqrt(qtz_hug_mix.P*1.e9*((1./(1000*qtz_hug_mix.rho0))-(1./(1000.*qtz_hug_mix.rho))))/1.e3
    qtz_hug_mix.phasefrac = hugfrac[ind]
    qtz_hug_mix.S = (stsmcsmooth.S_l[ind])*(1-hugfrac[ind])+(stsmcsmooth.S_s[ind])*(hugfrac[ind])
    
    qtz_hug_mix2= deepcopy(qtz_hug_mix)
    hugrho=[]
    hugfrac=[]
    ##-------------------------------------------------------------------------------
    ## QUARTZ LIQUID TO SFT
    
    hugcross = deepcopy(stsmcsmooth)
    
    #pind = np.flip(np.where((hugcross.P > 111) & (hugcross.P < 143))[0])
    pind = np.flip(np.where((hugcross.P > 90) & (hugcross.P < qtz_p_startliq))[0])
    # Latent heat is T*dS
    dUstart = hugcross.T[pind[0]]*(hugcross.S_l[pind[0]]-hugcross.S_s[pind[0]])
    print(dUstart/1.e6)
    karkipvt.isenval_liquid(hugcross.S_l[pind[0]])
    liq_uref = np.interp(hugcross.P[pind[0]],karkipvt.isenval_p,karkipvt.isenval_u)
    sftpvt.isenval_sft(hugcross.S_s[pind[0]])
    sft_uref = np.interp(hugcross.P[pind[0]],sftpvt.isenval_p,sftpvt.isenval_u)
    print(liq_uref,sft_uref)
    # U liquid is the new zero = 0.5*P*dV + E0 on the Hugoniot
    hugcross.U_l = np.zeros(len(hugcross.P))
    hugcross.U_s = np.zeros(len(hugcross.P))
    hugrho = np.zeros(len(hugcross.P))
    hugfrac = np.zeros(len(hugcross.P))
    hugcross.U_l[pind[0]] = 0
    hugcross.U_s[pind[0]] = -dUstart
    hugE0 = 0.5* hugcross.P[pind[0]]*1.e9 * (1./(2650)-1/(1000*hugcross.rho[pind[0]]))
    print('Hug E0 = ',hugE0)
    for iii in pind:
        karkipvt.isenval_liquid(hugcross.S_l[iii])
        liq_iii = np.interp(hugcross.P[iii],karkipvt.isenval_p,karkipvt.isenval_u)
        sftpvt.isenval_sft(hugcross.S_s[iii])
        sft_iii = np.interp(hugcross.P[iii],sftpvt.isenval_p,sftpvt.isenval_u)
        hugcross.U_l[iii] = liq_iii-liq_uref
        hugcross.U_s[iii] = sft_iii-sft_uref-dUstart
        sfrac = np.arange(11)/10
        rhonew = np.zeros(11)
        Unew = np.zeros(11)
        Ehug = np.zeros(11)
        for isf in range(11):
            rhonew[isf] = hugcross.rho[iii]*(1-sfrac[isf])+hugcross.rho_s[iii]*(sfrac[isf])
            Unew[isf] = hugcross.U_l[iii]*(1-sfrac[isf])+hugcross.U_s[iii]*(sfrac[isf])
            Ehug[isf] = 0.5*hugcross.P[iii]*1.e9*((1./2650)-1./(1000*rhonew[isf]))-hugE0
        res=hp.Intersection(rhonew,Unew,rhonew,Ehug)
        #print(res)
        if len(res[0])>0:
            hugrho[iii]=res[0][0]
            hugfrac[iii]=np.interp(hugrho[iii],rhonew,sfrac)
            #print(hugcross.P[iii],hugfrac[iii])
    
    qtz_hug_mix = hp.HPeos_Shock()
    ind = np.where(hugrho > 0)[0]
    qtz_hug_mix.rho = hugrho[ind] # g/cm3
    qtz_hug_mix.T = stsmcsmooth.T[ind] # K
    qtz_hug_mix.P = stsmcsmooth.P[ind] # GPa
    qtz_hug_mix.rho0 = 2.65 # g/cm3
    qtz_hug_mix.us = (1./(1000*qtz_hug_mix.rho0))*np.sqrt(qtz_hug_mix.P*1.e9/((1./(1000*qtz_hug_mix.rho0))-(1./(1000.*qtz_hug_mix.rho))))/1.e3
    qtz_hug_mix.up = np.sqrt(qtz_hug_mix.P*1.e9*((1./(1000*qtz_hug_mix.rho0))-(1./(1000.*qtz_hug_mix.rho))))/1.e3
    qtz_hug_mix.phasefrac = hugfrac[ind]
    qtz_hug_mix.S = (stsmcsmooth.S_l[ind])*(1-hugfrac[ind])+(stsmcsmooth.S_s[ind])*(hugfrac[ind])
    
    hugrho=[]
    hugfrac=[]
    
    ##-------------------------------------------------------------------------------
    ## QUARTZ CACL2 TO SFT
    if 0:
        # a phase is cacl2; b phase is sft
        # reference zero energy is the phase a at Hugoniot intersecting the phase boundary near Pt3
        hugcross23 = deepcopy(p23bnd)
        # this pind starts at end of CaCl2 Hugoniot and the moves up in pressure along phase boundary
        #pind = np.flip(np.where((hugcross23.P > Pt3) & (hugcross23.P < 120))[0])
        pind = np.flip(np.where((hugcross23.P > 113) & (hugcross23.P < 120))[0])
        #print('pind23 = ',hugcross23.P[pind])
        # Latent heat is T*dS
        dUstart = -hugcross23.T[pind[0]]*(hugcross23.S_b[pind[0]]-hugcross23.S_a[pind[0]])
        print('dustart = ',dUstart/1.e6)
        stvpvt.isenval_stv(hugcross23.S_a[pind[0]])
        stv_uref = np.interp(hugcross23.P[pind[0]],stvpvt.isenval_p,stvpvt.isenval_u)
        sftpvt.isenval_sft(hugcross23.S_b[pind[0]])
        sft_uref = np.interp(hugcross23.P[pind[0]],sftpvt.isenval_p,sftpvt.isenval_u)
        print('sft ref, stv ref = ',sft_uref,stv_uref)
        # U stv is the new zero = 0.5*P*dV + E0 on the Hugoniot
        hugcross23.U_a = np.zeros(len(hugcross23.P))
        hugcross23.U_b = np.zeros(len(hugcross23.P))
        hugrho23 = np.zeros(len(hugcross23.P))
        hugfrac23 = np.zeros(len(hugcross23.P))
        hugcross23.U_a[pind[0]] = 0
        hugcross23.U_b[pind[0]] = -dUstart
        hugE0 = 0.5* hugcross23.P[pind[0]]*1.e9 * (1./(2650)-1/(1000*hugcross23.rho_a[pind[0]]))
        print('Hug E0 = ',hugE0)
        for iii in pind:
            stvpvt.isenval_stv(hugcross23.S_a[iii])
            stv_iii = np.interp(hugcross23.P[iii],stvpvt.isenval_p,stvpvt.isenval_u)
            sftpvt.isenval_sft(hugcross23.S_b[iii])
            sft_iii = np.interp(hugcross23.P[iii],sftpvt.isenval_p,sftpvt.isenval_u)
            hugcross23.U_a[iii] = stv_iii-stv_uref
            hugcross23.U_b[iii] = sft_iii-sft_uref-dUstart
            sfrac = np.arange(11)/10 # 0 to 1, weights phase b
            rhonew = np.zeros(11)
            Unew = np.zeros(11)
            Ehug = np.zeros(11)
            for isf in range(11):
                rhonew[isf] = hugcross23.rho_a[iii]*(1-sfrac[isf])+hugcross23.rho_b[iii]*(sfrac[isf])
                Unew[isf] = hugcross23.U_a[iii]*(1-sfrac[isf])+hugcross23.U_b[iii]*(sfrac[isf])
                Ehug[isf] = 0.5*hugcross23.P[iii]*1.e9*((1./2650)-1./(1000*rhonew[isf]))-hugE0
            res=hp.Intersection(rhonew,Unew,rhonew,Ehug)
            if len(res[0])>0:
                #print(res)
                hugrho23[iii]=res[0][0]
                hugfrac23[iii]=np.interp(hugrho23[iii],rhonew,sfrac)
                print('mix23 = ',hugcross23.P[iii],hugfrac23[iii],hugcross23.T[iii])
        
        qtz_hug_mix23 = hp.HPeos_Shock()
        ind = np.where(hugrho23 > 0)[0]
        qtz_hug_mix23.rho = hugrho23[ind] # g/cm3
        qtz_hug_mix23.T = hugcross23.T[ind] # K # flipped because of initial index flip
        qtz_hug_mix23.P = hugcross23.P[ind] # GPa # flipped because of initial index flip
        qtz_hug_mix23.rho0 = 2.65 # g/cm3
        qtz_hug_mix23.us = (1./(1000*qtz_hug_mix23.rho0))*np.sqrt(qtz_hug_mix23.P*1.e9/((1./(1000*qtz_hug_mix23.rho0))-(1./(1000.*qtz_hug_mix23.rho))))/1.e3
        qtz_hug_mix23.up = np.sqrt(qtz_hug_mix23.P*1.e9*((1./(1000*qtz_hug_mix23.rho0))-(1./(1000.*qtz_hug_mix23.rho))))/1.e3
        qtz_hug_mix23.phasefrac = hugfrac23[ind] # fraction of phase b
        # do not need to flip entropy as long of frac is associated wiht phase b
        qtz_hug_mix23.S = hugcross23.S_a[ind]*(1-qtz_hug_mix23.phasefrac)+hugcross23.S_b[ind]*(qtz_hug_mix23.phasefrac)
        
        hugrho23=[]
        hugfrac23=[]
    
    #######################################################################################
    # Pressure bounds on the Omega data wiggles
    fs_wiggle_P1mix = np.interp(fs_usmax,fs_hug_mix1.us,fs_hug_mix1.P)
    fs_wiggle_P2mix = np.interp(fs_usmin,fs_hug_mix1.us,fs_hug_mix1.P)
    #######################################################################################

if 0:
    # Fused Silica Hugoniot Liquid Region
    # cubic from Root et al. 2018
    fs_hug_R19c = hp.HPeos_Shock()
    fs_hug_R19c.rho_ref = 2.204 # g/cm3
    fs_hug_R19c.up = np.arange(3,20,.005) # km/s
    fs_hug_R19c.us = 1.386+1.647*fs_hug_R19c.up-1.146e-2*fs_hug_R19c.up*fs_hug_R19c.up-0.6052e-4*fs_hug_R19c.up*fs_hug_R19c.up*fs_hug_R19c.up # km/s
    fs_hug_R19c.P = fs_hug_R19c.rho_ref*fs_hug_R19c.up*fs_hug_R19c.us
    fs_hug_R19c.rho = 1./((1./fs_hug_R19c.rho_ref)*(1.-fs_hug_R19c.up/fs_hug_R19c.us))
    
    # NEED TEMPERATURE REFIT
    # using Kaitlyn fit to the temperatures boslough + hicks -- only a small range in pressure; melting to XXX GPa 
    fs_hug_R19c.T = -1.99167591e+04 + (6.37525746e+03*(fs_hug_R19c.us)) + (-6.71164718e+02*(fs_hug_R19c.us)**2) + ((3.64336344e+01)*(fs_hug_R19c.us)**3) + ((-5.41405759e-01)*(fs_hug_R19c.us)**4)
    #
    # HICKS CALIBRATION FS
    # Calculate the Hicks Hugoniot points using our refit us-up
    #fs_hug_Hicks06.P = np.interp(fs_hug_Hicks06.us,fs_hug_liq.us,fs_hug_liq.P) # GPa
    #fs_hug_Hicks06.up = np.interp(fs_hug_Hicks06.us,fs_hug_liq.us,fs_hug_liq.up) # km/s
    #fs_hug_Hicks06.rho = np.interp(fs_hug_Hicks06.us,fs_hug_liq.us,fs_hug_liq.rho) # g/cm3
#
if 0:
    
    if qtz_liq_hug == 'NEW':
        # cubic fit STS May 7, 2025
        qtz_liq_usupfit = [ 4.97302828e-04, -3.17517552e-02,  1.84963433e+00,  1.75574087e+00]
        res=qtz_liq_usupfit
        qtz_hug_liq.us = res[3]+res[2]*qtz_hug_liq.up+res[1]*(qtz_hug_liq.up**2)+res[0]*(qtz_hug_liq.up**3)
    
    if qtz_liq_hug == 'KD13':
        # cubic fit from Knudson and Desjarlais PRB 2013
        # NOTE FIT INCLUDES STRANGE LOW TEMP HIGH ENTROPY PHASE DATA
        qtz_liq_usupfit = [1.754,1.862,-3.364e-2,5.666e-4]
        res=qtz_liq_usupfit
        qtz_hug_liq.us = res[0]+res[1]*qtz_hug_liq.up+res[2]*(qtz_hug_liq.up**2)+res[3]*(qtz_hug_liq.up**3)
    if qtz_liq_hug == 'KD09': 
        # exponential fit from Knudson and Desjarlais PRL 2009
        a = 6.26 # km/s
        b = 1.2
        c = 2.56
        d=0.37 # 1/(km/s)
        qtz_hug_liq.us   = (a + b*qtz_hug_liq.up - c*(qtz_hug_liq.up)*np.exp(-d*(qtz_hug_liq.up))) # km/s
    
    qtz_hug_liq.P = qtz_hug_liq.rho_ref*qtz_hug_liq.up*qtz_hug_liq.us
    qtz_hug_liq.rho = 1/((1/qtz_hug_liq.rho_ref)*(1-qtz_hug_liq.up/qtz_hug_liq.us))
    qtz_hug_liq.U = 0.5*qtz_hug_liq.P*1.e9*((1/(1000*qtz_hug_liq.rho_ref))-(1/(1000*qtz_hug_liq.rho))) # J/kg
    # temperature fit from Davies and Millot Omega data T = cubic in (Us km/s)
    qtz_hug_liq_tfit = [    3.51454387,   -17.48937584 ,  515.96055948, -1583.40664638]
    tres = qtz_hug_liq_tfit
    qtz_hug_liq.T = tres[3]+tres[2]*qtz_hug_liq.us+tres[1]*qtz_hug_liq.us**2+tres[0]*(qtz_hug_liq.us**3)
    
    # HICKS CALIBRATION FS
    # Calculate the Hicks Hugoniot points using our refit us-up
    #qtz_hug_Hicks06.P = np.interp(qtz_hug_Hicks06.us,qtz_hug_liq.us,qtz_hug_liq.P) # GPa
    #qtz_hug_Hicks06.up = np.interp(qtz_hug_Hicks06.us,qtz_hug_liq.us,qtz_hug_liq.up) # km/s
    #qtz_hug_Hicks06.rho = np.interp(qtz_hug_Hicks06.us,qtz_hug_liq.us,qtz_hug_liq.rho) # g/cm3
