
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

if 1:
    # calculate entropies
    # refpt is 4500 K; the density arrays are not the same for the Hugoniot and reference isentrope
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
        sarr[iii]=sisen+dS
    fs_hug_stv.S = deepcopy(sarr)

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
