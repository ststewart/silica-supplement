# load_solid_models.py
# Vinet parameters stv and cacl2 and pbo2 are in define_eos.py
# specific heat capacity prefactor solid_cvfactor is in define_eos.py
# liquid EOS model is assumed to be defined to generate EOS with the same density array
#
# quartz Hugoniot parameters set in define_eos.py
# qtz_cacl2_usupfit, qtz_stv_usupfit
#
#----------------
############################################################################
# reference curves for the model
# 4000 K isotherms from Driver et al. PNAS 2010
# Driver 4000 K isotherms
dat_file   = silicadatadir+'stvcacl2.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
vol=np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])/(np.power(1.8897161646321,3))
dat.rho      = sio2mol/vol*1.6605391 # amu/a^3 -> g/cm3
dat.rho_units = ['g/cm3']
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
#print(dat.P)
#print(dat.rho)
stvdat = deepcopy(dat)
# adjust density
stvdat.rho = stvdat.rho+stv_density_shift

#dat_file   = silicadatadir+'pbo2.txt'
dat_file   = silicadatadir+'pbo2-v2.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
vol=np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])/(np.power(1.8897161646321,3))
dat.rho      = sio2mol/vol*1.6605391 # amu/a^3 -> g/cm3
dat.rho_units = ['g/cm3']
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
#print(dat.P)
#print(dat.rho)
sftdat = deepcopy(dat)
# adjust density
sftdat.rho = sftdat.rho+sft_density_shift



# stishovite-CaCl2
stvpvt = hp.HPeos_Model()
stvpvt.rhoarr = deepcopy(karkipvt.rhoarr)
stvpvt.varr = 1/stvpvt.rhoarr
stvpvt.g0 = 1.22
stvpvt.q  = 2.22
# 1 Bohr radius = 0.529177249 angstrom
# 1.6605391 # amu/a^3 -> g/cm3
v0dft= 159 # bohr^3/sio2
v0 = v0dft/sio2mol*np.power(0.529177249,3)/1.6605391/1000.
stvpvt.rhoref = 1/v0/1.e3+stv_density_shift
stvpvt.garr = stvpvt.g0*np.power((stvpvt.varr*stvpvt.rhoref),stvpvt.q)

#print('stvrhoref =', stvpvt.rhoref)
# cv = stv_cvfactor*nR
stvpvt.cv = stv_cvfactor*3*8.314/sio2mol/1.e-3 # J/kg/mol * mol/g

# interpolate 4000 K isotherm
stvpvt.p4000 = np.zeros(len(stvpvt.rhoarr))
stvpvt.t4000 = np.zeros(len(stvpvt.rhoarr))+4000.
for iii in range(len(stvpvt.rhoarr)):
    stvpvt.p4000[iii] = np.interp(stvpvt.rhoarr[iii],np.flip(stvdat.rho),np.flip(stvdat.P))    


stvpvt.p0000 = np.zeros(len(stvpvt.rhoarr))
par = stv
rho0 = sio2mol*1.6605391/par[0]
eta = np.power(rho0/stvpvt.rhoarr,0.333)
parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
ind = np.where(parr <= 50)[0]
stvpvt.p0000[ind]=parr[ind]
par = cacl2
rho0 = sio2mol*1.6605391/par[0]
stvpvt.rho00 = rho0 # 0K reference density
eta = np.power(rho0/stvpvt.rhoarr,0.333)
parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
#ind = np.where(parr > 50)[0]
ind = np.where(parr > 49)[0]
stvpvt.p0000[ind]=parr[ind]

par=[]
parr=[]

# siefertite
sftpvt = hp.HPeos_Model()
sftpvt.rhoarr = deepcopy(karkipvt.rhoarr)
sftpvt.varr = 1/sftpvt.rhoarr
sftpvt.g0 = 1.27
sftpvt.q  = 2.05
# 1 Bohr radius = 0.529177249 angstrom
# 1.6605391 # amu/a^3 -> g/cm3
v0dft= 154.8 # bohr^3/sio2
v0 = v0dft/sio2mol*np.power(0.529177249,3)/1.6605391/1000.
sftpvt.rhoref = 1/v0/1.e3+sft_density_shift
sftpvt.garr = sftpvt.g0*np.power((sftpvt.varr*sftpvt.rhoref),sftpvt.q)

# cv = sft_cvfactor*nR
sftpvt.cv = sft_cvfactor*3*8.314/sio2mol/1.e-3 # J/kg/mol * mol/g

# interpolate 4000 K isotherm
sftpvt.p4000 = np.zeros(len(sftpvt.rhoarr))
sftpvt.t4000 = np.zeros(len(sftpvt.rhoarr))+4000.
for iii in range(len(sftpvt.rhoarr)):
    sftpvt.p4000[iii] = np.interp(sftpvt.rhoarr[iii],sftdat.rho,sftdat.P)    

sftpvt.p0000 = np.zeros(len(sftpvt.rhoarr))
par = pbo2
rho0 = sio2mol*1.6605391/par[0]
sftpvt.rho00 = rho0 # 0K reference density
eta = np.power(rho0/sftpvt.rhoarr,0.333)
parr = 3.*par[1]*((1-eta)/eta/eta)*np.exp(3./2.*(par[2]-1)*(1-eta))
sftpvt.p0000=parr

#print('sftrhoref =', sftpvt.rhoref)

par=[]
parr=[]

############################################################################
# Calculate a set of isotherms and isentropes to match the liquid model

stvpvt.stempsarr = np.asarray([3000,3130,3500,4000,4500,4800,5000,5075,5500,5900,6000,6100,6500,7000])
for iii in range(len(stvpvt.stempsarr)):
    temp = stvpvt.stempsarr[iii]
    dE = stvpvt.cv*(temp-4000) # 4000 is the reference curve # mks
    pth = dE*stvpvt.garr*(stvpvt.rhoarr*1000)/1.e9 # Pa->GPa
    #ind = np.where(stvpvt.rhoarr < stvpvt.rhoref)[0]
    #pth[ind]=0
    if temp == 3000:
        stvpvt.p3000 = stvpvt.p4000+pth
    if temp == 3130:
        stvpvt.p3130 = stvpvt.p4000+pth
    if temp == 3500:
        stvpvt.p3500 = stvpvt.p4000+pth
    if temp == 4000:
        stvpvt.p4000 = stvpvt.p4000+pth
    if temp == 4500:
        stvpvt.p4500 = stvpvt.p4000+pth
    if temp == 4800:
        stvpvt.p4800 = stvpvt.p4000+pth
    if temp == 5000:
        stvpvt.p5000 = stvpvt.p4000+pth
    if temp == 5075:
        stvpvt.p5075 = stvpvt.p4000+pth
    if temp == 5500:
        stvpvt.p5500 = stvpvt.p4000+pth
    if temp == 5900:
        stvpvt.p5900 = stvpvt.p4000+pth
    if temp == 6000:
        stvpvt.p6000 = stvpvt.p4000+pth
    if temp == 6100:
        stvpvt.p6100 = stvpvt.p4000+pth
    if temp == 6500:
        stvpvt.p6500 = stvpvt.p4000+pth
    if temp == 7000:
        stvpvt.p7000 = stvpvt.p4000+pth

# stv_reference_point = [-1,70,4500,3645] # (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg)       
rhofit = np.interp(stv_reference_point[1],stvpvt.p4500,stvpvt.rhoarr)
stv_reference_point[0]=rhofit
#print('stv reference density is (g/cm3)= ',rhofit)
stvpvt.refpt = np.asarray(stv_reference_point)

# calculate entropies
for iii in range(len(stvpvt.stempsarr)):
    # set conditions for this isentrope
    rhoref = stvpvt.refpt[0]
    rhoarr = deepcopy(stvpvt.rhoarr)
    garr = deepcopy(stvpvt.garr)
    tarr = np.zeros(len(rhoarr))
    parr = np.zeros(len(rhoarr))
    t0 = stvpvt.stempsarr[iii] # temp at rhoref

    ind = np.where(rhoarr > rhoref)[0]
    #print(iii,rhoref,ind[0])
    tarr[ind[0]]=t0
    for iii in range(ind[1],np.max(ind)+1):
        a = -garr[iii]*(np.log(1/rhoarr[iii])-np.log(1/rhoarr[iii-1]))
        lnt = a+np.log(tarr[iii-1])
        tarr[iii]=np.exp(lnt)
        #print(iii,a,tarr[iii-1],tarr[iii])
    rind = np.flip(np.where(rhoarr < rhoref)[0])
    for iii in rind:
        a = -garr[iii]*(np.log(1/rhoarr[iii])-np.log(1/rhoarr[iii+1]))
        lnt = a+np.log(tarr[iii+1])
        tarr[iii]=np.exp(lnt)
        #print(iii,a,tarr[iii+1],tarr[iii])
    # for each stv/cacl2 isentrope, calculate the pressure
    dE = stvpvt.cv*(tarr-4500) # 4500 is the reference curve FOR ENTROPY # mks
    parr = stvpvt.p4500+ dE*stvpvt.garr*(stvpvt.rhoarr*1000)/1.e9 # Pa->GPa   
    
    if t0 == 3000:
        stvpvt.st3000=deepcopy(tarr)
        stvpvt.sp3000=deepcopy(parr)
    if t0 == 3130:
        stvpvt.st3130=deepcopy(tarr)
        stvpvt.sp3130=deepcopy(parr)
    if t0 == 3500:
        stvpvt.st3500=deepcopy(tarr)
        stvpvt.sp3500=deepcopy(parr)
    if t0 == 4000:
        stvpvt.st4000=deepcopy(tarr)
        stvpvt.sp4000=deepcopy(parr)
    if t0 == 4500:
        stvpvt.st4500=deepcopy(tarr)
        stvpvt.sp4500=deepcopy(parr)
    if t0 == 4800:
        stvpvt.st4800=deepcopy(tarr)
        stvpvt.sp4800=deepcopy(parr)
    if t0 == 5000:
        stvpvt.st5000=deepcopy(tarr)
        stvpvt.sp5000=deepcopy(parr)
    if t0 == 5075:
        stvpvt.st5075=deepcopy(tarr)
        stvpvt.sp5075=deepcopy(parr)
    if t0 == 5500:
        stvpvt.st5500=deepcopy(tarr)
        stvpvt.sp5500=deepcopy(parr)
    if t0 == 5900:
        stvpvt.st5900=deepcopy(tarr)
        stvpvt.sp5900=deepcopy(parr)
    if t0 == 6000:
        stvpvt.st6000=deepcopy(tarr)
        stvpvt.sp6000=deepcopy(parr)
    if t0 == 6100:
        stvpvt.st6100=deepcopy(tarr)
        stvpvt.sp6100=deepcopy(parr)
    if t0 == 6500:
        stvpvt.st6500=deepcopy(tarr)
        stvpvt.sp6500=deepcopy(parr)
    if t0 == 7000:
        stvpvt.st7000=deepcopy(tarr)
        stvpvt.sp7000=deepcopy(parr)
    tarr=[]
    parr=[]

stvpvt.svalsarr = np.zeros(len(stvpvt.stempsarr))
for iii in range(len(stvpvt.stempsarr)):
    stvpvt.svalsarr[iii] = (np.log(stvpvt.stempsarr[iii])-np.log(stv_reference_point[2]))*stvpvt.cv+stvpvt.refpt[3]
    #print(iii,stvpvt.stempsarr[iii],stvpvt.svalsarr[iii])

############################################################################
# PbO2 Siefertite isotherms calculation

sftpvt.stempsarr = np.asarray([3000,3500,4000,4500,4800,5000,5500,5900,6000,6100,6500,7000])
for iii in range(len(sftpvt.stempsarr)):
    temp = sftpvt.stempsarr[iii]
    dE = sftpvt.cv*(temp-4000) # 4000 is the reference curve for ISOTHERMS # mks
    pth = dE*sftpvt.garr*(sftpvt.rhoarr*1000)/1.e9 # Pa->GPa
    ind = np.where(sftpvt.rhoarr < sftpvt.rhoref)[0]
    pth[ind] = 0
    if temp == 3000:
        sftpvt.p3000 = sftpvt.p4000+pth
    if temp == 3500:
        sftpvt.p3500 = sftpvt.p4000+pth
    if temp == 4000:
        sftpvt.p4000 = sftpvt.p4000+pth
    if temp == 4500:
        sftpvt.p4500 = sftpvt.p4000+pth
    if temp == 4800:
        sftpvt.p4800 = sftpvt.p4000+pth
    if temp == 5000:
        sftpvt.p5000 = sftpvt.p4000+pth
    if temp == 5500:
        sftpvt.p5500 = sftpvt.p4000+pth
    if temp == 5900:
        sftpvt.p5900 = sftpvt.p4000+pth
    if temp == 6000:
        sftpvt.p6000 = sftpvt.p4000+pth
    if temp == 6100:
        sftpvt.p6100 = sftpvt.p4000+pth
    if temp == 6500:
        sftpvt.p6500 = sftpvt.p4000+pth
    if temp == 7000:
        sftpvt.p7000 = sftpvt.p4000+pth

#####################################################################
## calculate reference points between CaCl2 and PbO2 and entropy tie point
##
## values needed for PV phase diagram
# find tie points for phase boundary CaCl2 to seifertite
rho0_23_sft = np.interp(Pt3_0k,sftpvt.p0000,sftpvt.rhoarr)
rhos_23_sft = np.interp(Pt3,sftpvt.p6100,sftpvt.rhoarr)
print(rho0_23_sft,rhos_23_sft)
rho0_23_stv = np.interp(Pt3_0k,stvpvt.p0000,stvpvt.rhoarr)
rhos_23_stv = np.interp(Pt3,stvpvt.p6100,stvpvt.rhoarr)
print(rho0_23_stv,rhos_23_stv)
# find tie points for phase boundary stishovite to CaCl2
rho0_12_stv = np.interp(Pt2_0k,stvpvt.p0000,stvpvt.rhoarr)
rhos_12_stv = np.interp(Pt2,stvpvt.p5075,stvpvt.rhoarr)
#print(rho0_12_stv,rhos_12_stv)
# stishovite to coesite tie point
rho0_01_stv = np.interp(Pt1_0k,stvpvt.p0000,stvpvt.rhoarr)
rhos_01_stv = np.interp(Pt1,stvpvt.p3130,stvpvt.rhoarr)

#####################################################################
# CALCULATE THE SOLIDUS
## calculate entropy for PbO2

# Siefertite tie point based on the triple point with CaCl2 T=5900 K TRIPLE POINT 
Ss_23_sft = np.interp(Pt3,stsmcsmooth.P,stsmcsmooth.S_l)
sft_reference_point[0] = rhos_23_sft # g/cm3
sft_reference_point[3] = Ss_23_sft+sft_dS # J/K/kg
sftpvt.refpt = sft_reference_point #  (rho,p,t,s)
print('SFT REFPT = ',sftpvt.refpt)
#print(sftpvt.stempsarr)
# calculate entropies
for iii in range(len(sftpvt.stempsarr)):
    # set conditions for this isentrope
    rhoref = sftpvt.refpt[0]
    rhoarr = deepcopy(sftpvt.rhoarr)
    garr = deepcopy(sftpvt.garr)
    tarr = np.zeros(len(rhoarr))
    t0 = sftpvt.stempsarr[iii] # temp at rhoref
    
    ind = np.where(rhoarr > rhoref)[0]
    #print(iii,rhoref,t0,ind[0])
    tarr[ind[0]]=t0
    for iii in range(ind[1],np.max(ind)+1):
        a = -garr[iii]*(np.log(1/rhoarr[iii])-np.log(1/rhoarr[iii-1]))
        lnt = a+np.log(tarr[iii-1])
        tarr[iii]=np.exp(lnt)
        #print(iii,a,tarr[iii-1],tarr[iii])
    rind = np.flip(np.where(rhoarr < rhoref)[0])
    for iii in rind:
        a = -garr[iii]*(np.log(1/rhoarr[iii])-np.log(1/rhoarr[iii+1]))
        lnt = a+np.log(tarr[iii+1])
        tarr[iii]=np.exp(lnt)
        #print(iii,a,tarr[iii+1],tarr[iii])
    
    # for each pbo2 isentrope, calculate the pressure
    dE = sftpvt.cv*(tarr-sftpvt.refpt[2]) # Temp on ref curve
#    parr = sftpvt.p5900+ dE*sftpvt.garr*(sftpvt.rhoarr*1000)/1.e9 # Pa->GPa   #5900 K is the reference for ENTROPY
    parr = sftpvt.p6100+ dE*sftpvt.garr*(sftpvt.rhoarr*1000)/1.e9 # Pa->GPa   #6100 K is the reference for ENTROPY
#    parr = sftpvt.p6000+ dE*sftpvt.garr*(sftpvt.rhoarr*1000)/1.e9 # Pa->GPa   #6000 K is the reference for ENTROPY

    if t0 == 3000:
        sftpvt.st3000=deepcopy(tarr)
        sftpvt.sp3000=deepcopy(parr)
    if t0 == 3500:
        sftpvt.st3500=deepcopy(tarr)
        sftpvt.sp3500=deepcopy(parr)
    if t0 == 4000:
        sftpvt.st4000=deepcopy(tarr)
        sftpvt.sp4000=deepcopy(parr)
    if t0 == 4500:
        sftpvt.st4500=deepcopy(tarr)
        sftpvt.sp4500=deepcopy(parr)
        #print(tarr)
    if t0 == 4800:
        sftpvt.st4800=deepcopy(tarr)
        sftpvt.sp4800=deepcopy(parr)
    if t0 == 5000:
        sftpvt.st5000=deepcopy(tarr)
        sftpvt.sp5000=deepcopy(parr)
    if t0 == 5500:
        sftpvt.st5500=deepcopy(tarr)
        sftpvt.sp5500=deepcopy(parr)
    if t0 == 5900:
        sftpvt.st5900=deepcopy(tarr)
        sftpvt.sp5900=deepcopy(parr)
    if t0 == 6000:
        sftpvt.st6000=deepcopy(tarr)
        sftpvt.sp6000=deepcopy(parr)
    if t0 == 6100:
        sftpvt.st6100=deepcopy(tarr)
        sftpvt.sp6100=deepcopy(parr)
    if t0 == 6500:
        sftpvt.st6500=deepcopy(tarr)
        sftpvt.sp6500=deepcopy(parr)
    if t0 == 7000:
        sftpvt.st7000=deepcopy(tarr)
        sftpvt.sp7000=deepcopy(parr)
    tarr=[]
    parr=[]

sftpvt.svalsarr = np.zeros(len(sftpvt.stempsarr))
for iii in range(len(sftpvt.stempsarr)):
    sftpvt.svalsarr[iii] = (np.log(sftpvt.stempsarr[iii])-np.log(sftpvt.refpt[2]))*sftpvt.cv+sftpvt.refpt[3]
    #print(iii,sftpvt.stempsarr[iii],sftpvt.svalsarr[iii])

#####################################################################
# CALCULATE THE SOLIDUS
#
# interpolate the stv-cacl2 isotherms for density on the melt curve
# then use isentropes to determine the entropy on the melt curve

default_drho = 1.01
stsmcsmooth.rho_s = stsmcsmooth.rho*default_drho # initial shift from liquidus at high pressures
stsmcsmooth.rho_s_low = stsmcsmooth.rho*default_drho*.9975 # initial shift from liquidus at high pressures
stsmcsmooth.rho_s_high = stsmcsmooth.rho*default_drho*1.0025 # initial shift from liquidus at high pressures
#stsmcsmooth.S_s = stsmcsmooth.S_l- 1.414*Rgas/(sio2mol/1000)-229-80+60
#stsmcsmooth.S_s = stsmcsmooth.S_l + sft_dS
stsmcsmooth.S_s = np.zeros(len(stsmcsmooth.S_l))

##stsmcsmooth
# use the Driver model to calculate the solidus densities for SiO2 based on stsmcsmooth array
# Stv and CaCl2 region for 14-Pt3 GPa
rhostvsmooth = np.zeros(len(stsmcsmooth.P))
entstvsmooth = np.zeros(len(stsmcsmooth.P))



for iii in range(0,len(stsmcsmooth.P)):
    # have P and T, calculate isotherm
    stvpvt.isotherm_stv(stsmcsmooth.T[iii])
    rhostvsmooth[iii] = np.interp(stsmcsmooth.P[iii],stvpvt.isotherm_p,stvpvt.rhoarr)
    entstvsmooth[iii] = np.interp(stsmcsmooth.P[iii],stvpvt.isotherm_p,stvpvt.isotherm_s)        
    if (stsmcsmooth.P[iii] <= Pt3) & (stsmcsmooth.P[iii] >= Pt1):
        stsmcsmooth.rho_s[iii]=rhostvsmooth[iii] 
        stsmcsmooth.S_s[iii]=entstvsmooth[iii]
    if ((stsmcsmooth.P[iii]< Pt1)):
        stsmcsmooth.rho_s[iii]=stsmcsmooth.rho[iii]*1.01 # Coesite/qtz - cst-L negative slope near stv TP - just offset slightly for plotting purposes Zhang et al. 1993
        stsmcsmooth.S_s[iii]=stsmcsmooth.S_l[iii]-300.
            
for iii in range(0,len(stsmcsmooth.P)):
    # have P and T, calculate isotherm
    stvpvt.isotherm_stv(stsmcsmooth.T[iii]+stsmcsmooth.T_std[iii])
    rhostvsmooth[iii] = np.interp(stsmcsmooth.P[iii],stvpvt.isotherm_p,stvpvt.rhoarr)
    if stsmcsmooth.P[iii] <= Pt3:
        stsmcsmooth.rho_s_low[iii]=rhostvsmooth[iii] 
for iii in range(0,len(stsmcsmooth.P)):
    # have P and T, calculate isotherm
    stvpvt.isotherm_stv(stsmcsmooth.T[iii]-stsmcsmooth.T_std[iii])
    rhostvsmooth[iii] = np.interp(stsmcsmooth.P[iii],stvpvt.isotherm_p,stvpvt.rhoarr)
    if stsmcsmooth.P[iii] <= Pt3:
        stsmcsmooth.rho_s_high[iii]=rhostvsmooth[iii] 


# interpolate the stv-cacl2 isotherms for density on the melt curve
# then use isentropes to determine the entropy on the melt curve

# Use Park et al. 2025 for Pyrite phase boundary
parkp = np.asarray([140,160,180,200,220,240])
parkdv = (np.asarray([.319,.276,.243,.215,.191,.172])*1.6605391/sio2mol)

##stsmcsmooth
# use the Driver model to calculate the solidus densities for SiO2 based on stsmcsmooth array
# PbO2 for P>Pt3 GPa
## uncertain transition to pyrite
## limit the change in density to 2%
rhosftsmooth = np.zeros(len(stsmcsmooth.P))
entsftsmooth = np.zeros(len(stsmcsmooth.P))
for iii in range(0,len(stsmcsmooth.P)):
    # have P and T, calculate isotherm
    sftpvt.isotherm_sft(stsmcsmooth.T[iii])
    rhosftsmooth[iii] = np.interp(stsmcsmooth.P[iii],sftpvt.isotherm_p,sftpvt.rhoarr)
    entsftsmooth[iii] = np.interp(stsmcsmooth.P[iii],sftpvt.isotherm_p,sftpvt.isotherm_s)
    if (stsmcsmooth.P[iii]> Pt3) & (stsmcsmooth.P[iii]<= Pt4) & (rhosftsmooth[iii] >= stsmcsmooth.rho[iii]*default_drho):
#    if (stsmcsmooth.P[iii]> Pt3) & (stsmcsmooth.P[iii]<= Pt4) :
        # sft region
        stsmcsmooth.rho_s[iii]=rhosftsmooth[iii] 
#        stsmcsmooth.S_s[iii]=entsftsmooth[iii] 
        stsmcsmooth.S_s[iii]=stsmcsmooth.S_l[iii]+sft_dS+4*(stsmcsmooth.P[iii]-165)

    if ((stsmcsmooth.P[iii]>= 165)&(stsmcsmooth.P[iii]<=240)):
        # pyrite region
        dvtmp = np.interp(stsmcsmooth.P[iii],parkp,parkdv)
        stsmcsmooth.rho_s[iii]=1/( (1/stsmcsmooth.rho[iii]) - dvtmp )
        stsmcsmooth.S_s[iii]=0
    if ((stsmcsmooth.P[iii]>= 240)):
        stsmcsmooth.rho_s[iii]=stsmcsmooth.rho[iii]*1.029 # last value from park et al. 2025
        stsmcsmooth.S_s[iii]=0    


# calculate low_high values; not implemented for pyrite phase
for iii in range(0,len(stsmcsmooth.P)):
    # have P and T, calculate isotherm
    sftpvt.isotherm_sft(stsmcsmooth.T[iii]+stsmcsmooth.T_std[iii])
    rhosftsmooth[iii] = np.interp(stsmcsmooth.P[iii],sftpvt.isotherm_p,sftpvt.rhoarr)
    if (stsmcsmooth.P[iii]> Pt3) & (stsmcsmooth.P[iii]<= Pt4) & (rhosftsmooth[iii] >= stsmcsmooth.rho[iii]*default_drho):
        stsmcsmooth.rho_s_low[iii]=rhosftsmooth[iii] 
for iii in range(0,len(stsmcsmooth.P)):
    # have P and T, calculate isotherm
    sftpvt.isotherm_sft(stsmcsmooth.T[iii]-stsmcsmooth.T_std[iii])
    rhosftsmooth[iii] = np.interp(stsmcsmooth.P[iii],sftpvt.isotherm_p,sftpvt.rhoarr)
    if (stsmcsmooth.P[iii]> Pt3) & (stsmcsmooth.P[iii]<= Pt4) & (rhosftsmooth[iii] >= stsmcsmooth.rho[iii]*default_drho):
        stsmcsmooth.rho_s_high[iii]=rhosftsmooth[iii] 



# check where entropy of the solid is zero and match the latest non-zero value
indzero = np.where(stsmcsmooth.S_s == 0.)[0]
indtie = np.flip(np.where(stsmcsmooth.S_s > 0.)[0])[0]
dShighP = stsmcsmooth.S_s[indtie]-stsmcsmooth.S_l[indtie]
stsmcsmooth.S_s[indzero] = stsmcsmooth.S_l[indzero] + dShighP+25 # pyrite a positive entropy shift Park et al. 2025

#stvatsft = np.interp(Pt3,stsmcsmooth.P,stsmcsmooth.S_s)
#print('stv at SFT REFPT = ',)
#dS_stvsft =  Ss_23_sft+sft_dS - stvatsft
#print('dS = ',dS_stvsft)

#######################################################################################
# Calculate densities along the CaCl2 to seifertite phase boundary
# these set in define_eos.py
#Pt3 = 132.   # GPa: CaCl2-type to seifertite on the melt curve
#Tt3 = 6100  # K
# solid-solid tie points
# CaCl2 to seifertite
#Pt23 = 120.2 # GPa
#Tt23 = 4000 # K

res23 = np.polyfit([Pt23,Pt3],[Tt23,Tt3],1)
if (Pt3 > Pt23):
    p23arr = np.arange(Pt23,Pt3,0.01)
else:
    p23arr = np.arange(Pt23,Pt3,-0.01)
t23arr = res23[1]+res23[0]*p23arr
#print(p23arr,t23arr)

p23bnd = hp.HPeos_PhaseBoundary()
p23bnd.P = p23arr
p23bnd.T = t23arr
p23bnd.rho_a = np.zeros(len(p23bnd.P))
p23bnd.rho_b = np.zeros(len(p23bnd.P))
p23bnd.S_a = np.zeros(len(p23bnd.P))
p23bnd.S_b = np.zeros(len(p23bnd.P))

# CACL2 IS PHASE A
# use the Driver model to calculate the boundary densities
for iii in range(0,len(p23bnd.P)):
    # have P and T, calculate isotherm
    stvpvt.isotherm_stv(p23bnd.T[iii])
    p23bnd.rho_a[iii] = np.interp(p23bnd.P[iii],stvpvt.isotherm_p,stvpvt.rhoarr)
    p23bnd.S_a[iii] = np.interp(p23bnd.P[iii],stvpvt.isotherm_p,stvpvt.isotherm_s)

# CACL2 IS PHASE B
# interpolate the stv-cacl2 isotherms for density on the boundary
for iii in range(0,len(p23bnd.P)):
    # have P and T, calculate isotherm
    sftpvt.isotherm_sft(p23bnd.T[iii])
    p23bnd.rho_b[iii]=np.interp(p23bnd.P[iii],sftpvt.isotherm_p,sftpvt.rhoarr)
    p23bnd.S_b[iii]=np.interp(p23bnd.P[iii],sftpvt.isotherm_p,sftpvt.isotherm_s)

#######################################################################################
# Calculate densities along the stv to CaCl2 phase boundary -- ERROR WIDTH in rho_a and rho_b
# these set in define_eos.py
# Pt2, Tt2, Pt12, Tt12
#
if 1:
    res12 = np.polyfit([Pt12,Pt2],[Tt12,Tt2],1)
    if (Pt2 > Pt12):
        p12arr = np.arange(Pt12,Pt2,0.01)
    else:
        p12arr = np.arange(Pt2,Pt12,0.01)
    t12arr = res12[1]+res12[0]*p12arr
    p12arr = np.concatenate(([63.1883],p12arr))
    t12arr = np.concatenate(([3000],t12arr))
if 0:
    res12 = np.polyfit([Pt2-1,Pt2],[4000,Tt2],1)
    p12arr = np.arange(Pt2-1,Pt2,0.01)
    t12arr = res12[1]+res12[0]*p12arr    

p12bnd = hp.HPeos_PhaseBoundary()
p12bnd.P = p12arr
p12bnd.T = t12arr
p12bnd.rho_a = np.zeros(len(p12bnd.P))
p12bnd.rho_b = np.zeros(len(p12bnd.P))
p12bnd.S_a = np.zeros(len(p12bnd.P))
p12bnd.S_b = np.zeros(len(p12bnd.P))

# STV IS PHASE A
# use the Driver model to calculate the boundary densities
for iii in range(0,len(p12bnd.P)):
    # have P and T, calculate isotherm
    stvpvt.isotherm_stv(p12bnd.T[iii])
    p12bnd.rho_a[iii] = np.interp(p12bnd.P[iii],stvpvt.isotherm_p,stvpvt.rhoarr)
    p12bnd.S_a[iii] = np.interp(p12bnd.P[iii],stvpvt.isotherm_p,stvpvt.isotherm_s)

# STV + 3 GPa IS PHASE B
dP = 0. # GPa
# interpolate the stv-cacl2 isotherms for density on the boundary
for iii in range(0,len(p12bnd.P)):
    # have P and T, calculate isotherm
    sftpvt.isotherm_sft(p12bnd.T[iii])
    p12bnd.rho_b[iii]=np.interp(p12bnd.P[iii]+dP,sftpvt.isotherm_p,sftpvt.rhoarr)
    p12bnd.S_b[iii]=np.interp(p12bnd.P[iii]+dP,sftpvt.isotherm_p,sftpvt.isotherm_s)

#######################################################################################
# Use Park et al. 2025 for Pyrite phase boundary
parkp = np.asarray([140,160,180,200,220,240])
parkdv = np.asarray([.319,.276,.243,.215,.191,.172])*1.6605391/sio2mol
ParkPt4 = 165 # GPa

