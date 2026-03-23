# Karki et al. 2007 Liquid Silica EOS Model with MODIFICATIONS
# not robust above 5 g/cm3 needs development

import intersection as itmp

karkireftemporiginal=3000

#--------------------------------------
# Load Karki isotherms 

dat_file   = silicadatadir+'Karki-2007-isotherms-v2.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
vol=np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.rho      = sio2mol/vol*1.6605391 # amu/a^3 -> g/cm3
dat.rho_units = ['g/cm3']
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[2])
dat.T_units      = ['K']
#print(dat.P)
#print(dat.rho)
#print(dat.T)
karki = deepcopy(dat)
dat=[]

#--------------------------------------
# Load Karki gammas 2007
v0 = 45.8 # A/SiO2
rho0 = sio2mol/v0*1.6605391 # g/cm3

# STS estimated melt curve line for data review draw as dotted above 155 GPa
dat_file   = silicadatadir+'karki-gamma.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
vol=np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])*v0
dat.rho      = sio2mol/vol*1.6605391 # amu/a^3 -> g/cm3
dat.rho_units = ['g/cm3']
dat.gamma=np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])

dat_file   = silicadatadir+'karki-cv.txt'
dat.cv=np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])*3*8.314/(sio2mol/1000) # J/K/kg

#print(dat.gamma)
#print(dat.rho)
#print(dat.cv)
karkigcvorig = deepcopy(dat)
dat=[]
vol=[]


#--------------------------------------
## Code a set of Karki original PVT isotherms 3000-7000 K

# Karki 3rd order polynomial fit
# 3000 Isotherm reference!
v0 = 45.8 # A/SiO2
rho0 = sio2mol/v0*1.6605391 # g/cm3
K0 = 5.2 # GPa
K1 = 22.5 # dimless
K1 = 7.4 # dimless
K0K2 = -452.3 # dimless
K0K2=-15
rhoarr = np.arange(.5,3.,.01)*rho0 # 2.8 to 6.5 g/cm3 covers range for phase diagram plots
vaarr = sio2mol*1.6605391/rhoarr

# 3000 K isotherm
ind = np.where(karki.T == 3000)[0]
# add a low density pressure point for the polynomial fit
xfit = np.concatenate((karki.rho[ind],np.asarray([.8])))
yfit = np.concatenate((karki.P[ind],np.asarray([0.])))    
res = np.polyfit(xfit,yfit,5)
parr3 = res[5]+res[4]*rhoarr+res[3]*rhoarr*rhoarr+res[2]*np.power(rhoarr,3)+res[1]*np.power(rhoarr,4)+res[0]*np.power(rhoarr,5)
print('first fit parr3 = ',np.interp(5.8,rhoarr,parr3))
if 0:
    # don't tweak for now
    xfit = np.concatenate((karki.rho[ind],np.asarray([.8,5.8])))
    yfit = np.concatenate((karki.P[ind],np.asarray([0.,210])))    
    res = np.polyfit(xfit,yfit,5)
    parr3 = res[5]+res[4]*rhoarr+res[3]*rhoarr*rhoarr+res[2]*np.power(rhoarr,3)+res[1]*np.power(rhoarr,4)+res[0]*np.power(rhoarr,5)
    print('second fit parr3 = ',np.interp(5.8,rhoarr,parr3))

T0=karkireftemporiginal # K Original reference temperature for Karki thermal pressure fit
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm6000 =  BV*(6000-T0)
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm5000 =  BV*(5000-T0)
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm4000 =  BV*(4000-T0)
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm7000 =  BV*(7000-T0)

BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm3130 =  BV*(3130-T0)


BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm3500 =  BV*(3500-T0)
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm4500 =  BV*(4500-T0)
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm5500 =  BV*(5500-T0)

BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm5900 =  BV*(5900-T0)

BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm6500 =  BV*(6500-T0)
BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
ptherm8000 =  BV*(8000-T0)

karkiorig = hp.HPeos_Model()
karkiorig.stempsarr = np.asarray([3000,3130,3500,4000,4500,5000,5500,5900,6000,6500,7000,8000]) # pinned by at refpt and cv
        
karkiorig.rhoref = rho0
karkiorig.rhoarr = deepcopy(rhoarr)
karkiorig.vaarr = deepcopy(vaarr)
karkiorig.p3000 = deepcopy(parr3)
karkiorig.t3000 = np.zeros(len(rhoarr))+3000.
karkiorig.p3130 = parr3+ptherm3130
karkiorig.t3130 = np.zeros(len(rhoarr))+3130.
karkiorig.p3500 = parr3+ptherm3500
karkiorig.t3500 = np.zeros(len(rhoarr))+3500.
karkiorig.p4000 = parr3+ptherm4000
karkiorig.t4000 = np.zeros(len(rhoarr))+4000.
karkiorig.p4500 = parr3+ptherm4500
karkiorig.t4500 = np.zeros(len(rhoarr))+4500.
karkiorig.p5000 = parr3+ptherm5000
karkiorig.t5000 = np.zeros(len(rhoarr))+5000.
karkiorig.p5500 = parr3+ptherm5500
karkiorig.t5500 = np.zeros(len(rhoarr))+5500.
karkiorig.p5900 = parr3+ptherm5900
karkiorig.t5900 = np.zeros(len(rhoarr))+5900.
karkiorig.p6000 = parr3+ptherm6000
karkiorig.t6000 = np.zeros(len(rhoarr))+6000.
karkiorig.p6500 = parr3+ptherm6500
karkiorig.t6500 = np.zeros(len(rhoarr))+6500.
karkiorig.p7000 = parr3+ptherm7000
karkiorig.t7000 = np.zeros(len(rhoarr))+7000.
karkiorig.p8000 = parr3+ptherm8000
karkiorig.t8000 = np.zeros(len(rhoarr))+8000.

parr3=[]
rhoarr=[]
vaarr=[]
BV=[]
rhoarr3=[]


#--------------------------------------
# fit the gamma

x = karkigcvorig.rho
y = karkigcvorig.gamma
#x = [2,2.5,3.,3.7,4.62,4.87,5.0]
#y = [.2,.2,.2,.5,.83,.95,1.03]

if 0:
    res = np.polyfit(x,y,4)
    karkiorig.garr = res[4]+res[3]*karkiorig.rhoarr+res[2]*np.power(karkiorig.rhoarr,2)+res[1]*np.power(karkiorig.rhoarr,3)+res[0]*np.power(karkiorig.rhoarr,4)
    #ind = np.where(karkiorig.rhoarr > 5.65)[0]
    #karkiorig.garr[ind]=0.66
    #ind = np.where(karkiorig.rhoarr > 4.5)[0]
    #karkiorig.garr[ind]=1.1
    ind = np.where(karkiorig.rhoarr < 3.0)[0]
    karkiorig.garr[ind]=0.18
if 1:    
    rho0=3.0
    varr = 1./karkiorig.rhoarr
    xxx_g0 = 0.17
    #xxx_q  = -2.2
    xxx_q  = -3.5
    xxx_garr = xxx_g0*np.power((varr*rho0),xxx_q)
    karkiorig.garr = xxx_garr
    ind = np.where(karkiorig.rhoarr < 3.)[0]
    karkiorig.garr[ind]=0.17
    #ind = np.where(karkiorig.rhoarr > 4.4)[0]
    #karkiorig.garr[ind]=1.1

if 0:
    rho0=3.5
    varr = 1./karkiorig.rhoarr
    xxx_g0 = .2
    #xxx_q  = -2.2
    xxx_q  = -4
    xxx_garr = xxx_g0*np.power((varr*rho0),xxx_q)
    karkiorig.garr = xxx_garr
    ind = np.where(karkiorig.rhoarr < rho0)[0]
    karkiorig.garr[ind]=0.17

# fit the cv at high density
cvavg = np.average(karkigcvorig.cv[3::])
karkiorig.cv = cvavg # J/K/kg

#--------------------------------------
## Apply constant density shift to the original Karki model

karkipvt = deepcopy(karkiorig)
# liquid_cvfactor set in define_eos.py
karkipvt.cv = liquid_cvfactor*3*Rgas/sio2mol*1000 # J/K/kg
karkigcv = deepcopy(karkigcvorig) 

karkipvt.rhoarr = deepcopy(karkiorig.rhoarr)+lda_liquid_density_shift
karkipvt.varr = 1./karkipvt.rhoarr
karkipvt.vaarr = sio2mol*1.6605391/karkiorig.rhoarr ## ORIG OTHERWISE WRONG PTH

if 0:
    # adjust the gamma model to match seifertite over the peak
    v0dft= 154.8 # bohr^3/sio2
    v0 = v0dft/sio2mol*np.power(0.529177249,3)/1.6605391/1000.
    sft_r0 = 1/v0/1.e3+sft_density_shift # g/cm3
    sft_g0 = 1.3
    sft_q  = 2.05
    sft_garr = sft_g0*np.power((karkipvt.varr*sft_r0),sft_q)+.1
    ind = np.where(karkipvt.rhoarr > 5.)[0]
    karkipvt.garr[ind]=sft_garr[ind]
    
    #test entropy
    #ind = np.where(karkipvt.rhoarr > 5.)[0]
    #karkipvt.garr[ind]=.1
if 0:
    # adjust the gamma model to MATCH QTZ LIQUID BRANCH CROSSING LIQUIDUS
    # this is a good match to the OCampo Results
    v0dft= 154.8 # bohr^3/sio2
    v0 = v0dft/sio2mol*np.power(0.529177249,3)/1.6605391/1000.
    sft_r0 = 1/v0/1.e3+sft_density_shift # g/cm3
    sft_g0 = 3.
    sft_q  = 10
    sft_garr = sft_g0*np.power((karkipvt.varr*sft_r0),sft_q)+.1
    ind = np.where(karkipvt.rhoarr > 4.9)[0]
    karkipvt.garr[ind]=sft_garr[ind]

if 1:
    # Using  OCampo et al. 2025 - very similar to empirical model found above
    A00=0.381
    B00=5.033
    C00=0.1056
    
    v00=0.4627 # cm3/g
    gammaO25low = 0.14*np.pow(karkipvt.varr/v00,-2.39)
    ind = np.where(karkipvt.rhoarr <= B00)[0]
    karkipvt.garr[ind]=gammaO25low[ind]

    gammaO25high = 2/3 + A00*np.exp(-np.pow(np.log(karkipvt.rhoarr/B00)/C00,2))
    ind = np.where(karkipvt.rhoarr>B00)[0]
    karkipvt.garr[ind]=gammaO25high[ind]
    #ind = np.where(karkipvt.rhoarr < 2.)[0]
    #karkipvt.garr[ind]=0.2
#------------------------------------------------
# calculate isentropes from the liquid_reference_point using a constant heat capacity, cvavg

# (rho,p,t,s) in (g/cm3,GPa,K,J/K/kg) # match Connolly vapor curve at 1 bar 
liquid_reference_point[0] = np.interp(liquid_reference_point[1],karkipvt.p4500,karkipvt.rhoarr)
karkipvt.refpt = liquid_reference_point

for iii in range(len(karkipvt.stempsarr)):
    # set conditions for this isentrope
    rhoref = karkipvt.refpt[0]
    rhoarr = deepcopy(karkipvt.rhoarr)
    garr = deepcopy(karkipvt.garr)
    tarr = np.zeros(len(rhoarr))
    parr = np.zeros(len(rhoarr))
    t0 = karkipvt.stempsarr[iii] # temp at rhoref

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
    # calculate the pressure on each isentrope
    # T, calculate the PV isotherm, find the P that matches the rho
    BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-karkipvt.vaarr)/1.86)))) # GPa
    parr = karkipvt.p3000+BV*(tarr-karkireftemporiginal) # karki thermal model reference temp 3000 K
    
    if t0 == 3000:
        karkipvt.st3000=deepcopy(tarr)
        karkipvt.sp3000=deepcopy(parr)
    if t0 == 3130:
        karkipvt.st3130=deepcopy(tarr)
        karkipvt.sp3130=deepcopy(parr)
    if t0 == 3500:
        karkipvt.st3500=deepcopy(tarr)
        karkipvt.sp3500=deepcopy(parr)
    if t0 == 4000:
        karkipvt.st4000=deepcopy(tarr)
        karkipvt.sp4000=deepcopy(parr)
    if t0 == 4500:
        karkipvt.st4500=deepcopy(tarr)
        karkipvt.sp4500=deepcopy(parr)
    if t0 == 5000:
        karkipvt.st5000=deepcopy(tarr)
        karkipvt.sp5000=deepcopy(parr)
    if t0 == 5500:
        karkipvt.st5500=deepcopy(tarr)
        karkipvt.sp5500=deepcopy(parr)
    if t0 == 5900:
        karkipvt.st5900=deepcopy(tarr)
        karkipvt.sp5900=deepcopy(parr)
    if t0 == 6000:
        karkipvt.st6000=deepcopy(tarr)
        karkipvt.sp6000=deepcopy(parr)
    if t0 == 6500:
        karkipvt.st6500=deepcopy(tarr)
        karkipvt.sp6500=deepcopy(parr)
    if t0 == 7000:
        karkipvt.st7000=deepcopy(tarr)
        karkipvt.sp7000=deepcopy(parr)
    if t0 == 8000:
        karkipvt.st8000=deepcopy(tarr)
        karkipvt.sp8000=deepcopy(parr)
    tarr=[]
    parr=[]

karkipvt.svalsarr = np.zeros(len(karkipvt.stempsarr))
for iii in range(len(karkipvt.stempsarr)):
    karkipvt.svalsarr[iii] = (np.log(karkipvt.stempsarr[iii])-np.log(karkipvt.refpt[2]))*karkipvt.cv+karkipvt.refpt[3]
    #karkipvt.svalsarr[iii] = (np.log(karkipvt.stempsarr[iii])-np.log(4500.))*karkipvt.cv+karkipvt.refpt[3]
    #print(iii,karkipvt.stempsarr[iii],karkipvt.svalsarr[iii])


#------------------------------------------------
## calculate melt curve V and S from Karki PVT

# use the Karki model to calculate the liquidus densities for SiO2 based on stsmcsmooth array
# Karki Model starts to diverge from DFT near density of 5.0 g/cm3
# replace with fit to DFT liquidus above this density
# the entropies are also incorrect; will need to be fixed later; for now, setting to zero

stsmcsmooth.S_l = np.zeros(len(stsmcsmooth.P))
stsmcsmooth.rho = np.zeros(len(stsmcsmooth.P))
for iii in range(0,len(stsmcsmooth.P)):
    #print(stsmcsmooth.T[iii])
    # calculate the isotherm first
    BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-karkipvt.vaarr)/1.86)))) # GPa
    pthermt =  BV*(stsmcsmooth.T[iii]-karkireftemporiginal) # Karki model reference temp
    parrt = karkipvt.p3000+pthermt
    # now we have P-V-T at this temperature; interpolate for density
    stsmcsmooth.rho[iii] = np.interp(stsmcsmooth.P[iii],parrt,karkipvt.rhoarr)
    # now that we have rho; find the entropy that goes with this point
    rrr = np.where(karkipvt.rhoarr >= stsmcsmooth.rho[iii])[0][0]
    # interpolate entsarr vs. T
    tempfit = np.asarray([karkipvt.st3000[rrr],karkipvt.st3130[rrr],\
                          karkipvt.st3500[rrr],karkipvt.st4000[rrr],\
                          karkipvt.st4500[rrr],karkipvt.st5000[rrr],\
                          karkipvt.st5500[rrr],karkipvt.st5900[rrr],karkipvt.st6000[rrr],\
                          karkipvt.st6500[rrr],karkipvt.st7000[rrr],karkipvt.st8000[rrr]])
    stsmcsmooth.S_l[iii] = np.interp(stsmcsmooth.T[iii],tempfit,karkipvt.svalsarr)
    # interpolate reference energy vs. T
    tempfit = np.asarray([karkipvt.st3000[rrr],karkipvt.st3130[rrr],\
                          karkipvt.st3500[rrr],karkipvt.st4000[rrr],\
                          karkipvt.st4500[rrr],karkipvt.st5000[rrr],\
                          karkipvt.st5500[rrr],karkipvt.st5900[rrr],karkipvt.st6000[rrr],\
                          karkipvt.st6500[rrr],karkipvt.st7000[rrr],karkipvt.st8000[rrr]])

# calculate uncertainty on rho associated with T_std
stsmcsmooth.rho_low = np.zeros(len(stsmcsmooth.rho))
stsmcsmooth.rho_high = np.zeros(len(stsmcsmooth.rho))

for iii in range(0,len(stsmcsmooth.P)):
    #print(stsmcsmooth.T[iii])
    # calculate the isotherm first
    BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-karkipvt.vaarr)/1.86)))) # GPa
    pthermt =  BV*(stsmcsmooth.T[iii]-stsmcsmooth.T_std[iii]-karkireftemporiginal) # Karki model reference temp
    parrt = karkipvt.p3000+pthermt
    # now we have P-V-T at this temperature; interpolate for density
    stsmcsmooth.rho_high[iii] = np.interp(stsmcsmooth.P[iii],parrt,karkipvt.rhoarr)
# calculate uncertainty on rho associated with T_std
for iii in range(0,len(stsmcsmooth.P)):
    #print(stsmcsmooth.T[iii])
    # calculate the isotherm first
    BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-karkipvt.vaarr)/1.86)))) # GPa
    pthermt =  BV*(stsmcsmooth.T[iii]+stsmcsmooth.T_std[iii]-karkireftemporiginal) # Karki model reference temp
    parrt = karkipvt.p3000+pthermt
    # now we have P-V-T at this temperature; interpolate for density
    stsmcsmooth.rho_low[iii] = np.interp(stsmcsmooth.P[iii],parrt,karkipvt.rhoarr)


## Load the DFT liquidus fit -- calcualted in melt-curve.ipynb
rho_l_fit = hp.HPeos_Fit()
rho_l_fit.outputfitpath = modelfitpath # in load_paths.py
rho_l_fit.filename = 'rho_l_fit_quadratic'
rho_l_fit.loadjson()
mc_new = hp.HPeos_MeltCurve()
mc_new.rho_l = np.arange(3.5,18,.01) # valid range for fit
mc_new.P_l = rho_l_fit.a[0]+rho_l_fit.a[1]*mc_new.rho_l+rho_l_fit.a[2]*mc_new.rho_l*mc_new.rho_l

# replace liquidus at intersection point
x,y=itmp.intersection(mc_new.rho_l,mc_new.P_l,stsmcsmooth.rho,stsmcsmooth.P)
print('LIQUIDUS TRANSITION RHO =',x,y)

if len(y)>1:
    # second intersection point at 5.02 g/cm3 and 121 GPa
    ind = np.where(stsmcsmooth.P > y[1])
    stsmcsmooth.rho[ind]=np.interp(stsmcsmooth.P[ind],mc_new.P_l,mc_new.rho_l)
    stsmcsmooth.rho_low[ind]=stsmcsmooth.rho[ind]*.995
    stsmcsmooth.rho_high[ind]=stsmcsmooth.rho[ind]*1.005
    
    # unsure about what to do with S when switching to high-pressure melt curve
    # this is near the transition from pyrite to cotunnite structure on the melt curve
    # for now fix a constant entropy until we learn more about the heat capacity in this phase space
    stsmcsmooth.S_l[ind]=stsmcsmooth.S_l[ind[0]]
