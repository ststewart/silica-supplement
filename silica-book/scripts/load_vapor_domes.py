# load_vapor_domes.py

## read in vapor dome

# Kraus 2012 vapor dome
dat_file   = silicadatadir+'Kraus-silica-vapor-dome-2.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.rho_l      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0])
dat.rho_v      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1])
dat.ent_v      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2])
dat.ent_l      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[3])
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[5])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[4])
dat.T_units      = ['K']
#print(dat.P)
#print(dat.ent_l)
#print(dat.T)
kdome = deepcopy(dat)

## read in Connolly data
C16highs = np.loadtxt(silicadatadir+'Connolly-high-cp.csv',delimiter=',',usecols=[0])
C16hight = np.loadtxt(silicadatadir+'Connolly-high-cp.csv',delimiter=',',usecols=[1])
C16bsrho  = np.loadtxt(silicadatadir+'Connolly-Fig4b-black-T-rho.csv',delimiter=',',usecols=[0])/1.e3
C16bstemp = np.loadtxt(silicadatadir+'Connolly-Fig4b-black-T-rho.csv',delimiter=',',usecols=[1])


C16lows = np.loadtxt(silicadatadir+'Connolly-low-cp.csv',delimiter=',',usecols=[0])
C16lowt = np.loadtxt(silicadatadir+'Connolly-low-cp.csv',delimiter=',',usecols=[1])
C16bdrho  = np.loadtxt(silicadatadir+'Connolly-Fig4b-black-dashed-T-rho.csv',delimiter=',',usecols=[0])/1.e3
C16bdtemp = np.loadtxt(silicadatadir+'Connolly-Fig4b-black-dashed-T-rho.csv',delimiter=',',usecols=[1])

C16sils = np.loadtxt(silicadatadir+'Connolly-SiL-VWRK.csv',delimiter=',',usecols=[0])
C16silt = np.loadtxt(silicadatadir+'Connolly-SiL-VWRK.csv',delimiter=',',usecols=[1])
C16rsrho  = np.loadtxt(silicadatadir+'Connolly-Fig4b-red-T-rho.csv',delimiter=',',usecols=[0])/1.e3
C16rstemp = np.loadtxt(silicadatadir+'Connolly-Fig4b-red-T-rho.csv',delimiter=',',usecols=[1])
