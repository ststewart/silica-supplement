# define-mc.py
print('################################################################')
print(' Read in fit for melt curve envelope ')

#book_path='./'
debug=False

# STS estimated melt curve wedge for data review
dat_file   = book_path+'scripts/melt-curve-wedge.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T_s      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])
dat.T_s_units      = ['K']
dat.T_l      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[2])
dat.T_l_units = ['K']
if debug:
    print(dat.P)
    print(dat.T_s)
    print(dat.T_l)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('America/Los_Angeles')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['judgement']
dat.groups = ['solid compounds and minerals']

stsmc = deepcopy(dat)



# STS estimated melt curve line for data review draw as dotted above 155 GPa
dat_file   = book_path+'scripts/melt-curve.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])
dat.T_units      = ['K']
dat.rho = np.zeros(len(dat.P)) # will be filled later
if debug:
    print(dat.P)
    print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('America/Los_Angeles')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['judgement']
dat.groups = ['solid compounds and minerals']

if 0:
    # interpolate to 0.5 GPa from 14.5 to 251 GPa
    mc_p = np.arange(min(dat.P),701.,0.1) # GPa
    mc_t = np.interp(mc_p,dat.P,dat.T)
    dat.P = mc_p
    dat.T = mc_t
dat.rho =  np.zeros(len(dat.P)) # will be filled later
dat.T_std = np.zeros(len(dat.P))+250 #K estimated 1sigma for the fiducial melt curve
stsmcsmooth = deepcopy(dat)
