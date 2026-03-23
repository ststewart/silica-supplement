# load_silica_data.py
# Updated 5/8/2025 STSM
#
# set path to data directories in load_paths.py
#rusbank_path = '../py-rusbank/'
#silicadatadir = '../data_silica_db/'
debug=False


################################################################
################################################################
#########       FUNCTIONS FOR DATABASE             #############
################################################################
################################################################

def bibkeyid(rlist, bibkey):
    for ii in range(len(rlist)):
    #print(ii,rlist[ii].bibkey)
        if rlist[ii].bibkey == [bibkey]:
            return ii


################################################################
################################################################
#########       LOAD RUSBANK DATASETS FIRST        #############
################################################################
################################################################
#
# If needed, open the file in binary mode
if True:
#if False:
    with open(rusbank_path+'rusbankhug_list.pkl', 'rb') as f:      
        # Call load method to deserialze
        rusbankhug_list = pickle.load(f)
    with open(rusbank_path+'rusbankref_list.pkl', 'rb') as f:      
        # Call load method to deserialze
        ref_list = pickle.load(f)

for item in ref_list:
    item.bibkey = []

## gather Rusbank quartz data
for item in rusbankhug_list:
    if item.matid == [174]: # Quartz = alpha-quartz
        rbqtz = deepcopy(item)
    if item.matid == [143]: # Quartz amorphous = fused silica|
        rbfs = deepcopy(item)
indqtz = np.where(np.asarray(rbqtz.m) == 1)[0]
qtzdf = hp.rusbanktodf([rbqtz])[0]
indfs = np.where(np.asarray(rbfs.m) == 1)[0]
fsdf = hp.rusbanktodf([rbfs])[0]

################################################################
################################################################
####            LOAD EOS STATE DATA                 ############
################################################################
################################################################
if 0: # move to hugoniot-calcs.py
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

#print'################################################################')
## Experimental Shock Hugoniot: Lyzenga and Ahrens GRL 1980
# The database lists alpha-quartz and amorphous quartz (fused silica) as different materials.
# So load into different database entries.

# add reference to database reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Lyzenga, G.A. and Ahrens, T.J., 1980. Shock temperature measurements in Mg2SiO4 and SiO2 at high pressures. Geophysical Research Letters, 7(2), pp.141-144.']
newref.doi = ['10.1029/GL007i002p00141']
bibstr = "@article{Lyzenga1980,\
  doi = {10.1029/gl007i002p00141},\
  url = {https://doi.org/10.1029/gl007i002p00141},\
  year = {1980},\
  month = feb,\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {7},\
  number = {2},\
  pages = {141--144},\
  author = {Gregory A. Lyzenga and Thomas J. Ahrens},\
  title = {Shock temperature measurements in {Mg$_2$SiO$_4$} and {SiO$_2$} at high pressures},\
  journal = {Geophysical Research Letters}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Lyzenga1980']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
#print(newref.reference[0])

# Lyzenga, G.A., Ahrens, T.J., 1980. Shock temperature measurements in Mg 2 SiO 4 and SiO 2 at high pressures. Geophys. Res. Lett. 7, 141–144. https://doi.org/10.1029/GL007i002p00141
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Lyzenga1980-qtz-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.650] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=7,usecols=[2],max_rows=7) 
hug.P_std  = np.loadtxt(hug_file,delimiter=',',skiprows=7,usecols=[3],max_rows=7) 
hug.P_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=7,usecols=[4],max_rows=7) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=7,usecols=[5],max_rows=7)
hug.T_units = ['K']
hug.m      = [1.0]*len(hug.P)
#rho = Q_(EOS.rho,config['tableunits']['density'])
#print(hug.P)
#print(hug.T)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Ta flyer, driver, impedance match Wackerle 1963 Hugoniot']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

# first new entry
shockdat_list = [hug]


# Lyzenga, G.A., Ahrens, T.J., 1980. Shock temperature measurements in Mg 2 SiO 4 and SiO 2 at high pressures. Geophys. Res. Lett. 7, 141–144. https://doi.org/10.1029/GL007i002p00141
# FUSED QUARTZ
hug_file     = silicadatadir+'Lyzenga1980-qtz-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[2],max_rows=7) 
hug.P_std  = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[3],max_rows=7) 
hug.P_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[4],max_rows=7) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[5],max_rows=7) 
hug.T_units = ['K']
hug.m      = [1.0]*len(hug.P)
#print(hug.P)
#print(hug.T)


hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Ta flyer, driver, impedance match Wackerle 1963 Hugoniot']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


################################################################
#print('################################################################')
## Lyzenga et al. 1983
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Lyzenga, G.A., Ahrens, T.J. and Mitchell, A.C., 1983. Shock temperatures of SiO2 and their geophysical implications. Journal of Geophysical Research: Solid Earth, 88(B3), pp.2431-2444.']
newref.doi = ['10.1029/JB088iB03p02431']
bibstr = "@article{Lyzenga1983,\
  doi = {10.1029/jb088ib03p02431},\
  url = {https://doi.org/10.1029/jb088ib03p02431},\
  year = {1983},\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {88},\
  number = {B3},\
  pages = {2431},\
  author = {Gregory A. Lyzenga and Thomas J. Ahrens and Arthur C. Mitchell},\
  title = {Shock temperatures of SiO₂ and their geophysical implications},\
  journal = {Journal of Geophysical Research}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Lyzenga1983']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

# Lyzenga, G.A., Ahrens, T.J., Mitchell, JGR 1983.
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Lyzenga1983-qtz-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.650] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[2],max_rows=9) 
hug.P_std  = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[3],max_rows=9) 
hug.P_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[4],max_rows=9) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[5],max_rows=9) 
hug.T_units = ['K']
hug.m      = [1.0]*len(hug.P)
#rho = Q_(EOS.rho,config['tableunits']['density'])
print(hug.P)
print(hug.T)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Ta flyer, driver, impedance match new qtz hugoniot in paper']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

# Lyzenga, G.A., Ahrens, T.J., Mitchell 1983.
# FUSED QUARTZ
hug_file     = silicadatadir+'Lyzenga1983-qtz-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[2],max_rows=7) 
hug.P_std  = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[3],max_rows=7) 
hug.P_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[4],max_rows=7) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[5],max_rows=7) 
hug.T_units = ['K']
hug.m      = [1.0]*len(hug.P)
print(hug.P)
print(hug.T)


hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Ta flyer, driver, impedance match new qtz hugoniot in paper']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


## Lyzenga et al. 1983 Quartz Hugoniot states - CRYSTALLINE QUARTZ
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Lyzenga1983-qtz-Hug.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.650] # g/cm3
hug.uflyer  = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2],max_rows=6) 
hug.uflyer_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3],max_rows=6) 
hug.uflyers_units = ['km/s']
hug.us       = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4],max_rows=6) 
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[5],max_rows=6) 
hug.us_units = ['km/s']
hug.up       = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[6],max_rows=6) 
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[7],max_rows=6) 
hug.up_units = ['km/s']
hug.P       = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[8],max_rows=6) 
hug.P_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[9],max_rows=6) 
hug.P_units = ['GPa']
hug.rho       = 1./((1./hug.rho_ref[0])*(1.-hug.up/hug.us))
hug.rho_std   = np.sqrt(hug.rho*hug.rho*(hug.us_std*hug.us_std/hug.us/hug.us+hug.up_std*hug.up_std/hug.up/hug.up))
hug.rho_units = ['g/cm3']
hug.m       = [1.0]*len(hug.P)
print(hug.P)
print(hug.rho)


hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Ta flyer, driver, impedance match new qtz hugoniot in paper']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


## Lyzenga et al. 1983 Quartz Hugoniot states - CRYSTALLINE QUARTZ
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Lyzenga1983-qtz-Hug.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.uflyer  = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[2],max_rows=5) 
hug.uflyer_std   = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[3],max_rows=5) 
hug.uflyers_units = ['km/s']
hug.us       = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[4],max_rows=5) 
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[5],max_rows=5) 
hug.us_units = ['km/s']
hug.up       = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[6],max_rows=5) 
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[7],max_rows=5) 
hug.up_units = ['km/s']
hug.P       = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[8],max_rows=5) 
hug.P_std   = np.loadtxt(hug_file,delimiter=',',skiprows=9,usecols=[9],max_rows=5) 
hug.P_units = ['GPa']
hug.rho       = 1./((1./hug.rho_ref[0])*(1.-hug.up/hug.us))
hug.rho_std   = np.sqrt(hug.rho*hug.rho*(hug.us_std*hug.us_std/hug.us/hug.us+hug.up_std*hug.up_std/hug.up/hug.up))
hug.rho_units = ['g/cm3']
hug.m       = [1.0]*len(hug.P)
print(hug.P)
print(hug.rho)


hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Ta flyer, driver, impedance match new qtz hugoniot in paper']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

################################################################
print('################################################################')
## Boslough JGR 1988 Postshock temperatures in silica

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Boslough, M.B., 1988. Postshock temperatures in silica. Journal of Geophysical Research: Solid Earth, 93(B6), pp.6477-6484.']
newref.doi = ['10.1029/JB093iB06p06477']
bibstr = "@article{Boslough1988,\
  doi = {10.1029/jb093ib06p06477},\
  url = {https://doi.org/10.1029/jb093ib06p06477},\
  year = {1988},\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {93},\
  number = {B6},\
  pages = {6477},\
  author = {M. B. Boslough},\
  title = {Postshock temperatures in silica},\
  journal = {Journal of Geophysical Research}\
}}"
newref.bibtex = [bibstr]
newref.bibkey = ['Boslough1988']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Boslough 1988 CRYSTALLINE Quartz
hug_file     = silicadatadir+'Boslough-1988-silica-shock-temp-release-withsigp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.65] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=4,usecols=[0],max_rows=9) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=4,usecols=[7],max_rows=9) 
hug.P_std_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=4,usecols=[1],max_rows=9) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=4,usecols=[2],max_rows=9) 
hug.T_units = ['K']
hug.T_rel      = np.loadtxt(hug_file,delimiter=',',skiprows=4,usecols=[3],max_rows=9) 
hug.T_rel_std  = np.loadtxt(hug_file,delimiter=',',skiprows=4,usecols=[4],max_rows=9) 
hug.T_rel_units = ['K']
hug.m      = [1.0]*len(hug.P)
print(hug.P)
print(hug.T)
print(hug.T_rel)


hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Reanalysis of Lyzenga et al. 1983 with added post-shock temperatures.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


# Boslough 1988 AMORPHOUS QUARTZ
hug_file     = silicadatadir+'Boslough-1988-silica-shock-temp-release-withsigp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[0],max_rows=7) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[7],max_rows=7) 
hug.P_std_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[1],max_rows=7) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[2],max_rows=7) 
hug.T_units = ['K']
hug.T_rel      = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[3],max_rows=7) 
hug.T_rel_std  = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[4],max_rows=7) 
hug.T_rel_units = ['K']
hug.m      = [1.0]*len(hug.P)
print(hug.P)
print(hug.T)
print(hug.T_rel)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Reanalysis of Lyzenga et al. 1983 with added post-shock temperatures.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


################################################################
print('################################################################')
## Luo et al. 2004

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Luo, S.N., Akins, J.A., Ahrens, T.J. and Asimow, P.D., 2004. Shock‐compressed MgSiO3 glass, enstatite, olivine, and quartz: Optical emission, temperatures, and melting. Journal of Geophysical Research: Solid Earth, 109(B5).']
newref.doi = ['10.1029/2003JB002860']
bibstr = "@article{Luo2004,\
  doi = {10.1029/2003jb002860},\
  url = {https://doi.org/10.1029/2003jb002860},\
  year = {2004},\
  month = may,\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {109},\
  number = {B5},\
  author = {Sheng-Nian Luo and Joseph A. Akins and Thomas J. Ahrens and Paul D. Asimow},\
  title = {Shock-compressed MgSiO₃ glass,  enstatite,  olivine,  and quartz: Optical emission,  temperatures,  and melting},\
  journal = {Journal of Geophysical Research: Solid Earth}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Luo2004']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Luo et al. JGR 2004
# CRYSTALLINE QUARTZ - X CUT
hug_file   = silicadatadir+'Luo2004-qtz-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.phase_ref = ['single crystal']
hug.orientation = ['x']
hug.rho_ref = [2.650] # g/cm3
hug.P      = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2])]
hug.P_std  = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3])]
hug.P_units = ['GPa']
hug.T      = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4]) ]
hug.T_std  = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[5]) ]
hug.T_units = ['K']
hug.us      = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[8]) ]
hug.us_std  = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[9]) ]
hug.us_units = ['km/s']
hug.up      = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[10]) ]
hug.up_std  = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[11]) ]
hug.up_units = ['km/s']
hug.rho      = [np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[12])] 
hug.rho_std  =[ np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[13]) ]
hug.rho_units = ['g/cm3']
hug.orientation = ['x']
print(hug.P)
print(hug.T)
#rho = Q_(EOS.rho,config['tableunits']['density'])

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]
hug.remlist = ['Ta flyer, driver, impedance match']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)



################################################################
if 0: # NOT USED OLD READ SET
    print('################################################################')
    ## Amodeo 2023 Thesis
    ## Not using this older folder; updated with additional shots in 2024; see subfolder imp-match-UCD
    
    # add reference to reference list
    
    maxrefid=0
    for item in ref_list:
        if int(item.refid[0]) > maxrefid:
            maxrefid = int(item.refid[0])
            
    print('Previous max refid # = ',maxrefid)
    newrefid = maxrefid+1
    print('New refid # = ',newrefid)
    
    newref = hp.HPeos_Reference()
    newref.refid = [newrefid]
    newref.reference = ['Amodeo 2023 Thesis']
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
    
    ref_list.append(newref)
    #pprint(list(newref))
    print(newref.reference[0])
    
    if 0:
        shockdat_list.append(hug) # replace amodeo qtz
        shockdat_list.append(hug) # replace amodeo fs
    
    if 1:
        result = subprocess.run(['whoami'],capture_output=True)
        #help(result)
        if result.stdout.decode() == 'sts\n':
            print('yes')
    
        print(result.stdout.decode())
        # Amodeo thesis 2023
        # CRYSTALLINE QUARTZ
        shot_ind_q = ["S44", "S27", "S28", "S32", "S45", "S26", "S33"]
        if result.stdout.decode() == 'sstewart\n':
            BOXPATH = r"/Users/sstewart/Library/CloudStorage/Box-Box/Quartz Fused Silica Plots - Shared/"
        if result.stdout.decode() == 'sts\n':
            BOXPATH = r"/Users/sts/Library/CloudStorage/Box-Box/Quartz Fused Silica Plots - Shared/"
        if result.stdout.decode() == 'kaitlyn\n':
            BOXPATH = r"/Users/kaitlyn/Library/CloudStorage/Box-Box/Quartz Fused Silica Plots - Shared/"
    
        hug_file   = BOXPATH+'Master plotting/IM_solutions_SiO2.txt'
    
        IM = pd.read_csv(hug_file, sep=" ")
        hug          = hp.HPeos_Shock()
        hug.matid    = [174]
        hug.material = ['Quartz']
        hug.formula  = ['SiO2']
        hug.rho_ref = [2.650] # g/cm3
        for i in range(len(shot_ind_q)):
            ind = np.where(IM["ID"] == shot_ind_q[i])
            ind = ind[0]
            hug.exp_id.append(shot_ind_q[i])
            hug.P.append(IM.loc[IM.ID==shot_ind_q[i],'P'].values[0]/1.e9)
            hug.P_std.append(IM.loc[IM.ID==shot_ind_q[i],'P_er'].values[0]/1.e9)
            hug.P_units = ['GPa']
            hug.us.append(IM.loc[IM.ID==shot_ind_q[i],'Us'].values[0]/1.e3)
            hug.us_std.append(IM.loc[IM.ID==shot_ind_q[i],'Us_er'].values[0]/1.e3)
            hug.us_units = ['km/s']
            hug.up.append(IM.loc[IM.ID==shot_ind_q[i],'Up'].values[0]/1.e3)
            hug.up_std.append(IM.loc[IM.ID==shot_ind_q[i],'Up_er'].values[0]/1.e3)
            hug.up_units = ['km/s']
            hug.rho.append(IM.loc[IM.ID==shot_ind_q[i],'rho'].values[0]/1.e3)
            hug.rho_std.append(IM.loc[IM.ID==shot_ind_q[i],'rho_er'].values[0]/1.e3)
            hug.rho_units = ['g/cm3']
            if shot_ind_q[i] == "S26":
                hug.T.append(5273)
                hug.T_std.append(253)
                hug.T_rel.append(4091)
                hug.T_rel_std.append(132)
            elif shot_ind_q[i] == "S44":
                hug.T.append(3576)
                hug.T_std.append(57)
                hug.T_rel.append(0)
                hug.T_rel_std.append(0)
            elif shot_ind_q[i] == "S45":
                hug.T.append(5837)
                hug.T_std.append(213)
                hug.T_rel.append(3688)
                hug.T_rel_std.append(78)
            else:
                T_file = BOXPATH+shot_ind_q[i]+"/output.txt"
                Tout = pd.read_csv(T_file, sep=" ")
                hug.T.append(Tout["fit_T_MC"][0])
                hug.T_std.append(Tout["fit_T_MC_error"][0])
                hug.T_rel.append(Tout["fit_T_MC_R"][0]) 
                hug.T_rel_std.append(Tout["fit_T_MC_R_error"][0]) 
            hug.T_units = ['K']
            hug.T_rel_units = ['K']
            P = IM.loc[IM.ID==shot_ind_q[i],'P'].values[0]/1.e9
            hug.S.append(1793+22.4*P-0.0426*P*P+4.319e-5*P*P*P-1.742e-8*np.power(P,4))
        hug.m      = [1.0]*len(hug.P)
    
        #rho = Q_(EOS.rho,config['tableunits']['density'])
        print(hug.rem)
        print(hug.exp_id)
        print(hug.P)
        print(hug.T)
        print(hug.T_rel)
        print(hug.T_rel_std)
    
        hug.collections = ['Silica']
        hug.ref = [1]*len(hug.P)
        hug.reflist = ['',newref.reference]
        hug.refidlist = [0,newref.refid[0]]
        hug.rem = [1]*len(hug.P)
        hug.remlist = ['SS304 flyer, driver, impedance match qtz hugoniot in paper']
        DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
        hug.actions = [['Kaitlyn M. Amodeo','create',DateTime_in_ISOFormat,'new entry']]
        hug.diagnostics = ['streak camera','bandpass pyrometer', 'PDV']
        hug.platform = ['gas gun']
        hug.quality = [0]*len(hug.P)
        hug.groups = ['solid compounds and minerals']
    
        shockdat_list.append(hug)
    
        # S44, S27, S28 -- all have wiggles, 44 no IR
        # S26 bad impact vel laser; 33 is OK shock vel + Knudson
    
    
        # Amodeo 2023 thesis
        # FUSED QUARTZ
    
        shot_ind_fs = ["S41", "S25", "S31", "S23", "S20", "S22", "S21"]
        if result == 'sstewart':
            BOXPATH = r"/Users/sstewart/Library/CloudStorage/Box-Box/Quartz Fused Silica Plots - Shared/"
        if result == 'sts':
            BOXPATH = r"/Users/sts/Library/CloudStorage/Box-Box/Quartz Fused Silica Plots - Shared/"
    #    if result == 'sts':
    #        BOXPATH = "./data_silica_db/"
        if result.stdout.decode() == 'kaitlyn\n':
            BOXPATH = r"/Users/kaitlyn/Library/CloudStorage/Box-Box/Quartz Fused Silica Plots - Shared/"
    
        hug_file   = BOXPATH+'Master plotting/IM_solutions_SiO2.txt'
    
        IM = pd.read_csv(hug_file, sep=" ")
        hug          = hp.HPeos_Shock()
        hug.matid    = [143]
        hug.material = ['Quartz amorphous']
        hug.formula  = ['SiO2']
        hug.rho_ref = [2.204] # g/cm3
        for i in range(len(shot_ind_fs)):
            ind = np.where(IM["ID"] == shot_ind_fs[i])
            ind = ind[0]
            hug.exp_id.append(shot_ind_fs[i])
            hug.P.append(IM.loc[IM.ID==shot_ind_fs[i],'P'].values[0]/1.e9)
            hug.P_std.append(IM.loc[IM.ID==shot_ind_fs[i],'P_er'].values[0]/1.e9)
            hug.P_units = ['GPa']
            hug.us.append(IM.loc[IM.ID==shot_ind_fs[i],'Us'].values[0]/1.e3)
            hug.us_std.append(IM.loc[IM.ID==shot_ind_fs[i],'Us_er'].values[0]/1.e3)
            hug.us_units = ['km/s']
            hug.up.append(IM.loc[IM.ID==shot_ind_fs[i],'Up'].values[0]/1.e3)
            hug.up_std.append(IM.loc[IM.ID==shot_ind_fs[i],'Up_er'].values[0]/1.e3)
            hug.up_units = ['km/s']
            hug.rho.append(IM.loc[IM.ID==shot_ind_fs[i],'rho'].values[0]/1.e3)
            hug.rho_std.append(IM.loc[IM.ID==shot_ind_fs[i],'rho_er'].values[0]/1.e3)
            hug.rho_units = ['g/cm3']
            T_file = BOXPATH+shot_ind_fs[i]+"/output.txt"
            Tout = pd.read_csv(T_file, sep=" ")
            hug.T.append(Tout["fit_T_MC"][0])
            hug.T_std.append(Tout["fit_T_MC_error"][0])
            hug.T_rel.append(Tout["fit_T_MC_R"][0]) 
            hug.T_rel_std.append(Tout["fit_T_MC_R_error"][0]) 
            hug.T_units = ['K']
            hug.T_rel_units = ['K']
            hug.S.append(2294+26.82*IM.loc[IM.ID==shot_ind_fs[i],'P'].values[0]/1.e9-\
                   0.0723*IM.loc[IM.ID==shot_ind_fs[i],'P'].values[0]/1.e9*IM.loc[IM.ID==shot_ind_fs[i],'P'].values[0]/1.e9)
            
    
        hug.m      = [1.0]*len(hug.P)
        print(hug.rem)
        print(hug.exp_id)
        print(hug.P)
        print(hug.T)
        print(hug.T_rel)
        print(hug.T_rel_std)
    
        hug.collections = ['Silica']
        hug.ref = [1]*len(hug.P)
        hug.reflist = ['',newref.reference]
        hug.refidlist = [0,newref.refid[0]]
        hug.rem = [1]*len(hug.P)
        hug.remlist = ['SS304 flyer, driver, impedance match qtz hugoniot in paper']
        DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
        hug.actions = [['Kaitlyn M. Amodeo','create',DateTime_in_ISOFormat,'new entry']]
        hug.diagnostics = ['streak camera','bandpass pyrometer', 'PDV']
        hug.platform = ['gas gun']
        hug.quality = [0]*len(hug.P)
        hug.groups = ['solid compounds and minerals']
    
        shockdat_list.append(hug)
    
        # shot 41 only one with a wiggle
        # shock 25 has poor IR data; good Si diode transit; no streak
    
    

################################################################
print('################################################################')
## Amodeo 2023 Thesis
## REVISED August 2024: updated with additional shots in 2024; see subfolder imp-match-UCD

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Amodeo, K. M., Davies, E. J., Spaulding, D. K., Millot, M., & Stewart, S. T. (n.d.). Transitional Melt in Shock Front: Reinterpreting SiO₂ Shock Temperatures and the Melt Curve. Journal of Geophysical Research: Solid Earth, in prep.']
bibstr = "@article{Amodeo2025, \
    title={Transitional Melt in Shock Front: Reinterpreting SiO₂ Shock Temperatures and the Melt Curve},\
    journal={Journal of Geophysical Research: Solid Earth}, \
    author={Amodeo, Kaitlyn M. and Davies, Erik J. and Spaulding, Dylan K. and Millot, Marius and Stewart, Sarah T.}, \
    pages={in prep.} \
    }"
newref.bibtex = [bibstr]
newref.bibkey = ['Amodeo2025']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

if 1:
    # load previously saved picked data from 
    with open(book_path+'imp-match-UCD/silicadata-N1000.pkl','rb') as f:
        silicadata_list=pickle.load(f)
    
    good_shots = ['shot44','shot27','shot28','shot48','shot14','shot47','shot45','shot33','shot26','shot41','shot50','shot31','shot25','shot23','shot60','shot20','shot22','shot21']
    qtz_good_list = []
    qtz_pgood_list = []
    qtz_good_shotlist = []
    fs_good_list = []
    fs_pgood_list = []
    fs_good_shotlist = []
    qtz_list = []
    fs_list = []
    qtz_p_list = []
    fs_p_list = []
    for i in range(len(silicadata_list)):
        #print(i,silicadata_list[i].material[0],silicadata_list[i].ref[0],silicadata_list[i].P[0],silicadata_list[i].P_std[0],silicadata_list[i].T[0])
        if (silicadata_list[i].material[0] == 'Fused Silica'):
            fs_list.append(silicadata_list[i])
            fs_p_list.append(silicadata_list[i].P[0])
            if silicadata_list[i].ref[0] in good_shots:
                fs_good_list.append(silicadata_list[i])
                fs_pgood_list.append(silicadata_list[i].P[0])
                fs_good_shotlist.append(silicadata_list[i].ref[0])
        else:
            qtz_list.append(silicadata_list[i])
            qtz_p_list.append(silicadata_list[i].P[0])
            if silicadata_list[i].ref[0] in good_shots:
                qtz_good_list.append(silicadata_list[i])
                qtz_pgood_list.append(silicadata_list[i].P[0])
                qtz_good_shotlist.append(silicadata_list[i].ref[0])

    
    #hug_file   = BOXPATH+'Master plotting/IM_solutions_SiO2.txt'

    #IM = pd.read_csv(hug_file, sep=" ")
    hug          = hp.HPeos_Shock()
    hug.matid    = [174]
    hug.material = ['Quartz']
    hug.formula  = ['SiO2']
    hug.phase_ref = ['single crystal']
    hug.orientation = ['z']
    hug.rho_ref = [2.650] # g/cm3
    for i in range(len(silicadata_list)):
        if (silicadata_list[i].material[0] == 'Z-cut Alpha-quartz'):
            #ind = np.where(IM["ID"] == shot_ind_q[i])
            #ind = ind[0]
            hug.exp_id.append(silicadata_list[i].ref[0])
            hug.P.append(silicadata_list[i].P[0])
            hug.P_std.append(silicadata_list[i].P_std[0])
            hug.P_units = ['GPa']
            hug.us.append(silicadata_list[i].us[0])
            hug.us_std.append(silicadata_list[i].us_std[0])
            hug.us_units = ['km/s']
            hug.up.append(silicadata_list[i].up[0])
            hug.up_std.append(silicadata_list[i].up_std[0])
            hug.up_units = ['km/s']
            hug.rho0.append(silicadata_list[i].rho0[0])
            hug.rho0_std.append(silicadata_list[i].rho0_std[0])
            hug.rho0_units = ['g/cm3']
            hug.rho.append(silicadata_list[i].rho[0])
            hug.rho_std.append(silicadata_list[i].rho_std[0])
            hug.rho_units = ['g/cm3']
            hug.T.append(silicadata_list[i].T[0])
            hug.T_std.append(silicadata_list[i].T_std[0])
            hug.T_source.append(silicadata_list[i].T_source[0]) 
            hug.T_units = ['K']
            hug.T_rel.append(silicadata_list[i].T_rel[0]) 
            hug.T_rel_std.append(silicadata_list[i].T_rel_std[0]) 
            hug.T_rel_source.append(silicadata_list[i].T_rel_source[0]) 
            hug.T_rel_units = ['K']
            P = silicadata_list[i].P[0]
            hug.S.append(1793+22.4*P-0.0426*P*P+4.319e-5*P*P*P-1.742e-8*np.power(P,4)) # need to use revised form ## XXXX STS


    hug.T = np.asarray(hug.T)
    hug.T_std = np.asarray(hug.T_std)
    hug.P = np.asarray(hug.P)
    hug.P_std = np.asarray(hug.P_std)
    hug.rho = np.asarray(hug.rho)
    hug.rho_std = np.asarray(hug.rho_std)
    hug.up = np.asarray(hug.up)
    hug.up_std = np.asarray(hug.up_std)
    hug.us = np.asarray(hug.us)
    hug.us_std = np.asarray(hug.us_std)
    hug.exp_id = np.asarray(hug.exp_id)
    
    hug.collections = ['Silica']
    hug.m   = [1.0]*len(hug.P) # nonporous
    hug.ref = [1]*len(hug.P)
    hug.reflist = ['',newref.reference]
    hug.refidlist = [0,newref.refid[0]]
    hug.rem = [1]*len(hug.P)
    hug.remlist = ['SS304 flyer, driver, impedance match qtz hugoniot in paper']
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    hug.actions = [['Kaitlyn M. Amodeo','create',DateTime_in_ISOFormat,'new entry']]
    hug.diagnostics = ['streak camera','bandpass pyrometer', 'PDV']
    hug.platform = ['gas gun']
    hug.quality = [0]*len(hug.P)
    hug.groups = ['solid compounds and minerals']

    shockdat_list.append(hug)

    print('Amodeo 2024 alpha-quartz:')
    print(hug.rem)
    print(hug.exp_id)
    print(hug.P)
    print(hug.T)
    print(hug.T_rel)
    print(hug.T_rel_std)
    print(hug.rho0)
    print(hug.rho0_std)

    # FUSED QUARTZ
    hug          = hp.HPeos_Shock()
    hug.matid    = [143]
    hug.material = ['Quartz amorphous']
    hug.formula  = ['SiO2']
    hug.rho_ref = [2.204] # g/cm3
    for i in range(len(silicadata_list)):
        if (silicadata_list[i].material[0] == 'Fused Silica'):
            hug.exp_id.append(silicadata_list[i].ref[0])
            hug.P.append(silicadata_list[i].P[0])
            hug.P_std.append(silicadata_list[i].P_std[0])
            hug.P_units = ['GPa']
            hug.us.append(silicadata_list[i].us[0])
            hug.us_std.append(silicadata_list[i].us_std[0])
            hug.us_units = ['km/s']
            hug.up.append(silicadata_list[i].up[0])
            hug.up_std.append(silicadata_list[i].up_std[0])
            hug.up_units = ['km/s']
            hug.rho0.append(silicadata_list[i].rho0[0])
            hug.rho0_std.append(silicadata_list[i].rho0_std[0])
            hug.rho0_units = ['g/cm3']
            hug.rho.append(silicadata_list[i].rho[0])
            hug.rho_std.append(silicadata_list[i].rho_std[0])
            hug.rho_units = ['g/cm3']
            hug.T.append(silicadata_list[i].T[0])
            hug.T_std.append(silicadata_list[i].T_std[0])
            hug.T_source.append(silicadata_list[i].T_source[0]) 
            hug.T_units = ['K']
            hug.T_rel.append(silicadata_list[i].T_rel[0]) 
            hug.T_rel_std.append(silicadata_list[i].T_rel_std[0]) 
            hug.T_rel_source.append(silicadata_list[i].T_rel_source[0]) 
            hug.T_rel_units = ['K']
            P = silicadata_list[i].P[0]
            hug.S.append(2294+26.82*P-0.0723*P*P) # need to use revised form ## XXXX STS
            if silicadata_list[i].ref[0] != 'shot50':
                hug.rem.append(1)      # SS304 IM
            else:
                hug.rem.append(2)      # Cu101 IM       

    hug.T = np.asarray(hug.T)
    hug.T_std = np.asarray(hug.T_std)
    hug.exp_id = np.asarray(hug.exp_id)
    hug.P = np.asarray(hug.P)
    hug.P_std = np.asarray(hug.P_std)
    hug.rho = np.asarray(hug.rho)
    hug.rho_std = np.asarray(hug.rho_std)
    hug.up = np.asarray(hug.up)
    hug.up_std = np.asarray(hug.up_std)
    hug.us = np.asarray(hug.us)
    hug.us_std = np.asarray(hug.us_std)
    hug.collections = ['Silica']
    hug.m      = [1.0]*len(hug.P)
    hug.ref = [1]*len(hug.P)
    hug.reflist = ['',newref.reference]
    hug.refidlist = [0,newref.refid[0]]
    hug.remlist = ['SS304 flyer, driver, impedance match qtz hugoniot in paper','Cu101 flyer, driver, impedance match qtz hugoniot in paper']
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    hug.actions = [['Kaitlyn M. Amodeo','create',DateTime_in_ISOFormat,'new entry']]
    hug.diagnostics = ['streak camera','bandpass pyrometer', 'PDV']
    hug.platform = ['gas gun']
    hug.quality = [0]*len(hug.P)
    hug.groups = ['solid compounds and minerals']

    shockdat_list.append(hug)
    
    print('Amodeo 2024 fused silica:')
    print(hug.rem)
    print(hug.exp_id)
    print(hug.P)
    print(hug.T)
    print(hug.T_rel)
    print(hug.T_rel_std)
    print(hug.rho0)
    print(hug.rho0_std)


################################################################
print('################################################################')
## Wackerle 1962 single wave shock Hugoniot data


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Wackerle, J., 1962. Shock‐wave compression of quartz. Journal of applied physics, 33(3), pp.922-937.']
newref.doi = ['10.1063/1.1777192']
bibstr = "@article{Wackerle1962,\
  doi = {10.1063/1.1777192},\
  url = {https://doi.org/10.1063/1.1777192},\
  year = {1962},\
  month = mar,\
  publisher = {{AIP} Publishing},\
  volume = {33},\
  number = {3},\
  pages = {922--937},\
  author = {Jerry Wackerle},\
  title = {Shock-Wave Compression of Quartz},\
  journal = {Journal of Applied Physics}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Wackerle1962']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Wackerle1962 single wave data only
# first wave shock vel 0.5% error
# second wave shock vel 2-4% error (improved over time)
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Wackerle-1962-hug.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.650] # g/cm3
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[0],max_rows=11) 
hug.us_std = 0.005*hug.us # estimated 0.5% error on 1st wave arrival velocity
hug.us_units = ['km/s']
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[1],max_rows=11) 
hug.up_std  = 0.01*hug.up # assumed 1% error on IM solution
hug.up_units = ['km/s']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2],max_rows=11)*0.1 # kbar->GPa
hug.P_std  = np.sqrt(hug.P*hug.P*(.01**2+.005**2))
hug.P_units = ['GPa']
v = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3],max_rows=11) # V/V0
hug.rho      = hug.rho_ref[0]/v
hug.rho_std  = np.sqrt(hug.rho*hug.rho*(.01**2+.005**2)) # estimated 0.5% error on 1st wave velocity; assume 1% on IM up
hug.rho_units = ['g/cm3']
hug.orientation = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4],max_rows=11,dtype='str')
print(hug.P)
print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera']
hug.platform = ['gas gun']
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

# Wackerle1962 single wave data only
# Fused Quartz
hug_file   = silicadatadir+'Wackerle-1962-hug.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[0],max_rows=6) 
hug.us_units = ['km/s']
hug.us_std = 0.005*hug.us # estimated 0.5% error on 1st wave arrival velocity
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[1],max_rows=6) 
hug.up_units = ['km/s']
hug.up_std  = 0.01*hug.up # assumed 1% error on IM solution
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[2],max_rows=6)*0.1 # kbar->GPa
hug.P_units = ['GPa']
hug.P_std  = np.sqrt(hug.P*hug.P*(.01**2+.005**2))
v = np.loadtxt(hug_file,delimiter=',',skiprows=14,usecols=[3],max_rows=6) # V/V0
hug.rho      = hug.rho_ref[0]/v
hug.rho_units = ['g/cm3']
hug.rho_std  = np.sqrt(hug.rho*hug.rho*(.01**2+.005**2)) # estimated 0.5% error on 1st wave velocity; assume 1% on IM up
print(hug.P)
print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

################################################################
print('################################################################')


## Luo et al. 2003 JGR Polymorphism silica
# Table 2: the Stishovite region model for alpha-quartz


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Luo, S.N., Ahrens, T.J. and Asimow, P.D., 2003. Polymorphism, superheating, and amorphization of silica upon shock wave loading and release. Journal of Geophysical Research: Solid Earth, 108(B9).']
newref.doi = ['10.1029/2002JB002317']
bibstr = "@article{Luo2003,\
  doi = {10.1029/2002jb002317},\
  url = {https://doi.org/10.1029/2002jb002317},\
  year = {2003},\
  month = sep,\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {108},\
  number = {B9},\
  author = {Sheng-Nian Luo and Thomas J. Ahrens and Paul D. Asimow},\
  title = {Polymorphism,  superheating,  and amorphization of silica upon shock wave loading and release},\
  journal = {Journal of Geophysical Research: Solid Earth}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Luo2003']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Luo et al. 2003 stishovite region calculation on alpha-quartz Hugoniot shock temp and release temp
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Luo2003-silica-model.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.650] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=8,usecols=[2])
hug.P_units = ['GPa']
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=8,usecols=[6])
hug.T_units = ['K']
hug.T_rel      = np.loadtxt(hug_file,delimiter=',',skiprows=8,usecols=[8])
hug.T_rel_units = ['K']
hug.S      = np.loadtxt(hug_file,delimiter=',',skiprows=8,usecols=[9])
hug.S_units = ['J/K/g']
v = np.loadtxt(hug_file,delimiter=',',skiprows=8,usecols=[3]) # V
hug.rho      = 1./v
hug.rho_units = ['g/cm3']
print(hug.P)
print(hug.T)
print(hug.T_rel)
print(hug.S)
print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.method = ['thermodynamic']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


################################################################
print('################################################################')


## Root et al. 2019 Fused Silica Hugoniot
# Root, S., Townsend, J.P. and Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. Journal of Applied Physics, 126(16), p.165901.


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Root, S., Townsend, J.P. and Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. Journal of Applied Physics, 126(16), p.165901.']
newref.doi = ['10.1063/1.5126205']
bibstr = "@article{Root2019,\
  doi = {10.1063/1.5126205},\
  url = {https://doi.org/10.1063/1.5126205},\
  year = {2019},\
  month = oct,\
  publisher = {{AIP} Publishing},\
  volume = {126},\
  number = {16},\
  pages = {165901},\
  author = {Seth Root and Joshua P. Townsend and Marcus D. Knudson},\
  title = {Shock compression of fused silica: An impedance matching standard},\
  journal = {Journal of Applied Physics}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Root2019']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Root, S., Townsend, J.P. and Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. Journal of Applied Physics, 126(16), p.165901.
# Fused Quartz
hug_file   = silicadatadir+'Root-2019-fused-silica.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.2] # g/cm3
hug.uflyer      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[2]) 
hug.uflyer_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[3]) 
hug.uflyer_units = ['km/s']
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[4]) 
hug.up_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[5]) 
hug.up_units = ['km/s']
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[6]) 
hug.us_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[7]) 
hug.us_units = ['km/s']
hug.rho      =  np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[8]) 
hug.rho_std      =  np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[9]) 
hug.rho_units = ['g/cm3']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[10])
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[11])
hug.P_units = ['GPa']
#print(hug.P)
#print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Al flyer']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR']
hug.platform = ['Z machine']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)



################################################################
print('################################################################')

## Tracy et al. 2018 XRD Fused Silica


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Tracy, S.J., Turneaure, S.J. and Duffy, T.S., 2018. In situ x-ray diffraction of shock-compressed fused silica. Physical review letters, 120(13), p.135702.']
newref.doi = ['10.1103/PhysRevLett.120.135702']
bibstr = "@article{Tracy2018,\
  doi = {10.1103/physrevlett.120.135702},\
  url = {https://doi.org/10.1103/physrevlett.120.135702},\
  year = {2018},\
  month = mar,\
  publisher = {American Physical Society ({APS})},\
  volume = {120},\
  number = {13},\
  author = {Sally June Tracy and Stefan J. Turneaure and Thomas S. Duffy},\
  title = {In situ X-Ray Diffraction of Shock-Compressed Fused Silica},\
  journal = {Physical Review Letters}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Tracy2018']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Tracy, S.J., Turneaure, S.J. and Duffy, T.S., 2018. In situ x-ray diffraction of shock-compressed fused silica. Physical review letters, 120(13), p.135702.
# Fused Quartz
hug_file   = silicadatadir+'Tracy2018-xrd-fused-silica.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.20] # g/cm3 # in SOM as 2.20+-0.01
hug.uflyer          = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[0]) 
hug.uflyer_std      = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[1]) 
hug.uflyer_units    = ['km/s']
hug.P               = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[2])
hug.P_std           = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[3])
hug.P_units         = ['GPa']
hug.rho             = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[4])
hug.rho_std         = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[5])
hug.rho_units       = ['g/cm3']
hug.phase_id        = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[6],dtype='str')
# calculate Us and Up
hug.up = np.sqrt(hug.P*((1/hug.rho_ref[0]) - (1/hug.rho)))
A = hug.P*((1/hug.rho_ref[0]) - (1/hug.rho))
hug.up_std = np.sqrt( hug.P_std*hug.P_std * np.pow(((1/hug.rho_ref[0]) - (1/hug.rho))/(2*np.sqrt(A)),2) + \
                      hug.rho_std*hug.rho_std * np.pow( (hug.P/(2*hug.rho*hug.rho*np.sqrt(A))), 2) )
hug.us = (1/hug.rho_ref[0])*np.sqrt(hug.P/((1/hug.rho_ref[0]) - (1/hug.rho)))
hug.us_std = np.sqrt( hug.P_std*hug.P_std*np.pow((1/(2.*np.sqrt(hug.P/((1/hug.rho_ref[0]) - (1/hug.rho))))),2) + \
                      hug.rho_std*hug.rho_std*np.pow( (1/hug.rho_ref[0])/(2.*hug.rho*hug.rho*np.sqrt(hug.P/((1/hug.rho_ref[0]) - (1/hug.rho)))),2) )
print(hug.P)
print(hug.rho)
print(hug.phase_id)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Al flyer']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['XRD','VISAR']
hug.platform = ['gas gun','APS-DCS']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

################################################################
print('################################################################')

## Tracy et al. 2020 shock XRD quartz

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Tracy, S.J., Turneaure, S.J. and Duffy, T.S., 2020. Structural response of α-quartz under plate-impact shock compression. Science Advances, 6(35), p.eabb3913.']
newref.doi = ['10.1126/sciadv.abb3913']
bibstr = "@article{Tracy2020,\
  doi = {10.1126/sciadv.abb3913},\
  url = {https://doi.org/10.1126/sciadv.abb3913},\
  year = {2020},\
  month = aug,\
  publisher = {American Association for the Advancement of Science ({AAAS})},\
  volume = {6},\
  number = {35},\
  author = {Sally June Tracy and Stefan J. Turneaure and Thomas S. Duffy},\
  title = {Structural response of $\alpha$-quartz under plate-impact shock compression},\
  journal = {Science Advances} }"
newref.bibtex = [bibstr]
newref.bibkey = ['Tracy2020']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])



# Tracy et al. 2020 Sci. Adv. APS XRD
# CRYSTALLINE QUARTZ
hug_file   = silicadatadir+'Tracy2020-xrd-quartz-v2.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref  = [2.650] # g/cm3
hug.exp_id        = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[0],dtype='str')
hug.upflyer       = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[2])
hug.upflyer_std   = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[3])
hug.upflyer_units = ['km/s']
hug.P        = np.asarray(np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[4]))
hug.P_std    = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[5])
hug.P_units  = ['GPa']
hug.rho       =  np.asarray(np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[6]))
hug.rho_std   = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[7])
hug.rho_units = ['g/cm3']
# calculate Us and Up
hug.up = np.sqrt(hug.P*((1/hug.rho_ref[0]) - (1/hug.rho)))
hug.us = (1/hug.rho_ref[0])*np.sqrt(hug.P/((1/hug.rho_ref[0]) - (1/hug.rho)))
print(hug.P)
print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.method = ['thermodynamic']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)



################################################################
print('################################################################')

## Sugiura et al. 1980 -- One fused silica shock temperature
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Sugiura, H., Kondo, K.I. and Sawaoka, A., 1980. Method for the measurement of temperature in shock compression of solids. Review of Scientific Instruments, 51(6), pp.750-752.']
newref.doi = ['10.1063/1.1136296']
bibstr = "@article{Sugiura1980,\
doi = {10.1063/1.1136296},\
url = {https://doi.org/10.1063%2F1.1136296},\
year = 1980,\
month = {jun},\
publisher = {{AIP} Publishing},\
volume = {51},\
number = {6},\
pages = {750--752},\
author = {Hiroshi Sugiura and Ken-ichi Kondo and Akira Sawaoka},\
title = {Method for the measurement of temperature in shock compression of solids},\
journal = {Review of Scientific Instruments}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Sugiura1980']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])
# Sugiura, H., Kondo, K.I. and Sawaoka, A., 1980. Method for the measurement of temperature in shock compression of solids. Review of Scientific Instruments, 51(6), pp.750-752
# Fused Quartz
hug_file   = silicadatadir+'Sugiura-1980-shock-temp-fs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P               = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[0])])
hug.P_std           = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[1])])
hug.P_units         = ['GPa']
hug.T               = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[2])])
hug.T_std           = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[3])])
hug.T_units         = ['K']
hug.rho             = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[4])])
hug.rho_std         = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[5])])
hug.rho_units       = ['g/cm3']
hug.us             = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[7])])
hug.us_std         = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[7])])*.01
hug.us_units       = ['g/cm3']

# calculate Us and Up
hug.up = hug.P/hug.rho_ref/hug.us
#hug.us = (1/hug.rho_ref[0])*np.sqrt(hug.P/((1/hug.rho_ref[0]) - (1/hug.rho)))

print(hug.P)
print(hug.T)
print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Tungsten flyer -> copper driver -> fused silica; reported copper shock 6.1 km/s and fs shock 6.1 km/s']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['pyrometer']
hug.platform = ['gas gun']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


################################################################
print('################################################################')
## Omega data on alpha-quartz T-Us from Millot and Davies
# add reference to reference list
omegapmin = 74.0 # GPa

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Millot and Davies, unpublished']
newref.bibkey = ['MDxx Omega']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

### SHOT 1

# Millot and Davies decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE QUARTZ
#hug_file   = silicadatadir+'Millot-Davies/70811_TandUs.txt'
hug_file   = silicadatadir+'Davies-hires/70811_TandUs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.phase_ref = ['single crystal']
hug.orientation = ['z']
hug.rho_ref  = [2.650] # g/cm3
hug.exp_id        = [70811]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[0])*1000. # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])*1000. # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[2])
hug.us_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[3])
hug.us_units  = ['km/s']
#hug.P = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.P) # GPa
#hug.up = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.up) # km/s
#hug.rho = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.rho) # g/cm3
#hug.P_units = ['GPa']
#hug.up_units = ['km/s']
#hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','Omega EPa']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

#M70811liq = deepcopy(hug) # keep a copy of the liquid Us-up calibration

## call function to converts Us-T to PVT taking into account solid phases
#hug.qtzhug(qtz_liq_hug)

shockdat_list.append(hug)
#M70811 = deepcopy(hug)


## shot 2

# Millot and Davies decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE QUARTZ
#hug_file   = silicadatadir+'Millot-Davies/70813_TandUs.txt'
hug_file   = silicadatadir+'Davies-hires/70813_TandUs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.phase_ref = ['single crystal']
hug.orientation = ['z']
hug.rho_ref  = [2.650] # g/cm3
hug.exp_id        = [70813]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[0])*1000. # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])*1000. # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[2])
hug.us_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[3])
hug.us_units  = ['km/s']
#hug.P = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.P) # GPa
#hug.up = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.up) # km/s
#hug.rho = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.rho) # g/cm3
#hug.P_units = ['GPa']
#hug.up_units = ['km/s']
#hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','Omega EPa']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

#M70813liq = deepcopy(hug) # keep a copy of the liquid Us-up calibration

## call function to converts Us-T to PVT taking into account solid phases
#hug.qtzhug(qtz_liq_hug)

shockdat_list.append(hug)
#M70813 = deepcopy(hug)

## Shot 3

# Millot and Davies decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE QUARTZ
#hug_file   = silicadatadir+'Millot-Davies/70817_TandUs.txt'
hug_file   = silicadatadir+'Davies-hires/70817_TandUs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.phase_ref = ['single crystal']
hug.orientation = ['z']
hug.rho_ref  = [2.650] # g/cm3
hug.exp_id        = [70817]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[0])*1000. # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])*1000. # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[2])
hug.us_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[3])
hug.us_units  = ['km/s']
#hug.P = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.P) # GPa
#hug.up = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.up) # km/s
#hug.rho = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.rho) # g/cm3
#hug.P_units = ['GPa']
#hug.up_units = ['km/s']
#hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','Omega EPa']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

#M70817liq = deepcopy(hug) # keep a copy of the liquid Us-up calibration

## call function to converts Us-T to PVT taking into account solid phases
#hug.qtzhug(qtz_liq_hug)

shockdat_list.append(hug)
#M70817  = deepcopy(hug)

## shot 4

# Millot and Davies decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE QUARTZ
#hug_file   = silicadatadir+'Millot-Davies/70820_TandUs.txt'
hug_file   = silicadatadir+'Davies-hires/70820_TandUs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.phase_ref = ['single crystal']
hug.orientation = ['z']
hug.rho_ref  = [2.650] # g/cm3
hug.exp_id        = [70820]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[0])*1000. # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])*1000. # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[2])
hug.us_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[3])
hug.us_units  = ['km/s']
#hug.P = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.P) # GPa
#hug.up = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.up) # km/s
#hug.rho = np.interp(hug.us,qtz_hug_KD.us,qtz_hug_KD.rho) # g/cm3
#hug.P_units = ['GPa']
#hug.up_units = ['km/s']
#hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','Omega EPa']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

#M70820liq = deepcopy(hug) # keep a copy of the liquid Us-up calibration
## call function to converts Us-T to PVT taking into account solid phases
#hug.qtzhug(qtz_liq_hug)

shockdat_list.append(hug)
#M70820  = deepcopy(hug)



################################################################
print('################################################################')

## Shaner 1981 McQueen and Fritz 1979, 1993 Shock Temperatures
#@article{mcqueen1979temperature,
#  title={The temperature of shocked {SiO$_2$}},
#  author={McQueen, RG and Fritz, JN and Hopson, JW},
#  journal={EOS},
#  volume={60},
#  pages={952},
#  year={1979}
#}
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Shaner, J.W., 1981. New ways of looking for phase transitions at multi-megabar dynamic pressures (No. LA-UR-81-2082; CONF-810828-4). Los Alamos National Lab.(LANL), Los Alamos, NM (United States).']
bibstr = '@techreport{Shaner1981,\
  title={New ways of looking for phase transitions at multi-megabar dynamic pressures},\
  author={Shaner, J. W.},\
  year={1981},\
  number={LA-UR-81-2082},\
  institution={Los Alamos National Lab.(LANL), Los Alamos, NM (United States)}\
}'
newref.bibtex = [bibstr]
newref.bibkey = ['Shaner1981']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# McQueen and Fritz
# Fused Quartz
hug_file   = silicadatadir+'Shaner-1981-sio2-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref  = [2.650] # g/cm3
hug.P               = np.asarray(np.loadtxt(hug_file,delimiter=',',skiprows=20,usecols=[0],max_rows=8))
hug.P_units         = ['GPa']
hug.T               = np.asarray(np.loadtxt(hug_file,delimiter=',',skiprows=20,usecols=[1],max_rows=8))
hug.T_units         = ['K']

print(hug.P)
print(hug.T)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Tungsten flyer -> copper driver -> fused silica; reported copper shock 6.1 km/s and fs shock 6.1 km/s']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['pyrometer']
hug.platform = ['gas gun']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

# McQueen and Fritz
# Fused Quartz
hug_file   = silicadatadir+'Shaner-1981-sio2-shock-temp.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P               = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[0],max_rows=15)
hug.P_units         = ['GPa']
hug.T               = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[1],max_rows=15)
hug.T_units         = ['K']
print(hug.P)
print(hug.T)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Tungsten flyer -> copper driver -> fused silica; reported copper shock 6.1 km/s and fs shock 6.1 km/s']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['pyrometer']
hug.platform = ['gas gun']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


########################
## Sugiura et al. 1982 -- Several fused silica shock temperature
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Sugiura, H., Kondo, K. and Sawaoka, A., 1982. Shock temperatures in fused silica measured by optical technique. Journal of Applied Physics, 53(6), pp.4512-4514.']
newref.doi = ['10.1063/1.331191']
bibstr = "@article{Sugiura1982,\
  title={Shock temperatures in fused silica measured by optical technique},\
  author={Sugiura, H and Kondo, K and Sawaoka, A},\
  journal={Journal of Applied Physics},\
  volume={53},\
  number={6},\
  pages={4512--4514},\
  year={1982},\
  doi={10.1063/1.331191},\
  publisher={American Institute of Physics}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Sugiura1982']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])
# Sugiura, H., Kondo, K. and Sawaoka, A., 1982. Shock temperatures in fused silica measured by optical technique. Journal of Applied Physics, 53(6), pp.4512-4514.
# Fused Quartz
hug_file   = silicadatadir+'Sugiura-1982-shock-temp-fs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P               = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[0])
hug.P_units         = ['GPa']
hug.T               = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[1])
hug.T_units         = ['K']
hug.T_b               = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[2])
hug.T_b_units         = ['K']
print(hug.P)
print(hug.T)
print(hug.T_b)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['No impact conditions reported; only pressure and temperature. T is 620 nm; T_b is 810 nm.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['pyrometer']
hug.platform = ['gas gun']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

########################
## Berryman et al. 4 fused silica Hugoniot Points
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Berryman, E.J., Winey, J.M., Gupta, Y.M. and Duffy, T.S., 2019. Sound Velocities in Shock‐Synthesized Stishovite to 72 GPa. Geophysical Research Letters, 46(23), pp.13695-13703.']
newref.doi = ['10.1029/2019GL085301']
bibstr = "@article{Berryman_2019,\
  title={Sound Velocities in Shock-Synthesized Stishovite to 72 GPa},\
  author={Berryman, Eleanor J and Winey, J Michael and Gupta, Yogendra M and Duffy, Thomas S},\
  journal={Geophysical Research Letters},\
  volume={46},\
  number={23},\
  pages={13695--13703},\
  year={2019},\
  publisher={Wiley Online Library},\
  doi={10.1029/2019GL085301}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Berryman_2019']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

# Fused Quartz
hug_file   = silicadatadir+'Berryman-2019-fused-silica-hug.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.201] # g/cm3
hug.us          = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3])
hug.us_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4])
hug.us_units    = ['km/s']
hug.up          = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[5])
hug.up_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[6])
hug.up_units    = ['km/s']
hug.P          = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[7])
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[8])
hug.P_units    = ['GPa']
hug.rho          = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[9])
hug.rho_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[10])
hug.rho_units    = ['g/cm3']

print(hug.P)
print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR']
hug.platform = ['gas gun']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

########################
## Knudson et al. 2006 Quartz Hugoniot Points from Z
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Knudson, M.D. and Desjarlais, M.P., 2013. Adiabatic release measurements in α-quartz between 300 and 1200 GPa: Characterization of α-quartz as a shock standard in the multimegabar regime. Physical Review B, 88(18), p.184107.']
newref.doi = ['10.1103/PhysRevB.88.184107']
bibstr = "@article{Knudson_2013,\
  title = {Adiabatic release measurements in α-quartz between 300 and 1200 GPa: Characterization of α-quartz as a shock standard in the multimegabar regime},\
  author = {Knudson, M. D. and Desjarlais, M. P.},\
  journal = {Phys. Rev. B},\
  volume = {88},\
  issue = {18},\
  pages = {184107},\
  numpages = {18},\
  year = {2013},\
  month = {Nov},\
  publisher = {American Physical Society},\
  doi = {10.1103/PhysRevB.88.184107},\
  url = {https://link.aps.org/doi/10.1103/PhysRevB.88.184107}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Knudson_2013']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

# Fused Quartz
hug_file   = silicadatadir+'Quartz_PRB88_and_JAP122.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.65] # g/cm3
hug.us          = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[4])
hug.us_std      = np.sqrt(np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[6]))
hug.us_units    = ['km/s']
hug.up          = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[5])
hug.up_std      = np.sqrt(np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[7]))
hug.up_units    = ['km/s']
hug.P          = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[9])
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[10])
hug.P_units    = ['GPa']
hug.rho          = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[11])
hug.rho_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[12])
hug.rho_units    = ['g/cm3']

#print(hug.P)
#print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR']
hug.platform = ['Z machine']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

#====================================================================================

################################################################
print('################################################################')
## Boslough JGR 1988 Postshock temperatures in silica WITH REVISED PRESSURES

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Boslough, M.B., 1988. Postshock temperatures in silica. Journal of Geophysical Research: Solid Earth, 93(B6), pp.6477-6484.']
newref.doi = ['10.1029/JB093iB06p06477']
bibstr = "@article{Boslough1988,\
  doi = {10.1029/jb093ib06p06477},\
  url = {https://doi.org/10.1029/jb093ib06p06477},\
  year = {1988},\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {93},\
  number = {B6},\
  pages = {6477},\
  author = {M. B. Boslough},\
  title = {Postshock temperatures in silica},\
  journal = {Journal of Geophysical Research}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Boslough1988']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'revised pressures']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Boslough 1988 CRYSTALLINE Quartz
hug_file     = silicadatadir+'Boslough-1988-silica-shock-temp-release-revised-P.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref = [2.65] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[0],max_rows=9) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[7],max_rows=9) 
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[1],max_rows=9) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[2],max_rows=9) 
hug.T_units = ['K']
hug.T_rel      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[3],max_rows=9) 
hug.T_rel_std  = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[4],max_rows=9) 
hug.T_rel_units = ['K']

hug.S      = np.zeros(len(hug.T)) 

hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[8],max_rows=9) 
hug.us_units = ['km/s']
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[9],max_rows=9) 
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[10],max_rows=9) 
hug.up_units = ['km/s']
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[11],max_rows=9) 
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[12],max_rows=9) /1000
hug.rho_std_units = ['g/cm3']
hug.orientation      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[13],max_rows=9,dtype='str') 
print(hug.P)
print(hug.T)
print(hug.T_rel)


hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Reanalysis of Lyzenga et al. 1983 with added post-shock temperatures.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


# Boslough 1988 AMORPHOUS QUARTZ
hug_file     = silicadatadir+'Boslough-1988-silica-shock-temp-release-revised-P.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[0],max_rows=7) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[7],max_rows=7) 
hug.T      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[1],max_rows=7) 
hug.T_std  = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[2],max_rows=7) 
hug.T_units = ['K']
hug.T_rel      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[3],max_rows=7) 
hug.T_rel_std  = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[4],max_rows=7) 
hug.T_rel_units = ['K']

hug.S      = np.zeros(len(hug.T)) 

hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[8],max_rows=7) 
hug.us_units = ['km/s']
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[9],max_rows=7) 
hug.us_std_units = ['km/s']
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[10],max_rows=7) 
hug.up_units = ['km/s']
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[11],max_rows=7) 
hug.up_std_units = ['km/s']
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[12],max_rows=7) /1000
hug.rho_std_units = ['g/cm3']
hug.orientation      = np.loadtxt(hug_file,delimiter=',',skiprows=10,usecols=[13],max_rows=7,dtype='str') 
print(hug.P)
print(hug.T)
print(hug.T_rel)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Reanalysis of Lyzenga et al. 1983 with added post-shock temperatures.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'revised pressures']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

############################################
print('Jackson, I., & Ahrens, T. J. (1979). Shock wave compression of single-crystal forsterite. Journal of Geophysical Research: Solid Earth, 84(B6), 3039–3048. https://doi.org/10.1029/JB084iB06p03039')


maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Jackson, I., & Ahrens, T. J. (1979). Shock wave compression of single-crystal forsterite. Journal of Geophysical Research: Solid Earth, 84(B6), 3039–3048. https://doi.org/10.1029/JB084iB06p03039']
newref.doi = ['10.1029/JB084iB06p03039']
bibstr = "@article{Jackson1979,\
  title={Shock wave compression of single-crystal forsterite},\
  author={Jackson, Ian and Ahrens, Thomas J},\
  journal={Journal of Geophysical Research: Solid Earth},\
  volume={84},\
  number={B6},\
  pages={3039--3048},\
  year={1979},\
  publisher={Wiley Online Library},\
  doi={10.1029/JB084iB06p03039}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Jackson1979']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])



# Jackson et al. 1979 
hug_file     = silicadatadir+'Jackson-1979-fs-hug.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.20] # g/cm3
hug.vi      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2],max_rows=6) 
hug.vi_units = ['km/s']
hug.vi_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3],max_rows=6) 
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[8],max_rows=6) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[9],max_rows=6) 
hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4],max_rows=6) 
hug.us_units = ['km/s']
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[5],max_rows=6) 
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[6],max_rows=6) 
hug.up_units = ['km/s']
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[7],max_rows=6) 
hug.up_units = ['km/s']
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[10],max_rows=6) 
hug.rho_std  = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[11],max_rows=6) 
hug.rho_units = ['g/cm3']
print(hug.P)
print(hug.T)
print(hug.T_rel)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

print('################################################################')
## NEW IMPEDANCE MATCH FOR JACKSON
# add reference to reference list
    
silicadata_list=[]
# load previously saved picked data from 
with open(book_path+'imp-match-JA79/jacksondata-N1000.pkl','rb') as f:
    silicadata_list=pickle.load(f)

# FUSED QUARTZ
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
for i in range(len(silicadata_list)):
    if (silicadata_list[i].material[0] == 'Fused Silica'):
        print('here Jackson FS')
        #hug.phase_ref.append(['fs'])
        hug.exp_id.append(silicadata_list[i].ref[0])
        hug.P.append(silicadata_list[i].P[0])
        hug.P_std.append(silicadata_list[i].P_std[0])
        hug.P_units = ['GPa']
        hug.us.append(silicadata_list[i].us[0])
        hug.us_std.append(silicadata_list[i].us_std[0])
        hug.us_units = ['km/s']
        hug.up.append(silicadata_list[i].up[0])
        hug.up_std.append(silicadata_list[i].up_std[0])
        hug.up_units = ['km/s']
        #hug.rho0.append(silicadata_list[i].rho0[0])
        #hug.rho0_std.append(silicadata_list[i].rho0_std[0])
        hug.rho0_units = ['g/cm3']
        hug.rho.append(silicadata_list[i].rho[0])
        hug.rho_std.append(silicadata_list[i].rho_std[0])
        hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.m      = [1.0]*len(hug.P)
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.remlist = ['TA flyer, driver']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','microphotometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

print('Jackson & Ahrens 1979 redo fused silica:')
print(hug.rem)
print(hug.exp_id)
print(hug.P)
print(hug.rho)


print('################################################################')
if 1:
    # add reference to reference list
    
    maxrefid=0
    for item in ref_list:
        if int(item.refid[0]) > maxrefid:
            maxrefid = int(item.refid[0])
            
    print('Previous max refid # = ',maxrefid)
    newrefid = maxrefid+1
    print('New refid # = ',newrefid)
    
    newref = hp.HPeos_Reference()
    newref.refid = [newrefid]
    newref.reference = ['Lyzenga et al. 1983 and Boslough 1988 redo impedance match in Amodeo et al. in prep.']
    newref.bibkey = ['L83+B88-STSredoIM']
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
    
    ref_list.append(newref)
    #pprint(list(newref))
    print(newref.reference[0])
    
    if 1:
        silicadata_list=[]
        # load previously saved picked data from 
        with open(book_path+'imp-match-Lyzenga/lyzengadata-N1000.pkl','rb') as f:
            silicadata_list=pickle.load(f)
        print('loaded lyzenga pickle file')
        good_shots = ['shotQ1','shotQ2','shotQ3','shotQ4','shotQ5','shotQ6','shotFS1','shotFS2','shotFS3','shotFS4','shotFS5']
        qtz_good_list = []
        qtz_pgood_list = []
        qtz_good_shotlist = []
        fs_good_list = []
        fs_pgood_list = []
        fs_good_shotlist = []
        qtz_list = []
        fs_list = []
        qtz_p_list = []
        fs_p_list = []
        for i in range(len(silicadata_list)):
            #print(i,silicadata_list[i].material[0],silicadata_list[i].ref[0],silicadata_list[i].P[0],silicadata_list[i].P_std[0],silicadata_list[i].T[0])
            if (silicadata_list[i].material[0] == 'Fused Silica'):
                fs_list.append(silicadata_list[i])
                fs_p_list.append(silicadata_list[i].P[0])
                if silicadata_list[i].ref[0] in good_shots:
                    fs_good_list.append(silicadata_list[i])
                    fs_pgood_list.append(silicadata_list[i].P[0])
                    fs_good_shotlist.append(silicadata_list[i].ref[0])
            else:
                qtz_list.append(silicadata_list[i])
                qtz_p_list.append(silicadata_list[i].P[0])
                if silicadata_list[i].ref[0] in good_shots:
                    qtz_good_list.append(silicadata_list[i])
                    qtz_pgood_list.append(silicadata_list[i].P[0])
                    qtz_good_shotlist.append(silicadata_list[i].ref[0])
    
        
        #hug_file   = BOXPATH+'Master plotting/IM_solutions_SiO2.txt'
    
        #IM = pd.read_csv(hug_file, sep=" ")
        hug          = hp.HPeos_Shock()
        hug.matid    = [174]
        hug.material = ['Quartz']
        hug.formula  = ['SiO2']
        hug.rho_ref = [2.650] # g/cm3
        for i in range(len(silicadata_list)):
            print(i,silicadata_list[i].material[0])
            if (silicadata_list[i].material[0] == 'Z-cut Alpha Quartz'): # or (silicadata_list[i].material[0] == 'Perp Z-cut Alpha-quartz'):
                print('here z')
                #ind = np.where(IM["ID"] == shot_ind_q[i])
                #ind = ind[0]
                hug.orientation.append('z')
                hug.exp_id.append(silicadata_list[i].ref[0])
                hug.P.append(silicadata_list[i].P[0])
                hug.P_std.append(silicadata_list[i].P_std[0])
                hug.P_units = ['GPa']
                hug.us.append(silicadata_list[i].us[0])
                hug.us_std.append(silicadata_list[i].us_std[0])
                hug.us_units = ['km/s']
                hug.up.append(silicadata_list[i].up[0])
                hug.up_std.append(silicadata_list[i].up_std[0])
                hug.up_units = ['km/s']
                hug.rho0_units = ['g/cm3']
                hug.rho.append(silicadata_list[i].rho[0])
                hug.rho_std.append(silicadata_list[i].rho_std[0])
                hug.rho_units = ['g/cm3']
                hug.T.append(silicadata_list[i].T[0])
                hug.T_std.append(silicadata_list[i].T_std[0])
                #hug.T_source.append(silicadata_list[i].T_source[0]) 
                hug.T_units = ['K']
                hug.T_rel.append(silicadata_list[i].T_rel[0]) 
                hug.T_rel_std.append(silicadata_list[i].T_rel_std[0]) 
                #hug.T_rel_source.append(silicadata_list[i].T_rel_source[0]) 
                hug.T_rel_units = ['K']
                P = silicadata_list[i].P[0]
                hug.S.append(1793+22.4*P-0.0426*P*P+4.319e-5*P*P*P-1.742e-8*np.power(P,4)) # need to use revised form ## XXXX STS
            if (silicadata_list[i].material[0] == 'Perp Z-cut Alpha Quartz'):
                print('here pz')
                #ind = np.where(IM["ID"] == shot_ind_q[i])
                #ind = ind[0]
                hug.orientation.append('pz')
                hug.exp_id.append(silicadata_list[i].ref[0])
                hug.P.append(silicadata_list[i].P[0])
                hug.P_std.append(silicadata_list[i].P_std[0])
                hug.P_units = ['GPa']
                hug.us.append(silicadata_list[i].us[0])
                hug.us_std.append(silicadata_list[i].us_std[0])
                hug.us_units = ['km/s']
                hug.up.append(silicadata_list[i].up[0])
                hug.up_std.append(silicadata_list[i].up_std[0])
                hug.up_units = ['km/s']
                hug.rho0_units = ['g/cm3']
                hug.rho.append(silicadata_list[i].rho[0])
                hug.rho_std.append(silicadata_list[i].rho_std[0])
                hug.rho_units = ['g/cm3']
                hug.T.append(silicadata_list[i].T[0])
                hug.T_std.append(silicadata_list[i].T_std[0])
                #hug.T_source.append(silicadata_list[i].T_source[0]) 
                hug.T_units = ['K']
                hug.T_rel.append(silicadata_list[i].T_rel[0]) 
                hug.T_rel_std.append(silicadata_list[i].T_rel_std[0]) 
                #hug.T_rel_source.append(silicadata_list[i].T_rel_source[0]) 
                hug.T_rel_units = ['K']
                P = silicadata_list[i].P[0]
                hug.S.append(1793+22.4*P-0.0426*P*P+4.319e-5*P*P*P-1.742e-8*np.power(P,4)) # need to use revised form ## XXXX STS
    
        hug.orientation = np.asarray(hug.orientation)
        hug.T = np.asarray(hug.T)
        hug.T_std = np.asarray(hug.T_std)
        hug.P = np.asarray(hug.P)
        hug.P_std = np.asarray(hug.P_std)
        hug.us = np.asarray(hug.us)
        hug.us_std = np.asarray(hug.us_std)
        hug.up = np.asarray(hug.up)
        hug.up_std = np.asarray(hug.up_std)
        hug.rho = np.asarray(hug.rho)
        hug.rho_std = np.asarray(hug.rho_std)
        
        hug.collections = ['Silica']
        hug.m   = [1.0]*len(hug.P) # nonporous
        hug.ref = [1]*len(hug.P)
        hug.reflist = ['',newref.reference]
        hug.refidlist = [0,newref.refid[0]]
        hug.rem = [1]*len(hug.P)
        hug.remlist = ['SS304 flyer, driver, impedance match qtz hugoniot in paper']
        DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
        hug.actions = [['Kaitlyn M. Amodeo','create',DateTime_in_ISOFormat,'new entry']]
        hug.diagnostics = ['streak camera','bandpass pyrometer', 'PDV']
        hug.platform = ['gas gun']
        hug.quality = [0]*len(hug.P)
        hug.groups = ['solid compounds and minerals']
    
        shockdat_list.append(hug)
    
        print('Boslough redo quartz:')
        print(hug.rem)
        print(hug.exp_id)
        print(hug.P)
        print(hug.T)
        print(hug.T_rel)
        print(hug.T_rel_std)
        print(hug.rho)
        print(hug.rho_std)
        print(hug.orientation)
    
        # FUSED QUARTZ
        hug          = hp.HPeos_Shock()
        hug.matid    = [143]
        hug.material = ['Quartz amorphous']
        hug.formula  = ['SiO2']
        hug.rho_ref = [2.204] # g/cm3
        for i in range(len(silicadata_list)):
            if (silicadata_list[i].material[0] == 'Fused Silica'):
                print('here FS')
                #hug.phase_ref.append(['fs'])
                hug.exp_id.append(silicadata_list[i].ref[0])
                hug.P.append(silicadata_list[i].P[0])
                hug.P_std.append(silicadata_list[i].P_std[0])
                hug.P_units = ['GPa']
                hug.us.append(silicadata_list[i].us[0])
                hug.us_std.append(silicadata_list[i].us_std[0])
                hug.us_units = ['km/s']
                hug.up.append(silicadata_list[i].up[0])
                hug.up_std.append(silicadata_list[i].up_std[0])
                hug.up_units = ['km/s']
                #hug.rho0.append(silicadata_list[i].rho0[0])
                #hug.rho0_std.append(silicadata_list[i].rho0_std[0])
                hug.rho0_units = ['g/cm3']
                hug.rho.append(silicadata_list[i].rho[0])
                hug.rho_std.append(silicadata_list[i].rho_std[0])
                hug.rho_units = ['g/cm3']
                hug.T.append(silicadata_list[i].T[0])
                hug.T_std.append(silicadata_list[i].T_std[0])
                #hug.T_source.append(silicadata_list[i].T_source[0]) 
                hug.T_units = ['K']
                hug.T_rel.append(silicadata_list[i].T_rel[0]) 
                hug.T_rel_std.append(silicadata_list[i].T_rel_std[0]) 
                #hug.T_rel_source.append(silicadata_list[i].T_rel_source[0]) 
                hug.T_rel_units = ['K']
                P = silicadata_list[i].P[0]
                hug.S.append(2294+26.82*P-0.0723*P*P) # need to use revised form ## XXXX STS

        hug.T = np.asarray(hug.T)
        hug.T_std = np.asarray(hug.T_std)
        hug.P = np.asarray(hug.P)
        hug.P_std = np.asarray(hug.P_std)
        hug.us = np.asarray(hug.us)
        hug.us_std = np.asarray(hug.us_std)
        hug.up = np.asarray(hug.up)
        hug.up_std = np.asarray(hug.up_std)
        hug.rho = np.asarray(hug.rho)
        hug.rho_std = np.asarray(hug.rho_std)
        
        hug.collections = ['Silica']
        hug.m      = [1.0]*len(hug.P)
        hug.ref = [1]*len(hug.P)
        hug.reflist = ['',newref.reference]
        hug.refidlist = [0,newref.refid[0]]
        hug.remlist = ['SS304 flyer, driver, impedance match qtz hugoniot in paper','Cu101 flyer, driver, impedance match qtz hugoniot in paper']
        DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
        hug.actions = [['Kaitlyn M. Amodeo','create',DateTime_in_ISOFormat,'new entry']]
        hug.diagnostics = ['streak camera','bandpass pyrometer', 'PDV']
        hug.platform = ['gas gun']
        hug.quality = [0]*len(hug.P)
        hug.groups = ['solid compounds and minerals']
    
        shockdat_list.append(hug)
        
        print('Boslough redo fused silica:')
        print(hug.rem)
        print(hug.exp_id)
        print(hug.P)
        print(hug.T)
        print(hug.T_rel)
        print(hug.T_rel_std)
        print(hug.rho)
        print(hug.rho_std)


############################################
# FENG ET AL GRL 2024 HUGONIOT POINTS
# this section creates the reference entry, which is reused below for decaying shocks
#
print('HUGONIOT POINTS')
print('Feng, X., Katagiri, K., Qu, J., Nonaka, K., Sun, L., Zhu, P., et al. (2024). Shock Compression of Coesite up to 950 GPa. Geophysical Research Letters, 51(24), e2024GL109873. https://doi.org/10.1029/2024GL109873')

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.doi = ['10.1029/2024GL109873']
newref.reference = ['Feng, X., Katagiri, K., Qu, J., Nonaka, K., Sun, L., Zhu, P., et al. (2024). Shock Compression of Coesite up to 950 GPa. Geophysical Research Letters, 51(24), e2024GL109873. https://doi.org/10.1029/2024GL109873']
bibstr = "@article{Feng2024, \
    title={Shock Compression of Coesite up to 950 GPa}, \
    volume={51}, \
    DOI={10.1029/2024GL109873}, \
    number={24}, \
    journal={Geophysical Research Letters}, \
    author={Feng, Xiaokang and Katagiri, Kento and Qu, Jia and Nonaka, Keita and Sun, Liang and Zhu, Pinwen and Ozaki, Norimasa and Sano, Takayoshi and Sekine, Toshimori and Yang, Wenge}, \
    year={2024}, \
    pages={e2024GL109873}\
    }"
newref.bibtex = [bibstr]
newref.bibkey = ['Feng2024']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

feng_refid = newref.refid[0]
feng_reference = newref.reference

# Feng et al. 2024
hug_file     = silicadatadir+'cst_hug2.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [149]
hug.material = ['Coesite']
hug.formula  = ['SiO2']
hug.rho_ref = [2.92] # g/cm3
hug.rho_ref_std = [0.03] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.phase_ref = ['polycrystalline']
hug.orientation = ['polycrystalline']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[7],max_rows=7) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[8],max_rows=7) 
hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[3],max_rows=7) 
hug.us_units = ['km/s']
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[4],max_rows=7) 
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[5],max_rows=7) 
hug.up_units = ['km/s']
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[6],max_rows=7) 
hug.up_units = ['km/s']
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[9],max_rows=7) 
hug.rho_std  = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[10],max_rows=7) 
hug.rho_units = ['g/cm3']
print(hug.P)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['line VISAR','bandpass pyrometer']
hug.platform = ['GEKKO laser']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


############################################
print('Luo, S.-N., Mosenfelder, J. L., Asimow, P. D., & Ahrens, T. J. (2002). Direct shock wave loading of Stishovite to 235 GPa: Implications for perovskite stability relative to an oxide assemblage at lower mantle conditions. Geophysical Research Letters, 29(14), 36-1-36–4. https://doi.org/10.1029/2002GL015627')


maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Luo, S.-N., Mosenfelder, J. L., Asimow, P. D., & Ahrens, T. J. (2002). Direct shock wave loading of Stishovite to 235 GPa: Implications for perovskite stability relative to an oxide assemblage at lower mantle conditions. Geophysical Research Letters, 29(14), 36-1-36–4. https://doi.org/10.1029/2002GL015627']
newref.doi = ['10.1029/2002GL015627']
bibstr = "@article{Luo_MAA_2002,\
    title={Direct shock wave loading of Stishovite to 235 GPa: Implications for perovskite stability relative to an oxide assemblage at lower mantle conditions}, \
    volume={29}, \
    DOI={10.1029/2002GL015627}, \
    number={14}, \
    journal={Geophysical Research Letters}, \
    author={Luo, Sheng-Nian and Mosenfelder, J. L. and Asimow, P. D. and Ahrens, Thomas J.}, \
    year={2002}, \
    pages={36-1-36–4}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Luo_MAA_2002']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Luo et al. GRL 2002
hug_file     = silicadatadir+'Luo2002-stv-cst.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [149]
hug.material = ['Coesite']
hug.formula  = ['SiO2']
hug.rho_ref = [2.92] # g/cm3
hug.rho_ref_std = [0.03] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[8],max_rows=4) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[9],max_rows=4) 
hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[4],max_rows=4) 
hug.us_units = ['km/s']
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[5],max_rows=4) 
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[6],max_rows=4) 
hug.up_units = ['km/s']
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[7],max_rows=4) 
hug.up_units = ['km/s']
hug.rho0      = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[3],max_rows=4) 
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[10],max_rows=4) 
hug.rho_std  = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[11],max_rows=4) 
hug.rho_units = ['g/cm3']
print(hug.P)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Likely polycrystalline.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

# SAME REFERENCE DIFFERENT MATERIAL - STV
print('Luo, S.-N., Mosenfelder, J. L., Asimow, P. D., & Ahrens, T. J. (2002). Direct shock wave loading of Stishovite to 235 GPa: Implications for perovskite stability relative to an oxide assemblage at lower mantle conditions. Geophysical Research Letters, 29(14), 36-1-36–4. https://doi.org/10.1029/2002GL015627')

# Luo et al. GRL 2002
hug_file     = silicadatadir+'Luo2002-stv-cst.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [672]
hug.material = ['Stishovite']
hug.formula  = ['SiO2']
hug.rho_ref = [4.31] # g/cm3
hug.rho_ref_std = [0.04] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[8],max_rows=4) 
hug.P_units = ['GPa']
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[9],max_rows=4) 
hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[4],max_rows=4) 
hug.us_units = ['km/s']
hug.us_std   = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[5],max_rows=4) 
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[6],max_rows=4) 
hug.up_units = ['km/s']
hug.up_std   = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[7],max_rows=4) 
hug.up_units = ['km/s']
hug.rho0      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[3],max_rows=4) 
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[10],max_rows=4) 
hug.rho_std  = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[11],max_rows=4) 
hug.rho_units = ['g/cm3']
print(hug.P)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


############################################
print('M. A. Podurets, G. V. Simakov, G. S. Telegin, R. F. Trunin, Polymorphism of silica in shock waves and equation of state of coesite and stishovite, Izv. Akad. Nauk SSSR. Fiz. Zemli 1, 16-25 (1981) [in Russian]')


maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['M. A. Podurets, G. V. Simakov, G. S. Telegin, R. F. Trunin, Polymorphism of silica in shock waves and equation of state of coesite and stishovite, Izv. Akad. Nauk SSSR. Fiz. Zemli 1, 16-25 (1981) [in Russian]']
newref.doi = ['']
bibstr = "@article{Podurets1981,\
  title={Polymorphism of silica in shock waves and equation of state of coesite and stishovite},\
  author={Podurets, MA and Simakov, GV and Telegin, GS and Trunin, RF},\
  journal={Izvestiya, Physics of the Solid Earth},\
  volume={17},\
  pages={9--15},\
  year={1981}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Podurets1981']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

hug_file     = silicadatadir+'Podurets-1981-cst-hug.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [149]
hug.material = ['Coesite']
hug.formula  = ['SiO2']
hug.rho_ref = [2.92] # g/cm3
hug.rho_ref_std = [0.03] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[3],max_rows=7) 
hug.P_units = ['GPa']
hug.m      = [1.0]*len(hug.P)
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[2],max_rows=7) 
hug.us_units = ['km/s']
hug.rho      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[5],max_rows=7) 
hug.rho_units = ['g/cm3']
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1],max_rows=7) 
hug.up_units = ['km/s']

hug.us_std  = hug.us*.01 # assumed error in velocity measurements
hug.up_std  = hug.up*.01
hug.P_std  = np.sqrt(hug.P*hug.P*(.01**2+.01**2))
hug.rho_std  = np.sqrt(hug.rho*hug.rho*(.01**2+.01**2))


print(hug.P)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['streak camera','bandpass pyrometer']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)


############################################
print('Jones, A. H., Isbell, W. M., Shipman, F. H., Perkins, R. D., Green, S. J., & Maiden, C. J. (1968). Material Properties Measurements for Selected Materials (Interim Report  Contract NAS2-3427 No. NASA NTRS 68N28631) (p. 56). Warren, Michigan: General Motors Technical Center.')


maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Jones, A. H., Isbell, W. M., Shipman, F. H., Perkins, R. D., Green, S. J., & Maiden, C. J. (1968). Material Properties Measurements for Selected Materials (Interim Report  Contract NAS2-3427 No. NASA NTRS 68N28631) (p. 56). Warren, Michigan: General Motors Technical Center.']
newref.doi = ['']
bibstr = "@report{Jones1968,\
   address={Warren, Michigan},\
    type={Interim Report  Contract NAS2-3427},\
   title={Material Properties Measurements for Selected Materials}, \
  number={NASA NTRS 68N28631},\
    institution={General Motors Technical Center},\
   author={Jones, A. H. and Isbell, W. M. and Shipman, F. H. and Perkins, R. D. and Green, S. J. and Maiden, C. J.},\
  year={1968},\
   month=apr,\
    pages={56} }"
newref.bibtex = [bibstr]
newref.bibkey = ['Jones1968']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

hug_file     = silicadatadir+'Jones1968-fs.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.20] # g/cm3
hug.rho_ref_std = [0.03] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.rho0      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[0])])
hug.rho0_units = ['g/cm3']
hug.P      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3])])
hug.P_units = ['GPa']
hug.P_std  = np.sqrt(hug.P*hug.P*(.01**2+.01**2))
hug.m      = [1.0]*len(hug.P)
hug.us      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[1])])
hug.us_std      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[1])*.01])
hug.us_units = ['km/s']
hug.rho      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4])])
hug.rho_units = ['g/cm3']
hug.rho_std  = np.sqrt(hug.rho*hug.rho*(.01**2+.01**2)) # estimated 0.5% error on 1st wave velocity; assume 1% on IM up
hug.up      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2])])
hug.up_std      = np.asarray([np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2])*.01])
hug.up_units = ['km/s']
print(hug.P)

hug.collections = ['Silica']
hug.ref = [1]*len(hug.P)
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['pins']
hug.platform = ['gas gun']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

################################################################
## COESITE DECAYING SHOCK TEMPERATURE Feng et al. GRL 2024
print('DECAYING SHOCK TEMPERATURES - 3 shots')
print('Feng, X., Katagiri, K., Qu, J., Nonaka, K., Sun, L., Zhu, P., et al. (2024). Shock Compression of Coesite up to 950 GPa. Geophysical Research Letters, 51(24), e2024GL109873. https://doi.org/10.1029/2024GL109873')

# REFERENCE ALREADY CREATED ABOVE: feng_reference, feng_refid

### SHOT 1    
# Feng et al. decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE COESITE
hug_file   = silicadatadir+'cst-43939.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [149]
hug.material = ['Coesite']
hug.formula  = ['SiO2']
hug.phase_ref = ['polycrystalline']
hug.orientation = ['polycrystalline']
hug.rho_ref  = [2.920] # g/cm3
hug.rho_ref_units  = ['g/cm3']
hug.rho_ref_std  = [0.03] # g/cm3
hug.exp_id        = [43941]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[9]) # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[10]) # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])
hug.us_std    = np.zeros(len(hug.us))
hug.us_units  = ['km/s']
if 0:
    hug.up = (hug.us-res[1])/res[0] # km/s
    ind = np.where(hug.us < ustrans1)
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-ress[1])/ress[0] # km/s
    ind = np.where((hug.us > ustrans1) & (hug.us < ustrans2))
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-rest[1])/rest[0] # km/s
    hug.P = hug.rho_ref*hug.up*hug.us # GPa
    hug.rho = 1/((1/hug.rho_ref[0])*(1-hug.up/hug.us)) # g/cm3
    hug.P_units = ['GPa']
    hug.up_units = ['km/s']
    hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',feng_reference]
hug.refidlist = [0,feng_refid]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','GEKKO']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

C43939 = deepcopy(hug) # keep a copy of the liquid Us-up calibration
shockdat_list.append(hug)

### SHOT 2  
# Feng et al. decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE COESITE
hug_file   = silicadatadir+'cst-43941.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [149]
hug.material = ['Coesite']
hug.formula  = ['SiO2']
hug.phase_ref = ['polycrystalline']
hug.orientation = ['polycrystalline']
hug.rho_ref  = [2.920] # g/cm3
hug.rho_ref_units  = ['g/cm3']
hug.rho_ref_std  = [0.03] # g/cm3
hug.exp_id        = [43939]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[9]) # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[10]) # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])
hug.us_std    = np.zeros(len(hug.us))
hug.us_units  = ['km/s']
if 0:
    hug.up = (hug.us-res[1])/res[0] # km/s
    ind = np.where(hug.us < ustrans1)
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-ress[1])/ress[0] # km/s
    ind = np.where((hug.us > ustrans1) & (hug.us < ustrans2))
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-rest[1])/rest[0] # km/s
    hug.P = hug.rho_ref*hug.up*hug.us # GPa
    hug.rho = 1/((1/hug.rho_ref[0])*(1-hug.up/hug.us)) # g/cm3
    hug.P_units = ['GPa']
    hug.up_units = ['km/s']
    hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',feng_reference]
hug.refidlist = [0,feng_refid]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','GEKKO']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

C43941 = deepcopy(hug) # keep a copy of the liquid Us-up calibration
shockdat_list.append(hug)

### SHOT 3
# Feng et al. decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# CRYSTALLINE COESITE
hug_file   = silicadatadir+'cst-43943.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [149]
hug.material = ['Coesite']
hug.formula  = ['SiO2']
hug.phase_ref = ['polycrystalline']
hug.orientation = ['polycrystalline']
hug.rho_ref  = [2.920] # g/cm3
hug.rho_ref_units  = ['g/cm3']
hug.rho_ref_std  = [0.03] # g/cm3
hug.exp_id        = [43943]
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[9]) # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[10]) # K
hug.T_units  = ['K']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=2,usecols=[1])
hug.us_std    = np.zeros(len(hug.us))
hug.us_units  = ['km/s']
if 0:
    hug.up = (hug.us-res[1])/res[0] # km/s
    ind = np.where(hug.us < ustrans1)
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-ress[1])/ress[0] # km/s
    ind = np.where((hug.us > ustrans1) & (hug.us < ustrans2))
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-rest[1])/rest[0] # km/s
    hug.P = hug.rho_ref*hug.up*hug.us # GPa
    hug.rho = 1/((1/hug.rho_ref[0])*(1-hug.up/hug.us)) # g/cm3
    hug.P_units = ['GPa']
    hug.up_units = ['km/s']
    hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',feng_reference]
hug.refidlist = [0,feng_refid]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','GEKKO']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

C43943 = deepcopy(hug) # keep a copy of the liquid Us-up calibration
shockdat_list.append(hug)


################################################################
# HICKS DECAYING SHOCK TEMPERATURE DATA WITH REFLECTIVITY-TEMPERATURE CORRECTION FROM CELLIERS ET AL. PRL 2010
print('Hicks, D. G., Boehly, T. R., Eggert, J. H., Miller, J. E., Celliers, P. M., & Collins, G. W. (2006). Dissociation of liquid silica at high pressures and temperatures. Physical Review Letters, 97(2), 025502. https://doi.org/10.1103/PhysRevLett.97.025502')

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.doi = ['10.1103/PhysRevLett.97.025502']
newref.reference = ['Hicks, D. G., Boehly, T. R., Eggert, J. H., Miller, J. E., Celliers, P. M., & Collins, G. W. (2006). Dissociation of liquid silica at high pressures and temperatures. Physical Review Letters, 97(2), 025502. https://doi.org/10.1103/PhysRevLett.97.025502']
bibstr = "@article{Hicks2006, \
    title={Dissociation of liquid silica at high pressures and temperatures}, \
    volume={97}, \
    DOI={10.1103/PhysRevLett.97.025502}, \
    number={2}, \
    journal={Physical Review Letters}, \
    author={Hicks, D. G. and Boehly, T. R. and Eggert, J. H. and Miller, J. E. and Celliers, P. M. and Collins, G. W.}, \
    year={2006},\
    pages={025502} \
    }"
newref.bibtex = [bibstr]
newref.bibkey = ['Hicks2006']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
#
# 

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

### FUSED SILICA DECAY
# use HPeos_Shock data structure
hug_file   = silicadatadir+'Hicks-2006-Data-Summary2.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref  = [2.204] # g/cm3
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[0],max_rows=26)
hug.us_std    = hug.us*.01
hug.us_units  = ['km/s']
hug.reflectivity = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[3],max_rows=26) # there is more reflectivity data than temperature data; not loaded here
hug.reflectivity_std = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[4],max_rows=26)
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[7],max_rows=26) # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=6,usecols=[8],max_rows=26) # K
hug.T_units  = ['K']

if 0: # update in hugoniot-calcs.py
    hug.up = (hug.us-res[1])/res[0] # km/s
    ind = np.where(hug.us < ustrans1)
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-ress[1])/ress[0] # km/s
    ind = np.where((hug.us > ustrans1) & (hug.us < ustrans2))
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-rest[1])/rest[0] # km/s
    hug.P = hug.rho_ref*hug.up*hug.us # GPa
    hug.rho = 1/((1/hug.rho_ref[0])*(1-hug.up/hug.us)) # g/cm3
    hug.P_units = ['GPa']
    hug.up_units = ['km/s']
    hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.remlist = ['Recalibrated using updated reflectivities from Celliers et al. 2010 doi:10.1103/PhysRevLett.104.184503']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

H06_fs = deepcopy(hug) 
shockdat_list.append(hug)


### QUARTZ DECAY
# use HPeos_Shock data structure
hug_file   = silicadatadir+'Hicks-2006-Data-Summary2.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [174]
hug.material = ['Quartz']
hug.formula  = ['SiO2']
hug.rho_ref  = [2.650] # g/cm3
hug.orientation = ['z']
hug.phase_ref = ['single crystal']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[0])
hug.us_std    = hug.us*.01
hug.us_units  = ['km/s']
hug.reflectivity = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[1])
hug.reflectivity_std = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[2])
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[5]) # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=5,usecols=[6]) # K
hug.T_units  = ['K']

if 0: # update in hugoniot-calcs.py
    hug.up = (hug.us-res[1])/res[0] # km/s
    ind = np.where(hug.us < ustrans1)
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-ress[1])/ress[0] # km/s
    ind = np.where((hug.us > ustrans1) & (hug.us < ustrans2))
    if len(ind)>0:
        hug.up[ind] = (hug.us[ind]-rest[1])/rest[0] # km/s
    hug.P = hug.rho_ref*hug.up*hug.us # GPa
    hug.rho = 1/((1/hug.rho_ref[0])*(1-hug.up/hug.us)) # g/cm3
    hug.P_units = ['GPa']
    hug.up_units = ['km/s']
    hug.rho_units = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.remlist = ['Recalibrated using updated reflectivities from Celliers et al. 2010 doi:10.1103/PhysRevLett.104.184503']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

H06_qtz = deepcopy(hug) 
shockdat_list.append(hug)

################################################################

################################################################
# MILLOT DECAYING SHOCK TEMPERATURE DATA STISHOVITE
print('Millot, M., Dubrovinskaia, N., Černok, A., Blaha, S., Dubrovinsky, L., Braun, D. G., et al. (2015). Shock compression of stishovite and melting of silica at planetary interior conditions. Science, 347(6220), 418–420. https://doi.org/10.1126/science.1261507')

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Millot, M., Dubrovinskaia, N., Černok, A., Blaha, S., Dubrovinsky, L., Braun, D. G., et al. (2015). Shock compression of stishovite and melting of silica at planetary interior conditions. Science, 347(6220), 418–420. https://doi.org/10.1126/science.1261507']
bibstr = "@article{Millot2015, \
    title={Shock compression of stishovite and melting of silica at planetary interior conditions}, \
    volume={347}, \
    DOI={10.1126/science.1261507}, \
    number={6220}, \
    journal={Science}, \
    author={Millot, M. and Dubrovinskaia, N. and Černok, A. and Blaha, S. and Dubrovinsky, L. and Braun, D. G. and Celliers, P. M. and Collins, G. W. and Eggert, J. H. and Jeanloz, R.}, \
    year={2015}, \
    pages={418–420} \
    }"
newref.bibtex = [bibstr]
newref.bibkey = ['Millot2015']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
#
# 

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])
    
### STISHOVITE DECAY
# use HPeos_Shock data structure; not reading in cv - need units
hug_file   = silicadatadir+'Millot2015_Stishovite_Temperature_Cv.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [672]
hug.material = ['Stishovite']
hug.phase_ref = ['single crystal']
hug.orientation = ['c']
hug.formula  = ['SiO2']
hug.rho_ref = [4.30] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[0],max_rows=311)
hug.us_std    = hug.us*.01 # estimated uncertainty in shock vel
hug.us_units  = ['km/s']
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[1],max_rows=311) # K
hug.T_std    = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[2],max_rows=311) # K
hug.T_units  = ['K']
hug.rho        = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[8],max_rows=311)
hug.rho_units  = ['g/cm3']
hug.P        = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[9],max_rows=311)
hug.P_units  = ['GPa']
hug.up        = np.loadtxt(hug_file,delimiter=',',skiprows=13,usecols=[10],max_rows=311)
hug.up_units  = ['km/s']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.remlist = ['Averaged data.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['line VISAR','SOP']
hug.platform = ['laser','OMEGA EP']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

#M15_stv = deepcopy(hug) 
shockdat_list.append(hug)


#------------------ STV HUGONIOT SAME REFERENCE 
### STISHOVITE HUGONIOT - xtl and polyxtl mixed
hug_file   = silicadatadir+'Millot2015_TableS2_converted.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [672]
hug.material = ['Stishovite']
hug.formula  = ['SiO2']
hug.rho_ref = [4.29] # g/cm3
hug.rho_ref_units = ['g/cm3']
hug.phase_ref = ['single xtl and poly-xtl targets']

hug.exp_id    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[0])
hug.target    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[1],dtype='str')
hug.phase_id  = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[2],dtype='str')
hug.lenergy   = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[3])
hug.lenergy_units = ['J']
hug.uflyer    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[4]) # holds quartz shock velocity
hug.uflyer_std = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[5]) 
hug.uflyer_units = ['km/s']
# shock state in sample
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[6])
hug.us_std    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[7])
hug.us_units  = ['km/s']
hug.up        = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[8])
hug.up_std    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[9])
hug.up_units  = ['km/s']
hug.P        = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[10])
hug.P_std    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[11])
hug.P_units  = ['GPa']
hug.rho        = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[12])
hug.rho_std    = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[13])
hug.rho_units  = ['g/cm3']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.remlist = ['target holds planar stack package description','uflyer holds quartz driver shock velocity']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['line VISAR']
hug.platform = ['laser','OMEGA EP']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)

################################################################
print('################################################################')
## Omega data on fused silica T-Us from Millot and Davies - uncalibratable for temperature, unpublished

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Millot and Davies, fused silica decay data, unpublished']
newref.bibkey = ['MDxx FS Omega']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

### ONE SHOT

# Millot and Davies decaying shock data
# use HPeos_Shock data structure: one for each experiment ID
# FUSED SILICA
hug_file   = silicadatadir+'69691_SOPvsUs.txt'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.204] # g/cm3
hug.exp_id        = [69691]
hug.us        = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[0])
hug.us_units  = ['km/s']
hug.T        = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[1])
hug.T_units  = ['counts']

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.remlist = ['Uncalibrated SOP temperature in counts. Used to determine shock velocity range of temperature inversion.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','SOP']
hug.platform = ['laser','Omega EP']
hug.method = ['experimental']
hug.quality = [0]*len(hug.us)
hug.groups = ['solid compounds and minerals']

M69691fs = deepcopy(hug) # keep a copy of the liquid Us-up calibration
shockdat_list.append(hug)

################################################################

## Root et al. 2019 Fused Silica Hugoniot


## Ocampo Fused Silica Hugoniot
# Ocampo, L. K., Winey, J. M., & Duffy, T. S. (2025). Sound Velocity and Grüneisen Parameter in Shock-Melted Silica at Deep Earth Pressures. Journal of Geophysical Research: Solid Earth, 130, e2025JB031366.

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Ocampo, L. K., Winey, J. M., & Duffy, T. S. (2025). Sound Velocity and Grüneisen Parameter in Shock-Melted Silica at Deep Earth Pressures. Journal of Geophysical Research: Solid Earth, 130, e2025JB031366.']
newref.doi = ['10.1029/2025JB031366']
bibstr = "@article{Ocampo2025, author = {L. K. Ocampo and J. M. Winey and T. S. Duffy}, doi = {10.1029/2025JB031366}, journal = {Journal of Geophysical Research: Solid Earth}, pages = {e2025JB031366}, title = {Sound velocity and Grüneisen parameter in shock-melted silica at deep Earth pressures}, volume = {130}, year = {2025}}"
newref.bibtex = [bibstr]
newref.bibkey = ['Ocampo2025']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Root, S., Townsend, J.P. and Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. Journal of Applied Physics, 126(16), p.165901.
# Fused Quartz
hug_file   = silicadatadir+'Ocampo2025-Table2.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.2] # g/cm3
hug.exp_id      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[0],dtype=str) 
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[1]) 
hug.us_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[2]) 
hug.us_units = ['km/s']
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[3]) 
hug.up_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[4]) 
hug.up_units = ['km/s']
hug.rho      =  np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[5]) 
hug.rho_std      =  np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[6]) 
hug.rho_units = ['g/cm3']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[7])
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[8])
hug.P_units = ['GPa']
hug.cs      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[9])
hug.cs_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[10])
hug.cs_units = ['km/s']
hug.gamma      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[11])
hug.gamma_std      = np.loadtxt(hug_file,delimiter=',',skiprows=3,usecols=[12])
hug.gamma_units = ['None']
#print(hug.P)
#print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = [' The impactor was either OFE C101 copper or platinum.']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR','PDV']
hug.platform = ['gas gun']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)



################################################################


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Root, S., Townsend, J.P. and Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. Journal of Applied Physics, 126(16), p.165901.']
newref.doi = ['10.1063/1.5126205']
bibstr = "@article{Root2019,\
  doi = {10.1063/1.5126205},\
  url = {https://doi.org/10.1063/1.5126205},\
  year = {2019},\
  month = oct,\
  publisher = {{AIP} Publishing},\
  volume = {126},\
  number = {16},\
  pages = {165901},\
  author = {Seth Root and Joshua P. Townsend and Marcus D. Knudson},\
  title = {Shock compression of fused silica: An impedance matching standard},\
  journal = {Journal of Applied Physics}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Root2019']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Root, S., Townsend, J.P. and Knudson, M.D., 2019. Shock compression of fused silica: An impedance matching standard. Journal of Applied Physics, 126(16), p.165901.
# Fused Quartz
hug_file   = silicadatadir+'Root-2019-fused-silica.csv'
hug          = hp.HPeos_Shock()
hug.matid    = [143]
hug.material = ['Quartz amorphous']
hug.formula  = ['SiO2']
hug.rho_ref = [2.2] # g/cm3
hug.uflyer      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[2]) 
hug.uflyer_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[3]) 
hug.uflyer_units = ['km/s']
hug.up      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[4]) 
hug.up_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[5]) 
hug.up_units = ['km/s']
hug.us      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[6]) 
hug.us_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[7]) 
hug.us_units = ['km/s']
hug.rho      =  np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[8]) 
hug.rho_std      =  np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[9]) 
hug.rho_units = ['g/cm3']
hug.P      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[10])
hug.P_std      = np.loadtxt(hug_file,delimiter=',',skiprows=1,usecols=[11])
hug.P_units = ['GPa']
#print(hug.P)
#print(hug.rho)

hug.collections = ['Silica']
hug.ref = [1]
hug.reflist = ['',newref.reference]
hug.refidlist = [0,newref.refid[0]]
hug.rem = [1]*len(hug.P)
hug.remlist = ['Al flyer']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
hug.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
hug.diagnostics = ['VISAR']
hug.platform = ['Z machine']
hug.method = ['experimental']
hug.quality = [0]*len(hug.P)
hug.groups = ['solid compounds and minerals']

shockdat_list.append(hug)



################################################################



################################################################

if 0:
    ### fused silica Omega data shock velocity vs. counts
    ### From Davies 2023
    omegafsus = np.loadtxt(silicadatadir+'69691_SOPvsUs.txt',usecols=[0],delimiter=',',dtype=float,skiprows=1)
    omegafsct = np.loadtxt(silicadatadir+'69691_SOPvsUs.txt',usecols=[1],delimiter=',',dtype=float,skiprows=1)
    ### Millot 2015 Stishovite decaying shock
    omegastvt = np.loadtxt(silicadatadir+'Stishovite_2015_Temperature_Cv.csv',usecols=[1],delimiter=',',dtype=float,skiprows=20,max_rows=250)
    omegastvt_err = np.loadtxt(silicadatadir+'Stishovite_2015_Temperature_Cv.csv',usecols=[2],delimiter=',',dtype=float,skiprows=20,max_rows=250)
    omegastvp = np.loadtxt(silicadatadir+'Stishovite_2015_Temperature_Cv.csv',usecols=[9],delimiter=',',dtype=float,skiprows=20,max_rows=250)




################################################################
################################################################
################################################################
######### PHASE BOUNDARY INFORMATION BELOW #####################
################################################################
################################################################
################################################################

print('################################################################')
print('################################################################')
print('################################################################')
print(' EOS entries for database eosdat_list ' )

## Melt curve Andrault 2022 corrected

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Andrault, D., Pison, L., Morard, G., Garbarino, G., Mezouar, M., Bouhifd, M.A. and Kawamoto, T., 2022. Comment on: Melting behavior of SiO2 up to 120 GPa (Andrault et al. 2020). Physics and Chemistry of Minerals, 49(2), p.3.']
newref.doi = ['10.1007/s00269-021-01174-2']
bibstr = "@article{Andrault2022,\
  doi = {10.1007/s00269-021-01174-2},\
  url = {https://doi.org/10.1007/s00269-021-01174-2},\
  year = {2022},\
  month = jan,\
  publisher = {Springer Science and Business Media {LLC}},\
  volume = {49},\
  number = {2},\
  author = {Denis Andrault and L. Pison and G. Morard and G. Garbarino and M. Mezouar and M. A. Bouhifd and T. Kawamoto},\
  title = {Comment on: Melting behavior of {SiO}2 up to 120~{GPa} (Andrault et al. 2020)},\
  journal = {Physics and Chemistry of Minerals}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Andrault2022']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Andrault et al. 2022 updated melt curve
# CRYSTALLINE QUARTZ
dat_file   = silicadatadir+'Andrault2022-silica-melt-updated.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
# interpolated/interpreted boundary points
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[2]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[1]) # K
dat.T_std  = np.zeros(len(dat.T))+250 # K error estimate based on Fig 2 grey band
dat.T_units = ['K']
dat.P_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[2]) # GPa
dat.P_solid_units = ['GPa']
dat.T_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[3]) # K
dat.T_solid_units = ['K']
dat.P_liquid      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[2]) # GPa
dat.P_liquid_units = ['GPa']
dat.T_liquid      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[4]) # K
dat.T_liquid_units = ['K']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.diagnostics = ['XRD']
dat.platform = ['laser heated DAC']
dat.remlist.append('Table 1: P and T are interpolated-interpreted boundaries between the solid and liquid XRD points')
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

# first entry in eosdat_list
eosdat_list = [dat]

# Andrault et al. 2022 updated melt curve
# CRYSTALLINE QUARTZ
dat_file   = silicadatadir+'Andrault-Fig2-sawtooth.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
# interpolated/interpreted boundary points
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1]) # K
dat.T_units = ['K']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.diagnostics = ['XRD']
dat.platform = ['laser heated DAC']
dat.remlist.append('P and T are interpolated-interpreted boundaries between the solid and liquid XRD points')
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

# first entry in eosdat_list
sawtooth = deepcopy(dat)

# Andrault et al. 2022 updated melt curve Fig 2 digitized
# CRYSTALLINE QUARTZ
dat_file   = silicadatadir+'Andrault-2022-liquid-solid-points.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
# interpolated/interpreted boundary points
dat.P_liquid      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0],max_rows=18) # GPa
dat.P_liquid_units = ['GPa']
dat.T_liquid      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1],max_rows=18) # K
dat.T_liquid_units = ['K']
dat.P_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=20,usecols=[0],max_rows=20) # GPa
dat.P_solid_units = ['GPa']
dat.T_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=20,usecols=[1],max_rows=20) # K
dat.T_solid_units = ['K']

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.diagnostics = ['XRD']
dat.platform = ['laser heated DAC']
dat.remlist.append('Solid and liquid points digitized from Fig. 2; pressures not tabulated for red and blue points')
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

# first entry in eosdat_list
eosdat_list.append(dat)

################################################################
print('################################################################')

## Melt curve Gonzales et al. Sci Rep 2016


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['González-Cataldo, F., Davis, S. and Gutiérrez, G., 2016. Melting curve of SiO2 at multimegabar pressures: implications for gas giants and super-Earths. Scientific Reports, 6(1), p.26537.']
newref.doi = ['10.1038/srep26537']
bibstr = "@article{GonzalezCataldo2016,\
  doi = {10.1038/srep26537},\
  url = {https://doi.org/10.1038/srep26537},\
  year = {2016},\
  month = may,\
  publisher = {Springer Science and Business Media {LLC}},\
  volume = {6},\
  number = {1},\
  author = {Felipe Gonz{\'{a}}lez-Cataldo and Sergio Davis and Gonzalo Guti{\'{e}}rrez},\
  title = {Melting curve of {SiO}2 at multimegabar pressures: implications for gas giants and super-Earths},\
  journal = {Scientific Reports}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['GonzalezCataldo2016']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Gonzales melt curve 2016
# CRYSTALLINE QUARTZ: after alpha-PbO2: calculates pyrite, Cotunnite and Fe2P melting
dat_file   = silicadatadir+'Gonzalez_SiO2_melting_curve-SciRep2016.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[4]) # K
dat.T_std      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[5]) # K
dat.T_units = ['K']
dat.rho_l      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[8]) 
dat.rho_l_units = ['g/cm3']
print(dat.P)
print(dat.T)
print(dat.rho_l)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.platform = ['ab initio']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

################################################################
print('################################################################')

## Belonoshko GCA 1995 Melt curve stishovite

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Belonoshko, A.B. and Dubrovinsky, L.S., 1995. Molecular dynamics of stishovite melting. Geochimica et Cosmochimica Acta, 59(9), pp.1883-1889.']
newref.doi = ['10.1016/0016-7037(95)00071-7']
bibstr = "@article{Belonoshko1995,\
  doi = {10.1016/0016-7037(95)00071-7},\
  url = {https://doi.org/10.1016/0016-7037(95)00071-7},\
  year = {1995},\
  month = may,\
  publisher = {Elsevier {BV}},\
  volume = {59},\
  number = {9},\
  pages = {1883--1889},\
  author = {Anatoly B. Belonoshko and Leonid S. Dubrovinsky},\
  title = {Molecular dynamics of stishovite melting},\
  journal = {Geochimica et Cosmochimica Acta}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Belonoshko1995']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])



# Belonoshko 1995 ab initio melt curve
# CRYSTALLINE QUARTZ
dat_file   = silicadatadir+'Belonoshko-1995-Stishovite-melt.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[0])*0.1 # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[1]) # K
dat.T_units = ['K']
dat.T_std = [150]*len(dat.P) # using the same error as reported for 140 GPa for all points
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.platform = ['ab initio']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


################################################################
print('################################################################')

## Fischer 2018 phase boundary laser-DAC exp Stv-CaCl2
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Fischer, R.A., Campbell, A.J., Chidester, B.A., Reaman, D.M., Thompson, E.C., Pigott, J.S., Prakapenka, V.B. and Smith, J.S., 2018. Equations of state and phase boundary for stishovite and CaCl2-type SiO2. American Mineralogist, 103(5), pp.792-802.']
newref.doi = ['10.2138/am-2018-6267']
bibstr = "@article{Fischer2018,\
  doi = {10.2138/am-2018-6267},\
  url = {https://doi.org/10.2138/am-2018-6267},\
  year = {2018},\
  month = may,\
  publisher = {Mineralogical Society of America},\
  volume = {103},\
  number = {5},\
  pages = {792--802},\
  author = {Rebecca A. Fischer and Andrew J. Campbell and Bethany A. Chidester and Daniel M. Reaman and Elizabeth C. Thompson and Jeffrey S. Pigott and Vitali B. Prakapenka and Jesse S. Smith},\
  title = {Equations of state and phase boundary for stishovite and {CaCl}2-type {SiO}2},\
  journal = {American Mineralogist}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Fischer2018']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Fischer et al. 2018 AM stv-CaCl2 phase boundary
dat_file   = silicadatadir+'Fischer-2018-trans-points.txt'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phasename_a = ['stishovite']  # short name of phase a
dat.phasename_b = ['CaCl2-type SiO2']  # short name of phase b
dat.structure_a = ['rutile']  # crystal structure
dat.structure_b = ['CaCl2']  # crystal structure
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=5,usecols=[0],max_rows=4) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=5,usecols=[1],max_rows=4) # K
dat.T_units = ['K']
dat.P_a      = np.loadtxt(dat_file,delimiter=',',skiprows=9,usecols=[0],max_rows=6) # GPa
dat.P_a_units = ['GPa']
dat.T_a      = np.loadtxt(dat_file,delimiter=',',skiprows=9,usecols=[1],max_rows=6) # K
dat.T_a_units = ['K']
dat.P_b      = np.loadtxt(dat_file,delimiter=',',skiprows=15,usecols=[0],max_rows=4) # GPa
dat.P_b_units = ['GPa']
dat.T_b      = np.loadtxt(dat_file,delimiter=',',skiprows=15,usecols=[1],max_rows=4) # K
dat.T_b_units = ['K']

print(dat.P, dat.P_a, dat.P_b)
print(dat.T, dat.T_a, dat.T_b)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
dat.remlist.append('Boundary fit T=64.6(49)*P-2830(350) in GPa and K')
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.platform = ['laser-heated DAC']
dat.method = ['experimental']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)



################################################################
print('################################################################')

## Tsuchiya et al. 2004 GRL Phase boundaries ab initio


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Tsuchiya, T., Caracas, R. and Tsuchiya, J., 2004. First principles determination of the phase boundaries of high‐pressure polymorphs of silica. Geophysical Research Letters, 31(11).']
newref.doi = ['10.1029/2004GL019649']
bibstr = "@article{Tsuchiya2004,\
  doi = {10.1029/2004gl019649},\
  url = {https://doi.org/10.1029/2004gl019649},\
  year = {2004},\
  month = jun,\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {31},\
  number = {11},\
  pages = {304003 (8pp)},\
  author = {Taku Tsuchiya and Razvan Caracas and Jun Tsuchiya},\
  title = {First principles determination of the phase boundaries of high-pressure polymorphs of silica},\
  journal = {Geophysical Research Letters}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Tsuchiya2004']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Tsuchiya, T., Caracas, R. and Tsuchiya, J., 2004. First principles determination of the phase boundaries of high‐pressure polymorphs of silica.
dat_file   = silicadatadir+'Tsuchiya-2004-phase-boundaries.txt'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phasename_a = ['stishovite']  # short name of phase a
dat.phasename_b = ['CaCl2-SiO2']  # short name of phase b
dat.structure_a = ['rutile']  # crystal structure
dat.structure_b = ['CaCl2']  # crystal structure
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=5,usecols=[0],max_rows=14) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=5,usecols=[1],max_rows=14) # K
dat.T_units = ['K']

print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
dat.remlist.append('Boundary fit P = 56.1 + 0.00585 T (K) GPa.')
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


# Tsuchiya, T., Caracas, R. and Tsuchiya, J., 2004. First principles determination of the phase boundaries of high‐pressure polymorphs of silica.
dat_file   = silicadatadir+'Tsuchiya-2004-phase-boundaries.txt'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phasename_a = ['CaCl2-type SiO2']  # short name of phase a
dat.phasename_b = ['alpha-PbO2-type SiO2']  # short name of phase b
dat.structure_a = ['CaCl2']  # crystal structure
dat.structure_b = ['alpha-PbO2']  # crystal structure
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=19,usecols=[0],max_rows=14) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=19,usecols=[1],max_rows=14) # K
dat.T_units = ['K']

print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
dat.remlist.append('Boundary fit P = 56.1 + 0.00585 T (K) GPa.')
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

################################################################
print('################################################################')

## Oganov et al. 2005 phase boundaries


# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Oganov, A.R., Gillan, M.J. and Price, G.D., 2005. Structural stability of silica at high pressures and temperatures. Physical Review B, 71(6), p.064104.']
newref.doi = ['10.1103/PhysRevB.71.064104']
bibstr = "@article{Oganov2005,\
  doi = {10.1103/physrevb.71.064104},\
  url = {https://doi.org/10.1103/physrevb.71.064104},\
  year = {2005},\
  month = feb,\
  publisher = {American Physical Society ({APS})},\
  volume = {71},\
  number = {6},\
  author = {Artem R. Oganov and Michael J. Gillan and G. David Price},\
  title = {Structural stability of silica at high pressures and temperatures},\
  journal = {Physical Review B}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Oganov2005']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Oganov, A.R., Gillan, M.J. and Price, G.D., 2005. Structural stability of silica at high pressures and temperatures. Physical Review B, 71(6), p.064104.
dat_file   = silicadatadir+'Oganov-2005-silica-phase-boundaries.txt'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phasename_a = ['stishovite']  # short name of phase a
dat.phasename_b = ['CaCl2-SiO2']  # short name of phase b
dat.structure_a = ['rutile']  # crystal structure
dat.structure_b = ['CaCl2']  # crystal structure
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=5,usecols=[0],max_rows=10) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=5,usecols=[1],max_rows=10) # K
dat.T_units = ['K']

print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

# Oganov, A.R., Gillan, M.J. and Price, G.D., 2005. Structural stability of silica at high pressures and temperatures. Physical Review B, 71(6), p.064104.
dat_file   = silicadatadir+'Oganov-2005-silica-phase-boundaries.txt'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phasename_a = ['CaCl2-type SiO2']  # short name of phase a
dat.phasename_b = ['alpha-PbO2-type SiO2']  # short name of phase b
dat.structure_a = ['CaCl2']  # crystal structure
dat.structure_b = ['alpha-PbO2']  # crystal structure
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=15,usecols=[0],max_rows=9) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=15,usecols=[1],max_rows=9) # K
dat.T_units = ['K']

print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

# Oganov, A.R., Gillan, M.J. and Price, G.D., 2005. Structural stability of silica at high pressures and temperatures. Physical Review B, 71(6), p.064104.
dat_file   = silicadatadir+'Oganov-2005-silica-phase-boundaries.txt'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phasename_a = ['alpha-PbO2-type SiO2']  # short name of phase a
dat.phasename_b = ['pyrite-type SiO2']  # short name of phase b
dat.structure_a = ['alpha-PbO2']  # crystal structure
dat.structure_b = ['pyrite']  # crystal structure
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=24,usecols=[0],max_rows=10) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=24,usecols=[1],max_rows=10) # K
dat.T_units = ['K']

print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)



################################################################
print('################################################################')


## Usui 2010 phase boundaries
#Says stv->CaCl2 is 2nd order phase trans with no dV
# Says CaCl2 -> alpha PbO2 is 1st order 
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Usui, Y. and Tsuchiya, T., 2010. Ab initio two-phase molecular dynamics on the melting curve of SiO2. Journal of Earth Science, 21(5), pp.801-810.']
newref.doi = ['10.1007/s12583-010-0126-9']
bibstr = "@article{Usui2010,\
  doi = {10.1007/s12583-010-0126-9},\
  url = {https://doi.org/10.1007/s12583-010-0126-9},\
  year = {2010},\
  month = oct,\
  publisher = {Springer Science and Business Media {LLC}},\
  volume = {21},\
  number = {5},\
  pages = {801--810},\
  author = {Yusuke Usui and Taku Tsuchiya},\
  title = {Ab initio two-phase molecular dynamics on the melting curve of {SiO}2},\
  journal = {Journal of Earth Science}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Usui2010']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Usui & Tsuchiya 2010 phase boundaries
dat_file   = silicadatadir+'Usui-2010-Table1.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

Pall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0]) # GPa
Tall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[4]) # K
rhoall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[6]) # g/cm3
phaseall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[8],dtype='str') # Liquid or Solid
#remarkall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[8]) # joint plane

ind = np.where(phaseall == 'Liquid')
dat.P_liquid      = Pall[ind]
dat.P_liquid_units = ['GPa']
dat.T_liquid      = Tall[ind]
dat.T_liquid_units      = ['K']
dat.rho_liquid      = rhoall[ind]
dat.rho_liquid_units      = ['g/cm3']
print(dat.P_liquid)
print(dat.T_liquid)
print(dat.rho_liquid)
ind = np.where(phaseall == 'Solid')
dat.P_solid      = Pall[ind]
dat.P_solid_units      = ['GPa']
dat.T_solid      = Tall[ind]
dat.T_solid_units = ['K']
dat.rho_solid      = rhoall[ind]
dat.rho_solid_units = ['g/cm3']
print(dat.P_solid)
print(dat.T_solid)
print(dat.rho_solid)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

#print(phaseall)


###############################################################
### Usui Melt Curve digitized from Geng & Mohn Figure 1

dat_file   = silicadatadir+'Geng2024-original-data/Usui-2010-mc-digitized.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.P = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0]) # GPa
dat.T = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1]) # K
dat.P_units = ['GPa']
dat.T_units      = ['K']

dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

u10mc = deepcopy(dat) # not in eosdat_list at this time 10/10/25 - not useful - needs to be redigitized
###############################################################


################################################################
print('################################################################')

## Gasparik Phase Diagram for low-pressures, 2nd Edition book 2014

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Gasparik, T., 2003. Phase diagrams for geoscientists (p. 462). Berlin: Springer.']
newref.doi = ['10.1007/978-1-4614-5776-3']
newref.isbn = ['978-1-4614-5775-6']
bibstr = "@book{Gasparik2014,\
  title={Phase diagrams for geoscientists},\
  author={Gasparik, Tibor},\
  doi = {10.1007/978-1-4614-5776-3},\
  url = {https://doi.org/10.1007/978-1-4614-5776-3},\
  edition = {2nd},\
  year={2014},\
  publisher={Springer}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Gasparik2014']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])



# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-melt-curve.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units      = ['K']
#print(dat.P)
#print(dat.T)

# modify the stv melt curve to better match lab data
ind = np.where(dat.T > 2700)[0]
#dat.T[ind]=5000
dat.T[ind]=dat.T[ind]-(dat.T[ind]-2700)*.5

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)
# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-low-high-qtz.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['Low Quartz']  # short name of phase a
dat.phasename_b = ['High Quartz']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
#print(dat.P)
#print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)
# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-crs-trd-tconst.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['Tridymite']  # short name of phase a
dat.phasename_b = ['Cristobalite']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
#print(dat.P)
#print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-coesite-stv.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['Coesite']  # short name of phase a
dat.phasename_b = ['Stishovite']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
ind = np.where(dat.T < 2900)[0]
dat.P = dat.P[ind]
dat.T = dat.T[ind]
#print(dat.P)
#print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-trd-hqtz.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['Tridymite']  # short name of phase a
dat.phasename_b = ['High Quartz']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-cr-hqtz.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['Crystobalite']  # short name of phase a
dat.phasename_b = ['High Quartz']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-lqtz-coesite.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['Low Quartz']  # short name of phase a
dat.phasename_b = ['Coesite']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
#print(dat.P)
#print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


# Gasparik
dat_file   = silicadatadir+'Gasparik-Silica-Phase-Diagram-hqtz-coesite.csv'
dat         = hp.HPeos_PhaseBoundary()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phasename_a = ['High Quartz']  # short name of phase a
dat.phasename_b = ['Coesite']  # short name of phase b

dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1])+273.15 # C to K
dat.T_units  = ['K']
#print(dat.P)
#print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['fitted']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)




################################################################
print('################################################################')
## Shen and Lazor Stv melt curve DAC

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Shen, G. and Lazor, P., 1995. Measurement of melting temperatures of some minerals under lower mantle pressures. Journal of Geophysical Research: Solid Earth, 100(B9), pp.17699-17713.']
newref.doi = ['10.1029/95jb01864']
bibstr = "@article{Shen1995,\
  doi = {10.1029/95jb01864},\
  url = {https://doi.org/10.1029/95jb01864},\
  year = {1995},\
  month = sep,\
  publisher = {American Geophysical Union ({AGU})},\
  volume = {100},\
  number = {B9},\
  pages = {17699--17713},\
  author = {Guoyin Shen and Peter Lazor},\
  title = {Measurement of melting temperatures of some minerals under lower mantle pressures},\
  journal = {Journal of Geophysical Research: Solid Earth}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Shen1995']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
list(newref)

# Shen & Lazer stv melt curve
dat_file   = silicadatadir+'Shen-1995-stv-melt-curve.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.phase_solid = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[0],dtype='str')
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[2])
dat.P_std      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[3])
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[4])
dat.T_std      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[5])
dat.T_units      = ['K']
print(dat.P)
print(dat.T)
print(dat.phase_solid)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.platform = ['laser-heated DAC']
dat.method = ['experimental']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

################################################################

## Luo et al. EPSL 2002 Melt curve stishovite

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Luo, S.N., Çaǧin, T., Strachan, A., Goddard III, W.A. and Ahrens, T.J., 2002. Molecular dynamics modeling of stishovite. Earth and Planetary Science Letters, 202(1), pp.147-157.']
newref.doi = ['10.1016/S0012-821X(02)00749-5']
bibstr = "@article{Luo2002,\
  title={Molecular dynamics modeling of stishovite},\
  author={Luo, Sheng-Nian and Çaǧin, Tahir and Strachan, Alejandro and Goddard III, William A and Ahrens, Thomas J},\
  journal={Earth and Planetary Science Letters},\
  volume={202},\
  number={1},\
  pages={147--157},\
  year={2002},\
  doi={10.1016/S0012-821X(02)00749-5},\
  publisher={Elsevier}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Luo2002']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])



# Luo et al. 2002 stv melt curve
# CRYSTALLINE QUARTZ
dat_file   = silicadatadir+'Luo-2002-stv-melt.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1]) # K
dat.T_units = ['K']
dat.T_std = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2])
dat.T_l      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[3]) # K
print(dat.P)
print(dat.T)
print(dat.T_std)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
dat.remlist = ['T is two-phase model; T_l is single-phase model']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.platform = ['ab initio']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

################################################################


## Geng et al. 2024 melt curve
# Says CaCl2 -> alpha PbO2 is 1st order 
# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Geng, M., and Mohn, C. E. (2024). Ab initio constraints on the melting of silica at high pressures up to 500 GPa. Physical Review B, 109(2), 024106.'] 
newref.doi = ['10.1103/PhysRevB.109.024106']
bibstr = "@article{Geng2024,\
    title={Ab initio constraints on the melting of silica at high pressures up to 500 GPa},\
    doi={10.1103/PhysRevB.109.024106},\
    journal={Physical Review B},\
    volume={109},\
    number={2},\
    url={https://doi.org/10.1103/PhysRevB.109.024106},\
    publisher={American Physical Society},\
    author={Geng, Ming and Mohn, Chris E.},\
    year={2024},\
    month=jan,\
    pages={024106}\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Geng2024']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Geng and Mohn melt curve
dat_file   = silicadatadir+'Geng-2024-data.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

Pall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[2]) # GPa
Tall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1]) # K
rhoall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[3]) # g/cm3
phaseall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[5],dtype='str') # Liquid or Solid
#remarkall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[8]) # joint plane

print(phaseall)
ind = np.where(phaseall == 'L')
print('Geng liquid')
print(phaseall[ind])
dat.P_liquid      = Pall[ind]
dat.P_liquid_units = ['GPa']
dat.T_liquid      = Tall[ind]
dat.T_liquid_units      = ['K']
dat.rho_liquid      = rhoall[ind]
dat.rho_liquid_units      = ['g/cm3']
print(dat.P_liquid)
print(dat.T_liquid)
print(dat.rho_liquid)
ind = np.where(phaseall != 'L')
print('Geng solid')
print(phaseall[ind])
dat.P_solid      = Pall[ind]
dat.P_solid_units      = ['GPa']
dat.T_solid      = Tall[ind]
dat.T_solid_units = ['K']
dat.rho_solid      = rhoall[ind]
dat.rho_solid_units = ['g/cm3']
print(dat.P_solid)
print(dat.T_solid)
print(dat.rho_solid)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

#print(phaseall)


################################################################

if 0:
    ## OLD USING DIGITIZED POINTS
    ## Geng et al. 2024 melt curve Figure 2 - liquid densities
    # reference already in reference list
    
    # Geng and Mohn melt curve liquid density Figure 2A
    dat_file   = silicadatadir+'Geng2024-melt-density-on-melt-curve.txt'
    dat         = hp.HPeos_MeltCurve()
    dat.matid    = [174]
    dat.material = ['Quartz']
    dat.formula  = ['SiO2']
    dat.rho_ref = [2.650] # g/cm3
    
    Pall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0]) # GPa
    rhoall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1]) # g/cm3
    dat.P_liquid      = Pall
    dat.P_liquid_units = ['GPa']
    dat.rho_liquid      = rhoall
    dat.rho_liquid_units      = ['g/cm3']
    print(dat.P_liquid)
    print(dat.rho_liquid)
    
    dat.ref = [1]
    dat.reflist = ['',newref.reference]
    dat.refidlist = [0,newref.refid[0]]
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
    dat.method = ['computational']
    dat.groups = ['solid compounds and minerals']
    
    eosdat_list.append(dat)

## NEW USING AUTHOR PROVIDED POINTS 7/31/2025
## Geng et al. 2024 melt curve Figure 2 - liquid densities
# reference already in reference list

# Geng and Mohn melt curve liquid density Figure 2A
dat_file   = silicadatadir+'Geng2024-original-data/Fig2-melt-density.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

Pall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[0]) # GPa
Tall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[1]) # K
rhoall = np.loadtxt(dat_file,delimiter=',',skiprows=1,usecols=[2]) # g/cm3
dat.P_liquid      = Pall
dat.P_liquid_units = ['GPa']
dat.T_liquid      = Tall
dat.T_liquid_units = ['K']
dat.rho_liquid      = rhoall
dat.rho_liquid_units      = ['g/cm3']
print(dat.P_liquid)
print(dat.rho_liquid)

dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)



################################################################

## Geng et al. 2024 melt curve Figure 2 - CN on melt curve
# reference already in reference list

# Geng and Mohn melt curve CN of liquid Figure 2A - original data from authors
dat_file   = silicadatadir+'Geng2024-original-data/Fig2-melt-CN.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

Pall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0]) # GPa
CNall = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1]) # CN@
dat.P_liquid      = Pall
dat.P_liquid_units = ['GPa']
dat.CN_liquid      = CNall
print(dat.P_liquid)
print(dat.CN_liquid)

dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

###############################################################
### Geng & Mohn 2024 Melt Curve Figure 1 caption: text is correct; Figure 1 plot is not
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.P      = np.arange(500)+50 # GPa PREFERRED FIT BY GENG Email 8/5/2024 2-phase fit
dat.P_units = ['GPa']
dat.T      = 6.018*np.power((1+(dat.P+202.338)/4.736e-04),0.5126)
dat.T_std      = np.zeros(len(dat.T))+200 # estimated 1 sigma from author emailed figure
dat.T_units      = ['K']
dat.rho = np.zeros(len(dat.P)) # will be filled later

dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

gm24mc = deepcopy(dat) # not in eosdat_list at this time 10/10/25

### Geng & Mohn 2024 Melt Curve Figure 1 caption: text is correct; Figure 1 plot is not
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3

dat.P      = np.arange(500)+50 # GPa 
dat.P_units = ['GPa']
dat.T      = 6.426*np.power((1+(dat.P+256.011)/6.474e-04),0.5141) # z-method fit (not preferred)
dat.T_units      = ['K']
dat.rho = np.zeros(len(dat.P)) # will be filled later

dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.method = ['computational']
dat.groups = ['solid compounds and minerals']

gm24mcz = deepcopy(dat) # not in eosdat_list at this time 10/10/25
###############################################################

################################################################

#=============================================
## Melt curve Andrault 2020 original

# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Andrault, D., Morard, G., Garbarino, G., Mezouar, M., Bouhifd, M. A., & Kawamoto, T. (2020). Melting behavior of SiO2 up to 120 GPa. Physics and Chemistry of Minerals, 47(2), 10. https://doi.org/10.1007/s00269-019-01077-3']
newref.doi = ['10.1007/s00269-019-01077-3']
bibstr = "@article{Andrault2020,\
	title = {Melting behavior of {SiO$_2$} up to 120 {GPa}},\
	volume = {47},\
	issn = {1432-2021},\
	url = {https://doi.org/10.1007/s00269-019-01077-3},\
	doi = {10.1007/s00269-019-01077-3},\
	number = {2},\
	journal = {Physics and Chemistry of Minerals},\
	author = {Andrault, D. and Morard, G. and Garbarino, G. and Mezouar, M. and Bouhifd, M. A. and Kawamoto, T.},\
	month = jan,\
	year = {2020},\
	keywords = {Lower mantle pressures, Melt structure, Melting diagram, SiO2},\
	pages = {10},\
}"
newref.bibtex = [bibstr]
newref.bibkey = ['Andrault2020']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])


# Andrault et al. 2020 original melt curve
# CRYSTALLINE QUARTZ
dat_file   = silicadatadir+'Andrault-2020-Table-1-only.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
# interpolated/interpreted boundary points
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1]) # K
dat.T_units = ['K']
dat.P_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2]) # GPa
dat.P_solid_units = ['GPa']
dat.T_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[3]) # K
dat.T_solid_units = ['K']
dat.P_liquid      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2]) # GPa
dat.P_liquid_units = ['GPa']
dat.T_liquid      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[4]) # K
dat.T_liquid_units = ['K']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.diagnostics = ['XRD']
dat.platform = ['laser heated DAC']
dat.remlist.append('Table 1: P and T are interpolated-interpreted boundaries between the solid and liquid XRD points')
dat.quality = [0]
dat.groups = ['solid compounds and minerals']


eosdat_list.append(dat)

################################################################

#=============================================
## Melt curve Gonzalez et al. 2016 addition pyrite points


## Melt curve Gonzales et al. Sci Rep 2016
# look up reference from above
xx = bibkeyid(ref_list,'GonzalezCataldo2016')
print(ref_list[xx].reference[0])

# additional pyrite points
pyritedat = np.loadtxt(silicadatadir+'Pyrite-type_SiO2_melting_curve_FGonzalez-SciRep2016.txt',delimiter=',',skiprows=2)

dat_file   = silicadatadir+'Pyrite-type_SiO2_melting_curve_FGonzalez-SciRep2016.txt'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[0]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[1]) # K
dat.T_units = ['K']
dat.U      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[2]) 
dat.U_units = ['eV/SiO2']
dat.rho_l      = np.loadtxt(dat_file,delimiter=',',skiprows=2,usecols=[3]) 
dat.rho_l_units = ['g/cm3']
print(dat.P)
print(dat.T)
print(dat.rho_l)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',ref_list[xx].reference]
dat.refidlist = [0,ref_list[xx].refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
dat.platform = ['ab initio','PBE','Z method']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

################################################################

#=============================================
## Park et al. 2025 sft-pyr-melt curve



# add reference to reference list

maxrefid=0
for item in ref_list:
    if int(item.refid[0]) > maxrefid:
        maxrefid = int(item.refid[0])
        
print('Previous max refid # = ',maxrefid)
newrefid = maxrefid+1
print('New refid # = ',newrefid)

newref = hp.HPeos_Reference()
newref.refid = [newrefid]
newref.reference = ['Park, D., Deng, X., & Deng, J. (2025, October 15). Melting phase relation of seifertite and pyrite-type SiO2 determined by machine learning potentials. arXiv. https://doi.org/10.48550/arXiv.2510.13119']
newref.doi = ['xxx']
bibstr = "@article{Park_Deng_Deng_2025, title={Melting phase relation of seifertite and pyrite-type SiO2 determined by machine learning potentials}, url={http://arxiv.org/abs/2510.13119}, DOI={10.48550/arXiv.2510.13119}, abstractNote={Silica (SiO2) is fundamental to both industrial technology and planetary science, yet the phase relations of its high-pressure polymorphs remain poorly constrained. Here, we develop two machine learning potentials (MLPs) for SiO2 that faithfully represent the SCAN and PBEsol exchange-correlation functionals over a wide temperature (1000-10000 K) and pressure (100-400 GPa) range using deep neural networks. With large-scale two-phase simulations powered by these potentials, we determine the melting curves of seifertite and pyrite-type SiO2 and infer the solid-solid phase boundary between these two phases. The SCAN functional, which captures intermediate-range van der Waals interactions, reproduces structural and thermodynamic properties with high fidelity, predicting melting temperatures 6-10 % higher and a seifertite to pyrite-type transition pressure 22 % higher than the PBEsol. The strongly negative Clapeyron slope (-6.1 MPa/K) of this transition suggests that mantle convection could be highly layered in super-Earth exoplanets, potentially affecting their long-term thermal evolution and habitability.}, note={arXiv:2510.13119 [cond-mat]}, number={arXiv:2510.13119}, publisher={arXiv}, author={Park, Doyoon and Deng, Xin and Deng, Jie}, year={2025}, month=oct}"
newref.bibtex = [bibstr]
newref.bibkey = ['Park2025']
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
newref.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]

ref_list.append(newref)
#pprint(list(newref))
print(newref.reference[0])

dat_file   = silicadatadir+'Park2025-Fig5-pb.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[0]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[1]) # K
dat.T_units = ['K']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry sft-pyr boundary']]
dat.platform = ['ab initio','PBE']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)

#-------------- PBE melt curve --------------------
# PBE solidus, liquidus
dat_file   = silicadatadir+'Park2025-densities.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phase_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[0],dtype=str) 
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[1]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[3]) # K
dat.T_units = ['K']
dat.rho_s      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[5]) 
dat.rho_s_units = ['g/cm3']
dat.rho_l      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[7]) 
dat.rho_l_units = ['g/cm3']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'Park PBE melt curve']]
dat.platform = ['ab initio','PBE']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)


#-------------- SCAN melt curve --------------------
# SCAN solidus, liquidus
dat_file   = silicadatadir+'Park2025-densities.csv'
dat         = hp.HPeos_MeltCurve()
dat.matid    = [174]
dat.material = ['Quartz']
dat.formula  = ['SiO2']
dat.rho_ref = [2.650] # g/cm3
dat.phase_solid      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[0],dtype=str) 
dat.P      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[1]) # GPa
dat.P_units = ['GPa']
dat.T      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[2]) # K
dat.T_units = ['K']
dat.rho_s      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[4]) 
dat.rho_s_units = ['g/cm3']
dat.rho_l      = np.loadtxt(dat_file,delimiter=',',skiprows=3,usecols=[6]) 
dat.rho_l_units = ['g/cm3']
print(dat.P)
print(dat.T)

dat.collection = ['Silica']
dat.ref = [1]
dat.reflist = ['',newref.reference]
dat.refidlist = [0,newref.refid[0]]
DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'Park SCAN melt curve']]
dat.platform = ['ab initio','PBE']
dat.quality = [0]
dat.groups = ['solid compounds and minerals']

eosdat_list.append(dat)



################################################################
################################################################
# moved to define-mc.py
if 0:
    print('################################################################')
    print(' Read in sts fit for melt curve envelope ')
    
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
    print(dat.P)
    print(dat.T_s)
    print(dat.T_l)
    
    dat.collection = ['Silica']
    dat.ref = [1]
    dat.reflist = ['',newref.reference]
    dat.refidlist = [0,newref.refid[0]]
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
    dat.method = ['judgement']
    dat.groups = ['solid compounds and minerals']
    
    stsmc = deepcopy(dat)
    
    
    
    # STS estimated melt curve line for data review draw as dotted above 155 GPa
    dat_file   = bookpath+'scripts/melt-curve.txt'
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
    print(dat.P)
    print(dat.T)
    
    dat.collection = ['Silica']
    dat.ref = [1]
    dat.reflist = ['',newref.reference]
    dat.refidlist = [0,newref.refid[0]]
    DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #
    dat.actions = [['Sarah T. Stewart','create',DateTime_in_ISOFormat,'new entry']]
    dat.method = ['judgement']
    dat.groups = ['solid compounds and minerals']
    
    if 1:
        # interpolate to 0.5 GPa from 14.5 to 251 GPa
        mc_p = np.arange(14.,501.,0.1) # GPa
        mc_t = np.interp(mc_p,dat.P,dat.T)
        dat.P = mc_p
        dat.T = mc_t
    dat.rho =  np.zeros(len(dat.P)) # will be filled later
    stsmcsmooth = deepcopy(dat)
    


################################################################
print('################################################################')
print(' Add some decaying shock data == not in database yet! ')

# MOVE TO hugoniot-calcs.py
if 0:
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


################################################################
print('################################################################')
print('################################################################')
print('################################################################')

# make shorter names for the database list
shockdat = shockdat_list
eosdat = eosdat_list

################################################################
print('################################################################')
print('################################################################')
print('################################################################')

###
# find reference IDs for plots
print('Number of shockdat_list entries', len(shockdat_list))
for iii in range(0,len(shockdat_list)):
    print(iii,shockdat_list[iii].material,shockdat_list[iii].reflist[1])


# find reference IDs for plots
print('Number of EOS entries', len(eosdat_list))
for iii in range(0,len(eosdat_list)):
    print(iii,eosdat_list[iii].reflist)

################################################################
## SHORTCUTS FOR SHOCKDAT_LIST
################################################################
# dictionary labels for plots and indices indices for shockdat_list

L80qtz = {'lab' : 'L80 qtz', 'id' : 0}
L80fs = {'lab' : 'L80 FS', 'id' : 1}
L83qtz = {'lab' : 'L83 qtz', 'id' : 2}
L83fs = {'lab' : 'L83 FS', 'id' : 3}
L83qtzhug = {'lab' : 'L83 qtz', 'id' : 4}
L83fshug = {'lab' : 'L83 FS', 'id' : 5}

B88qtz = {'lab' : 'B88 qtz', 'id' : 6}
B88fs = {'lab' : 'B88 FS', 'id' : 7}
L04qtz = {'lab' : 'L04 qtz', 'id' : 8}
A23qtz = {'lab' : 'A23 qtz', 'id' : 9}
A23fs = {'lab' : 'A23 FS', 'id' : 10}
W62qtz = {'lab' : 'W62 qtz', 'id' : 11}
W62fs = {'lab' : 'W62 FS', 'id' : 12}
L03qtz = {'lab' : 'L03 qtz model', 'id' : 13}
R19fs = {'lab' : 'R19 FS', 'id' : 14}
T18fs = {'lab' : 'T18 FS', 'id' : 15}
T20qtz = {'lab' : 'T20 qtz', 'id' : 16}
S80fs = {'lab' : 'S80 FS', 'id' : 17}
# millot and davies decaying shocks 18-21
MDqtz70811 = {'lab' : '70811', 'id' : 18}
MDqtz70813 = {'lab' : '70813', 'id' : 19}
MDqtz70817 = {'lab' : '70817', 'id' : 20}
MDqtz70820 = {'lab' : '70820', 'id' : 21}
S81qtz = {'lab' : 'S81 qtz', 'id' : 22}
S81fs = {'lab' : 'S81 FS', 'id' : 23}
S82fs = {'lab' : 'S82 FS', 'id' : 24} # Sugiura et al. JAP 1982
B19fs = {'lab' : 'B19 FS', 'id' : 25} # Berryman et al. GRL 2019
K13qtz = {'lab' : 'K13 qtz', 'id' : 26} # Knudson & Desjarlais PRB 2013
# do not use
#B88Rqtz = {'lab' : 'B88R qtz', 'id' : 27} # Boslough 1988 revised shock Pressures STS
#B88Rfs = {'lab' : 'B88R FS', 'id' : 28} # Boslough 1988 revised shock Pressures STS
J79fs = {'lab' : 'J79 FS', 'id' : 29} # Jackson & Ahrens JGR FS hug
J79Rfs = {'lab' : 'J79R FS', 'id' : 30} # Jackson & Ahrens JGR FS hug
B88RRqtz = {'lab' : 'B88R qtz', 'id' : 31} # Boslough 1988 revised shock Pressures STS redo IM
B88RRfs = {'lab' : 'B88R FS', 'id' : 32} # Boslough 1988 revised shock Pressures STS redo IM
F24cst = {'lab' : 'F24 cst', 'id' : 33} # Feng et al. GRL 2024 Coesite Hugoniot
L02cst = {'lab' : 'L02 cst', 'id' : 34} # Luo et al. GRL 2002 Coesite Hugoniot
L02stv = {'lab' : 'L02 stv', 'id' : 35} # Luo et al. GRL 2002 Stishovite Hugoniot
P81cst = {'lab' : 'P81 cst', 'id' : 36} # Podurets Izv. Akad. 1981 Coesite Hugoniot
J68fs = {'lab' : 'J68 FS', 'id' : 37} # Jones et al. 1968 NASA NTRS 19680019159
# Feng et al. GRL 2024 3 decaying Us-T 38, 39, 40
F24cst43939 = {'lab' : '43939', 'id' : 38}
F24cst43941 = {'lab' : '43941', 'id' : 39}
F24cst43943 = {'lab' : '43943', 'id' : 40}
H06fs = {'lab' : 'H06 FS', 'id' : 41} # Hicks et al. 2006; recalibrated with Celliers et al. PRL 2010
H06qtz = {'lab' : 'H06 qtz', 'id' : 42} # Hicks et al. 2006; recalibrated with Celliers et al. PRL 2010
M15stv = {'lab' : 'M15 stv', 'id' : 43} # Millot et al. Science 2015 stishovite decaying shock data
M15stvhug = {'lab' : 'M15 stv', 'id' : 44} # Millot et al. Science 2015 stishovite steady shock points
MDfs69691 = {'lab' : '69691', 'id' : 45} # Millot & Davies unpublished fused silica decaying shock
O25fs = {'lab' : 'O25 fs', 'id' : 46}  # O'Campo et al. JGR 2025 fused silica hug, sound speeds, gamma

# alternatively
# hp.HPeos_ID data class to hold these shortcuts
# e.g., L80qtz.label, L80qtz.i for the index in shockdat_list
iL80qtz = hp.HPeos_ID('L80 qtz', 0)
iL80fs = hp.HPeos_ID('L80 FS', 1)
iL83qtz = hp.HPeos_ID('L83 qtz', 2)
iL83fs = hp.HPeos_ID('L83 FS', 3)
iL83qtzhug = hp.HPeos_ID('L83 qtz', 4)
iL83fshug = hp.HPeos_ID('L83 FS', 5)
iB88qtz = hp.HPeos_ID('B88 qtz', 6)
iB88fs = hp.HPeos_ID('B88 FS', 7)
iL04qtz = hp.HPeos_ID('L04 qtz', 8)
iA23qtz = hp.HPeos_ID('A23 qtz', 9)
iA23fs = hp.HPeos_ID('A23 FS', 10)
iW62qtz = hp.HPeos_ID('W62 qtz', 11)
iW62fs = hp.HPeos_ID('W62 FS', 12)
iL03qtz = hp.HPeos_ID('L03 qtz model', 13)
iR19fs = hp.HPeos_ID('R19 FS', 14)
iT18fs = hp.HPeos_ID('T18 FS', 15)
iT20qtz = hp.HPeos_ID('T20 qtz', 16)
iS80fs = hp.HPeos_ID('S80 FS', 17)
iMDqtz70811 = hp.HPeos_ID('70811', 18)
iMDqtz70813 = hp.HPeos_ID('70813', 19)
iMDqtz70817 = hp.HPeos_ID('70817', 20)
iMDqtz70820 = hp.HPeos_ID('70820', 21)
iS81qtz = hp.HPeos_ID('S81 qtz', 22)
iS81fs = hp.HPeos_ID('S81 FS', 23)
iS82fs = hp.HPeos_ID('S82 FS', 24) # Sugiura et al. JAP 1982
iB19fs = hp.HPeos_ID('B19 FS', 25) # Berryman et al. GRL 2019
iK13qtz = hp.HPeos_ID('K13 qtz', 26) # Knudson & Desjarlais PRB 2013
iJ79fs = hp.HPeos_ID('J79 FS', 29) # Jackson & Ahrens JGR FS hug
iJ79Rfs = hp.HPeos_ID('J79R FS', 30) # Jackson & Ahrens JGR FS hug
iB88RRqtz = hp.HPeos_ID('B88R qtz', 31) # Boslough 1988 revised shock Pressures STS redo IM
iB88RRfs = hp.HPeos_ID('B88R FS', 32) # Boslough 1988 revised shock Pressures STS redo IM
iF24cst = hp.HPeos_ID('F24 cst', 33) # Feng et al. GRL 2024 Coesite Hugoniot
iL02cst = hp.HPeos_ID('L02 cst', 34) # Luo et al. GRL 2002 Coesite Hugoniot
iL02stv = hp.HPeos_ID('L02 stv', 35) # Luo et al. GRL 2002 Stishovite Hugoniot
iP81cst = hp.HPeos_ID('P81 cst', 36) # Podurets Izv. Akad. 1981 Coesite Hugoniot
iJ68fs = hp.HPeos_ID('J68 FS', 37) # Jones et al. 1968 NASA NTRS 19680019159
iF24cst43939 = hp.HPeos_ID('43939', 38)
iF24cst43941 = hp.HPeos_ID('43941', 39)
iF24cst43943 = hp.HPeos_ID('43943', 40)
iH06fs = hp.HPeos_ID('H06 FS', 41) # Hicks et al. 2006; recalibrated with Celliers et al. PRL 2010
iH06qtz = hp.HPeos_ID('H06 qtz', 42) # Hicks et al. 2006; recalibrated with Celliers et al. PRL 2010
iM15stv = hp.HPeos_ID('M15 stv', 43) # Millot et al. Science 2015 stishovite decaying shock data
iM15stvhug = hp.HPeos_ID('M15 stv', 44) # Millot et al. Science 2015 stishovite steady shock
iMDfs69691 = hp.HPeos_ID('69691', 45) # Millot & Davies unpublished fused silica decaying shock
iO25fs = hp.HPeos_ID('O25 FS', 46) # O'Campo et al. JGR 2025 fused silica hug, sound speeds, gamma


################################################################
## SHORTCUTS FOR EOSDAT_LIST
################################################################
# dictionary labels for plots and indices indices for eosdat_list
A22mc = {'lab' : 'A22 melt LDAC-XRD', 'id' : 0} # laser-heated DAC+XRD
A22mcp = {'lab' : 'A22 melt exp plot', 'id' : 1} # laser-heated DAC+XRD
G16mc = {'lab' : 'G16 melt pyrite PBE', 'id' : 2} # PBE DFT
G16mcb = {'lab' : 'G16 melt cotunnite PBE', 'id' : 2} # PBE DFT
G16mcc = {'lab' : 'G16 melt Fe₂P PBE', 'id' : 2} # PBE DFT
B95mc = {'lab' : 'B95 melt stv TTAM', 'id' : 3} # TTAM-IP MD sims 
F18pb = {'lab' : 'F18 stv-CaCl₂ LDAC-XRD', 'id' : 4} # laser-heated DAC+XRD
T04pba = {'lab' : 'T04 stv-CaCl₂ PBE', 'id' : 5} # GGA-PBE DFT
T04pbb = {'lab' : 'T04 CaCl₂-sft PBE', 'id' : 6} # GGA-PBE DFT
O05pba = {'lab' : 'O05 stv-CaCl₂ LDA', 'id' : 7} # LDA DFT
O05pbb = {'lab' : 'O05 CaCl₂-sft LDA', 'id' : 8} # LDA DFT
O05pbc = {'lab' : 'O05 sft-pyrite LDA', 'id' : 9} # LDA DFT
U10mcl = {'lab' : 'U10 liquid LDA', 'id' : 10} # LDA DFT
U10mcs = {'lab' : 'U10 solid LDA', 'id' : 10} # LDA DFT
G14mc = {'lab' : 'G14 melt curve', 'id' : 11} # Gasparik book phase boundaries for low-pressure phases (still some issues with Coesite-Stv)
G14lhqtz = {'lab' : 'G14 low qtz-high qtz phase boundary', 'id' : 12}
G14trcr = {'lab' : 'G14 tridymite-cristobalite phase boundary', 'id' : 13}
G14crstv = {'lab' : 'G14 coesite-stishovite phase boundary', 'id' : 14}
G14trhq = {'lab' : 'G14 tridymite-high qtz phase boundary', 'id' : 15}
G14crhq = {'lab' : 'G14 cristobalite-high qtz phase boundary', 'id' : 16}

G14lqcs = {'lab' : 'G14 low qtz-coesite phase boundary', 'id' : 17}
G14hqcs = {'lab' : 'G14 high qtz-coesite phase boundary', 'id' : 18}

S95mcstv = {'lab' : 'S95 melt LDAC', 'id' : 19} # laser-heated DAC+optical

L02mcstv = {'lab' : 'L02 melt stv CMD', 'id' : 20} # Classical MD; does not agree with other data
G24mcl = {'lab' : 'G24 liquid PBE', 'id' : 21} # GGA-PBE DFT 
G24mcs = {'lab' : 'G24 solid PBE', 'id' : 21} # GGA-PBE DFT # same entry in list different fields
G24mclrho = {'lab' : 'G24 liquidus PBE', 'id' : 22} # GGA-PBE DFT # new entry
G24mclcn = {'lab' : 'G24 CN PBE', 'id' : 23} # GGA-PBE DFT # new entry
A20mc = {'lab' : 'A20 melt LDAC-XRD', 'id' : 24} # laser-heated DAC+XRD
G16mcpyr = {'lab' : 'G16 melt pyrite PBE', 'id' : 25} # PBE DFT z method additional points
P25pb = {'lab' : 'P25 sft-pyrite PBEsol', 'id' : 26} # PBEsol sft-pyr phase boundary
P25mcpbe = {'lab' : 'P25 melt PBEsol', 'id' : 27} # PBEsol 
P25mcscan = {'lab' : 'P25 melt SCAN', 'id' : 28} # PBEsol 

# alternatively
# hp.HPeos_ID data class to hold these shortcuts
# e.g., A22mc.label, A22mc.i for the index in eosdat_list
iA22mc = hp.HPeos_ID( 'A22 melt LDAC-XRD', 0) # laser-heated DAC+XRD
iA22mcp = hp.HPeos_ID( 'A22 melt exp plot', 1) # laser-heated DAC+XRD
iG16mc = hp.HPeos_ID( 'G16 melt pyrite PBE', 2) # PBE DFT
iG16mcb = hp.HPeos_ID( 'G16 melt cotunnite PBE', 2) # PBE DFT
iG16mcc = hp.HPeos_ID( 'G16 melt Fe₂P PBE', 2) # PBE DFT
iB95mc = hp.HPeos_ID( 'B95 melt stv TTAM', 3) # TTAM-IP MD sims 
iF18pb = hp.HPeos_ID( 'F18 stv-CaCl₂ LDAC-XRD', 4) # laser-heated DAC+XRD
iT04pba = hp.HPeos_ID( 'T04 stv-CaCl₂ PBE', 5) # GGA-PBE DFT
iT04pbb = hp.HPeos_ID( 'T04 CaCl₂-sft PBE', 6) # GGA-PBE DFT
iO05pba = hp.HPeos_ID( 'O05 stv-CaCl₂ LDA', 7) # LDA DFT
iO05pbb = hp.HPeos_ID( 'O05 CaCl₂-sft LDA', 8) # LDA DFT
iO05pbc = hp.HPeos_ID( 'O05 sft-pyrite LDA', 9) # LDA DFT
iU10mcl = hp.HPeos_ID( 'U10 liquid LDA', 10) # LDA DFT
iU10mcs = hp.HPeos_ID( 'U10 solid LDA', 10) # LDA DFT
iG14mc = hp.HPeos_ID( 'G14 melt curve', 11) # Gasparik book phase boundaries for low-pressure phases (still some issues with Coesite-Stv)
iG14lhqtz = hp.HPeos_ID( 'G14 low qtz-high qtz phase boundary', 12)
iG14trcr = hp.HPeos_ID( 'G14 tridymite-cristobalite phase boundary', 13)
iG14crstv = hp.HPeos_ID( 'G14 coesite-stishovite phase boundary', 14)
iG14trhq = hp.HPeos_ID( 'G14 tridymite-high qtz phase boundary', 15)
iG14crhq = hp.HPeos_ID( 'G14 cristobalite-high qtz phase boundary', 16)
iG14lqcs = hp.HPeos_ID( 'G14 low qtz-coesite phase boundary', 17)
iG14hqcs = hp.HPeos_ID( 'G14 high qtz-coesite phase boundary', 18)
iS95mcstv = hp.HPeos_ID( 'S95 melt LDAC', 19) # laser-heated DAC+optical
iL02mcstv = hp.HPeos_ID( 'L02 melt stv CMD', 20) # Classical MD; does not agree with other data
iG24mcl = hp.HPeos_ID( 'G24 liquid PBE', 21) # GGA DFT PBE
iG24mcs = hp.HPeos_ID( 'G24 solid PBE', 21) # GGA DFT PBE
iG24mclrho = hp.HPeos_ID( 'G24 liquidus PBE', 22) # GGA DFT PBE
iG24mclcn = hp.HPeos_ID( 'G24 CN PBE', 23) # GGA DFT PBE
iA20mc = hp.HPeos_ID( 'A20 melt LDAC-XRD', 24) # laser-heated DAC+XRD
iG16mcpyr = hp.HPeos_ID( 'G16 melt pyrite PBE', 25) # additional pyrite z method melt curve points not in file loaded above
iP25pb = hp.HPeos_ID( 'P25 PBEsol', 26) # PBEsol sft-pyr phase boundary
iP25mcpbe = hp.HPeos_ID( 'P25 mc PBEsol', 27) # PBEsol sft-pyr phase boundary
iP25mcscan = hp.HPeos_ID( 'P25 mc SCAN', 28) # PBEsol sft-pyr phase boundary

#for iii in range(0,len(ref_list)):
#    print(iii,ref_list[iii].bibkey,ref_list[iii].doi)#ref_list[iii].reference,)

#print('Number of shockdat_list entries', len(shockdat_list))
#for iii in range(0,len(shockdat_list)):
#    print(iii,shockdat_list[iii].material,shockdat_list[iii].refidlist,shockdat_list[iii].quality)

if 0: # move dup
    ### fused silica Omega data shock velocity vs. counts
    ### From Davies 2023
    omegafsus = np.loadtxt(silicadatadir+'69691_SOPvsUs.txt',usecols=[0],delimiter=',',dtype=float,skiprows=1)
    omegafsct = np.loadtxt(silicadatadir+'69691_SOPvsUs.txt',usecols=[1],delimiter=',',dtype=float,skiprows=1)
    
    ### Millot 2015 Stishovite decaying shock
    omegastvt = np.loadtxt(silicadatadir+'Stishovite_2015_Temperature_Cv.csv',usecols=[1],delimiter=',',dtype=float,skiprows=20,max_rows=250)
    omegastvt_err = np.loadtxt(silicadatadir+'Stishovite_2015_Temperature_Cv.csv',usecols=[2],delimiter=',',dtype=float,skiprows=20,max_rows=250)
    omegastvp = np.loadtxt(silicadatadir+'Stishovite_2015_Temperature_Cv.csv',usecols=[9],delimiter=',',dtype=float,skiprows=20,max_rows=250)



