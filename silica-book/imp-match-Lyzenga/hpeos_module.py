## Functions to support High Pressure EOS Database
## include this command in python: import hpeos_module as hp
## S. T. Stewart
## ststewart@ucdavis.edu
#
# Needs: a way to report fits with covariances; this can be a separate structure
# but it needs a customized format? ok for MongoDB -- a text string that provides the function
# parameters list and covariance list; problems can be ordering of parameters
#
__version__='0.2.1' # 7/8/2023 - set independent repo version control
#__version__='0.1.1' # 5/13/2023 adding more data types; meltcurve; phase boundaries; NEED generic EOS/thermo info structure
#__version__='0.1' # 4/27/2023
#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import deepcopy
#
# every entry is a list because empty list is a valid object for regexp
# when different materials have different parameters; empty lists do not throw errors
#
class HPeos_Material:
    """Class database materials and associated database entries."""
    def __init__(self):
        # primary database material information
        self.matid         = []     # Database Material ID number = IHED ID for those in IHED
        self.material      = []     # Database Material name
        self.aliases       = []     # Other names used for the same material
        self.formula       = []     # Database Material chemical formula
        self.phase         = []     # list of short name(s) for the material phase/crystal structure
        self.db_items      = []     # database object IDs for material properties data associated with this material
        self.db_refs       = []     # database object IDs for references associated with this material
        self.groups        = []     # Rusbank material groups
        self.collections   = []     # Collection name for developer (and future user) customized groupings
        self.actions       = []     # append actions: 'create', 'review', 'add', 'remove', 'revise'
                                    # action format: [name, action, date, note]                                
    def __iter__(self):
        '''Needed to convert class object into dictionary.'''
        for key in self.__dict__:
            yield key, getattr(self,key)
#
class HPeos_Shock:
    """Class for shock data: principal Hugoniot; impedance matched reshock and release. Use dynamic PT class for ramps. 
        Uncertainties:
            _std is one standard deviation
            add asymmetric uncertainties as needed
            _ustd is upper standard deviation
            _lstd is lower standard deviation
    """
    def __init__(self):
        # primary database material information
        self.matid = []         # Database Material ID number = IHED ID for those in IHED
        self.material = []      # Database Material name
        self.formula = []       # Database Material chemical formula
        self.rho_ref = []       # Database material reference density
        self.rho_ref_std = []   
        self.rho_ref_units = []
        self.phase_ref = []     # short name of phase/crystal structure of reference state
        self.orientation = []
        # initial state information
        self.m = []            # initial porosity [rhoref/initial density]
        self.rho0 = []         # initial density [for gases in ihed; for precompressed or high-pressure phases for other work]
        self.rho0_std = []          
        self.rho0_units = []
        self.P0 = []           # initial pressure [for gases in ihed]
        self.P0_std = []          
        self.P0_units = [] 
        self.T0 = []           # initial temperature
        self.T0_std = []        
        self.T0_units = []
        self.phase0 = []       # initial phase if not the same as reference state material
        # first shock information
        self.exp_id = []       # experiment ID (e.g., shot number)
        self.phase_id = []     # short name of phase ID in this portion of the Hugoniot (e.g., XRD info in shock)
        self.uflyer = []       # impact velocity of the flyer plate
        self.uflyer_std = []
        self.uflyer_units = []
        self.up = []           # particle velocity
        self.up_std = []        
        self.up_units = []
        self.us = []           # shock velocity
        self.us_std = []       
        self.us_units = []
        self.P = []            # shock pressure
        self.P_std = []        
        self.P_units = []
        self.rr0 = []          # R/R0
        self.rr0_std = []       
        self.rr0_units = []    # dimensionless
        self.rho = []          # shock density
        self.rho_std = []       
        self.rho_units = []
        self.U = []            # sp. internal energy
        self.U_std = []        
        self.U_units = []
        self.T = []            # temperature in shock
        self.T_std = []          
        self.T_units = []
        self.T_b = []            # alternate temperature; make sure explained in remarks
        self.T_b_std = []          
        self.T_b_units = []
        self.emissivity = [] # emissivity for shock temperature
        self.reflectivity = [] # reflectivity for shock temperature
        self.cs = []           # sound speed
        self.cs_std = []      
        self.cs_units = []
        self.S = []            # sp. entropy in shock
        self.S_std = []          
        self.S_units = []
        # off Hugoniot states
        # release=rel
        self.T_rel = []            # temperature 
        self.T_rel_std = []    
        self.T_rel_units = []
        self.up_rel = []           # particle velocity
        self.up_rel_std = []        
        self.up_rel_units = []
        self.P_rel = []            # shock pressure
        self.P_rel_std = []        
        self.P_rel_units = []
        self.rho_rel = []          # shock density
        self.rho_rel_std = []       
        self.rho_rel_units = []
        # mixed phases
        self.phasefrac = []
        # off Hugoniot states
        # reshock=rs
        self.T_rs = []            # temperature
        self.T_rs_std = []    
        self.T_rs_units = []
        self.up_rs = []           # particle velocity
        self.up_rs_std = []        
        self.up_rs_units = []
        self.P_rs = []            # shock pressure
        self.P_rs_std = []        
        self.P_rs_units = []
        self.rho_rs = []          # shock density
        self.rho_rs_std = []       
        self.rho_rs_units = []
        # metadata
        self.ref = []          # reference for data point = index for reflist
        self.rem = []          # remark footnote number
        self.quality = []      # quality flag: 0=no review; x=number of reviews; reviews are entered in action; -1: suspect/problem point
        self.reflist = []      # list of references (Rusbank)
        self.refidlist = []    # list of references ID numbers (Rusbank)
        self.remlist = []      # list of remarks (Rusbank); notes on calibration go here
        self.denotions = []    # list of denotions (Rusbank)
        self.groups = []       # Rusbank material groups
        self.collections = []  # Collection name for grouping data of different types for a single material or technique
        self.actions = []      # append actions: 'create', 'review', 'add', 'remove', 'revise'
                               # action format: [name, action, date, note]                                
        self.datatype = []     # data type in this MongoDB document
        self.platform = []     # 'gas gun', 'Z machine', 'NIF', etc.
        self.diagnostics = []  # 'streak camera','VISAR','PDV', 'XRD', etc.
        self.calcinfo = []     # description of first principles, thermodynamic calculation, or empirical/fitted calculation
        self.methods = []      # 'experimental', 'first principles', 'thermodynamic', 'fitted'
    def qtzhug(self,qtz_liq_hug):
        # pass in the Hugoniot object with us and T fields populated
        # define liquid Hugoniot
        
        qtz_hug_KD      = HPeos_Shock()
        qtz_hug_KD.rho_ref = 2.650 # g/cm3
        qtz_hug_KD.up   = np.arange(2.5,10.,.001) # km/s
        if qtz_liq_hug == 'KD13':
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

        self.P = np.interp(self.us,qtz_hug_KD.us,qtz_hug_KD.P) # GPa
        self.up = np.interp(self.us,qtz_hug_KD.us,qtz_hug_KD.up) # km/s
        self.rho = np.interp(self.us,qtz_hug_KD.us,qtz_hug_KD.rho) # g/cm3
        self.P_units = ['GPa']
        self.up_units = ['km/s']
        self.rho_units = ['g/cm3']

        if 1:
            # metastable region
            # us-up fit Us = res[1]+res[0]*up
            #resm = [1.30952823, 3.72635724] # old liquidus
            resm = [1.1975675,  4.21474804] # Boslough zcut
            #ind = np.where((self.us < 10.6) & (self.us > 9.5))[0]
            ind = np.where((self.T < 6200) & (self.us > 9.5))[0]
            self.up[ind] = (self.us[ind]-resm[1])/resm[0]
            # now solve for the rest of the variables
            self.P[ind] = self.rho_ref[0]*self.us[ind]*self.up[ind]
            self.rho[ind] = 1./( (1./self.rho_ref[0])*(1-self.up[ind]/self.us[ind]) )
        
        # cacl2 region
        # us-up fit Us = res[1]+res[0]*up
        res = [1.68838078, 1.79714477]
        ind = np.where(self.us <= 9.5)
        self.up[ind] = (self.us[ind]-res[1])/res[0]
        # now solve for the rest of the variables
        self.P[ind] = self.rho_ref[0]*self.us[ind]*self.up[ind]
        self.rho[ind] = 1./( (1./self.rho_ref[0])*(1-self.up[ind]/self.us[ind]) )
    def __iter__(self):
        '''Needed to convert class object into dictionary. Only prints non-empty lists.'''
        for key in self.__dict__:
            if len(getattr(self,key)):
                yield key, getattr(self,key)
#
class HPeos_MeltCurve:
    """Class for Melt Curves information."""
    def __init__(self):
        # primary database material information
        self.matid       = []  # Database Material ID number = IHED ID for those in IHED
        self.material    = []  # Database Material name
        self.formula     = []  # Database Formula
        # melt curve parameters
        self.P = []            # pressure single component system (usually BOUNDED by other data in this structure)
        self.P_std = []         
        self.P_units = []
        self.T = []            # temperature single component system (usually BOUNDED by other data in this structure)
        self.T_std = []         
        self.T_units = []
        self.P_l = []          # pressure liquidus
        self.P_l_std = []         
        self.P_l_units = []
        self.P_s = []          # pressure solidus
        self.P_s_std = []          
        self.P_s_units = []
        self.rho_l = []        # liquidus density
        self.rho_l_std = []        
        self.rho_l_units = []
        self.rho_s = []        # solidus density
        self.rho_s_std = []        
        self.rho_s_units = []
        self.U_l = []          # sp. internal energy liquidus
        self.U_l_std = []         
        self.U_l_units = []
        self.U_s = []          # sp. internal energy solidus
        self.U_s_std = []         
        self.U_s_units = []
        self.T_l = []          # temperature liquidus
        self.T_l_std = []          
        self.T_l_units = []
        self.T_s = []          # temperature solidus
        self.T_s_std = []         
        self.T_s_units = []
        self.S_l = []          # sp. entropy liquidus
        self.S_l_std = []         
        self.S_l_units = []
        self.S_s = []          # sp. entropy solidus
        self.S_s_std = []       
        self.S_s_units = []
        self.P_liquid = []          # pressure liquid field phase ID
        self.P_liquid_std = []         
        self.P_liquid_units = []
        self.T_liquid = []          # temperature liquid field phase ID
        self.T_liquid_std = []         
        self.T_liquid_units = []
        self.P_solid = []          # pressure solid field phase ID
        self.P_solid_std = []         
        self.P_solid_units = []
        self.T_solid = []          # temperature solid field phase ID
        self.T_solid_std = []         
        self.T_solid_units = []
        self.phase_solid = []  # short name of solid phase in this portion of boundary
        # metadata
        self.ref = []        # reference for data point = index for reflist
        self.rem = []        # remark footnote number
        self.quality = []    # quality flag: 0=no review; x=number of reviews; reviews are entered in action
        self.reflist = []    # list of references
        self.refidlist = []  # list of references ID numbers
        self.remlist = []    # list of remarks; include fits here
        self.denotions = []  # list of denotions; define phases here
        self.groups = []     # Rusbank material groups
        self.collections = []
        self.actions = []     # append actions: 'create', 'review', 'add', 'remove', 'revise'
                              # action format: [name, action, date, note]                                
        self.datatype = [] # data type in this MongoDB document
        self.platform = [] # 'gas gun', 'Z machine', 'NIF', etc.
        self.diagnostics = [] # 'XRD','Pyrometry','SVS','VISAR','PDV', etc.
        self.calcinfo = []     # description of first principles, thermodynamic calculation, or empirical/fitted calculation
        self.methods = []      # 'experimental', 'first principles', 'thermodynamic', 'fitted'
    def __iter__(self):
        '''Needed to convert class object into dictionary. Only prints non-empty lists.'''
        for key in self.__dict__:
            if len(getattr(self,key)):
                yield key, getattr(self,key)
#
class HPeos_PhaseBoundary:
    """Class for generic Phase Boundary a<->b experimental data/fit/model."""
    def __init__(self):
        # primary database material information
        self.matid       = []  # Database Material ID number = IHED ID for those in IHED
        self.material    = []  # Database Material name
        self.formula     = []  # Database Formula
        # phase boundary parameters
        self.phasename_a = []  # short name of phase a
        self.phasename_b = []  # short name of phase b
        self.structure_a = []  # crystal structure
        self.structure_b = []  # crystal structure
        self.P = []            # pressure boundary a-b USUALLY BOUNDED BY OTHER DATA IN THIS STRUCTURE
        self.P_std = []         
        self.P_units = []
        self.T = []           # temperature boundary a-b USUALLY BOUNDED BY OTHER DATA IN THIS STRUCTURE
        self.T_std = []         
        self.T_units = []
        self.rho_a = []        # density phase a
        self.rho_a_std = []        
        self.rho_a_units = []
        self.rho_b = []        # density phase b
        self.rho_b_std = []        
        self.rho_b_units = []
        self.U_a = []          # sp. internal energy phase a
        self.U_a_std = []         
        self.U_a_units = []
        self.U_b = []          # sp. internal energy phase b
        self.U_b_std = []         
        self.U_b_units = []
        self.S_a = []          # sp. entropy phase a
        self.S_a_std = []         
        self.S_a_units = []
        self.S_b = []          # sp. entropy phae b
        self.S_b_std = []       
        self.S_b_units = []
        self.P_a = []          # pressure field phase a ID
        self.P_a_std = []         
        self.P_a_units = []
        self.T_a = []          # temperature field phase a ID
        self.T_a_std = []         
        self.T_a_units = []
        self.P_b = []          # pressure field phase b ID
        self.P_b_std = []         
        self.P_b_units = []
        self.T_b = []          # temperature field phase b ID
        self.T_b_std = []         
        self.T_b_units = []
        # metadata
        self.ref = []        # reference for data point = index for reflist
        self.rem = []        # remark footnote number
        self.quality = []    # quality flag: 0=no review; x=number of reviews; reviews are entered in action
        self.reflist = []    # list of references
        self.refidlist = []  # list of references ID numbers
        self.remlist = []    # list of remarks; include fits here
        self.denotions = []  # list of denotions; define phases here
        self.groups = []     # Rusbank material groups
        self.collections = []
        self.actions = []     # append actions: 'create', 'review', 'add', 'remove', 'revise'
                              # action format: [name, action, date, note]                                
        self.datatype = [] # data type in this MongoDB document
        self.platform = [] # 'gas gun', 'Z machine', 'NIF', etc.
        self.diagnostics = [] # 'XRD','Pyrometry','SVS','VISAR','PDV', etc.
        self.calcinfo = []     # description of first principles, thermodynamic calculation, or empirical/fitted calculation
        self.methods = []      # 'experimental', 'first principles', 'thermodynamic', 'fitted'
    def __iter__(self):
        '''Needed to convert class object into dictionary. Only prints non-empty lists.'''
        for key in self.__dict__:
            if len(getattr(self,key)):
                yield key, getattr(self,key)
#
class HPeos_Reference:
    """Reference class for ShockDB"""
    def __init__(self):
        self.refid = []     # Reference ID number; integer to keep track of Rusbank references without doi numbers
        self.reference = [] # List with a single text string of full reference in plain ascii as given in Rusbank
        self.doi = []       # List of single doi string for journal articles
        self.isbn = []      # List of single isbn string for books
        self.pid = []       # if no doi or isbn, make a DOE persistent identifier (PID) for this object https://www.osti.gov/pids/
        self.bibtex = []    # List of single string containing full bibtex entry for this object
        self.bibkey = []    # bibtex key in bibtex entry; useful if collating the database bibtex into one file; make sure unique
        self.actions = []   # append actions: 'create', 'review', 'add', 'remove', 'revised'
                            # action format: [name, action, date, note]                                
    def __iter__(self):
        '''Needed to convert class object into dictionary. Only prints non-empty lists.'''
        for key in self.__dict__:
            if len(getattr(self,key)):
                yield key, getattr(self,key)
#
def printhpeoshug(data):
    '''Print Rusbank-like ascii file from mongodb dictionary.'''
    print(data['material']+', R0 = '+str(data['rhoref'][0])+' g/cm3')
    print('     m         up         us          P        rho   ref')
    print('     -       km/s       km/s        GPa      g/cm3     -')
    print('--------------------------------------------------------')
    for i in range(len(data['up'])):
        print(f"{data['m'][i]:6.4}",f"{data['up'][i]:10.6}",f"{data['us'][i]:10.6}",f"{data['P'][i]:10.6}",f"{data['rho'][i]:10.6}",f"{data['ref'][i]:5}")
    print('\nReferences:')
    for i in range(1,len(data['reflist'])):
        print(str(i)+'. '+data['reflist'][i])
#
def plothpeoshug(data):
    '''Plot from mongodb dictionary using matplotlib.'''
    datasub = {'m':data['m'],'up':data['up'],'us':data['us'],'ref':data['ref']}
    df = pd.DataFrame.from_dict(datasub)
    ind = df['m']==1.
    plt.scatter(df['up'][ind],df['us'][ind])
    plt.title(data['material']+' rhoref= '+str(data['rhoref'][0])+' g/cm3')
    plt.xlabel('Particle Velocity ['+data['up_units'][0]+']')
    plt.ylabel('Shock Velocity ['+data['us_units'][0]+']')
#
def hpeoshugtodf(data):
    '''Convert mongodb dictionary to pandasdataframe.'''
    data_dict = {'m':np.asarray(data['m']),'up':np.asarray(data['up']),'us':np.asarray(data['us']),'ref':np.asarray(data['ref'])}
    df = pd.DataFrame.from_dict(data_dict)
    return df
#
def rusbanktodf(rusbankhug_list):
    '''Make a big dataframe of all the data in the rusbankhug_list for plotting.'''
    count=0
    materialall = []
    for data in rusbankhug_list:
        #print(data.m)
        if len(data.m)>0:
            refidlist = []
            for item in data.ref:
                refidlist.append(data.refidlist[int(item)])
            referencelist = []
            for item in data.ref:
                referencelist.append(data.reflist[int(item)])
            if count == 0:
                matidarr = np.full(len(data.m),data.matid)
                matnamelist = data.material * len(data.m)
                if len(data.formula) > 0:
                    formulalist = data.formula * len(data.m)
                else:
                    formulalist = np.full(len(data.m),'')
                grouplist = [data.groups[0]] * len(data.m) # only selecting 1st group for now
                #print(data.matid,len(data.m),len(referencelist),len(formulalist),len(refidlist))
                data_dict = \
                    {'material':matnamelist,'matid':matidarr,'m':np.asarray(data.m),\
                     'up':np.asarray(data.up),'us':np.asarray(data.us),\
                     'P':np.asarray(data.P),'rho':np.asarray(data.rho),\
                     'formula':formulalist,'group':grouplist,\
                     'refid':refidlist, 'reference':referencelist,\
                     'ref':np.asarray(data.ref)}
                dfall = pd.DataFrame.from_dict(data_dict)
                materialall.append(data.material[0])
                count = count+1
            else:
                matidarr = np.full(len(data.m),data.matid)
                matnamelist = data.material * len(data.m)
                if len(data.formula) > 0:
                    formulalist = data.formula * len(data.m)
                else:
                    formulalist = np.full(len(data.m),'')
                grouplist = [data.groups[0]] * len(data.m) # only selecting 1st group for now
                #print(data.matid,len(data.m),len(referencelist),len(formulalist),len(refidlist))
                data_dict = \
                    {'material':matnamelist,'matid':matidarr,'m':np.asarray(data.m),\
                     'up':np.asarray(data.up),'us':np.asarray(data.us),\
                     'P':np.asarray(data.P),'rho':np.asarray(data.rho),\
                     'formula':formulalist,'group':grouplist,\
                     'refid':refidlist, 'reference':referencelist,\
                     'ref':np.asarray(data.ref)}
                df = pd.DataFrame.from_dict(data_dict)
                dfall = pd.concat([dfall,df],ignore_index=True)
                materialall.append(data.material[0])
    return dfall,materialall
#
def hpeoshugtodf(shockdbhug_cursor):
    '''Make a big dataframe of all the data from shockdb for plotting. Pymongo objects are dictionaries.'''
    count=0
    materialall = []
    for data in shockdbhug_cursor:
        #print(data.m)
        if len(data['m'])>0:
            refidlist = []
            for item in data['ref']:
                refidlist.append(data['refidlist'][int(item)])
            referencelist = []
            for item in data['ref']:
                referencelist.append(data['reflist'][int(item)])
            if count == 0:
                matidarr = np.full(len(data['m']),data['matid'])
                matnamelist = data['material'] * len(data['m'])
                formulalist = data['formula'] * len(data['m'])
                grouplist = [data['groups'][0]] * len(data['m']) # only selecting 1st group for now
                #print(data.matid,len(data.m),len(referencelist),len(formulalist),len(refidlist))
                data_dict = \
                    {'material':matnamelist,'matid':matidarr,'m':np.asarray(data['m']),\
                     'up':np.asarray(data['up']),'us':np.asarray(data['us']),\
                     'P':np.asarray(data['P']),'rho':np.asarray(data['rho']),\
                     'formula':formulalist,'group':grouplist,\
                     'refid':refidlist, 'reference':referencelist}
                dfall = pd.DataFrame.from_dict(data_dict)
                materialall.append(data['material'][0])
                count = count+1
            else:
                matidarr = np.full(len(data['m']),data['matid'])
                matnamelist = data['material'] * len(data['m'])
                formulalist = data['formula'] * len(data['m'])
                grouplist = [data['groups'][0]] * len(data['m']) # only selecting 1st group for now
                #print(data.matid,len(data.m),len(referencelist),len(formulalist),len(refidlist))
                data_dict = \
                    {'material':matnamelist,'matid':matidarr,'m':np.asarray(data['m']),\
                     'up':np.asarray(data['up']),'us':np.asarray(data['us']),\
                     'P':np.asarray(data['P']),'rho':np.asarray(data['rho']),\
                     'formula':formulalist,'group':grouplist,\
                     'refid':refidlist, 'reference':referencelist}
                df = pd.DataFrame.from_dict(data_dict)
                dfall = pd.concat([dfall,df],ignore_index=True)
                materialall.append(data['material'][0])
    return dfall,materialall
#
# End of shockdb_module.py file
#
# Additional classes for Silica Melt Curve Paper
class HPeos_Model:
    """Class for EOS model isotherms and isentropes
    """
    def __init__(self):
        self.rhoref = []
        self.K0 = []
        self.Kp = []
        self.g0 = []
        self.q = []
        self.rhoarr = []
        self.rho00 = [] # 0K reference density
        self.varr = [] # = 1/rhoarr cgs
        self.vaarr = [] # volume in mol/anstrom^3 for Karki
        self.garr = []
        self.p0000 = []
        self.p3000 = []
        self.p3130 = []
        self.p3500 = []
        self.p4000 = []
        self.p4500 = []
        self.p4800 = []
        self.p5000 = []
        self.p5075 = []
        self.p5500 = []
        self.p5900 = []
        self.p6000 = []
        self.p6100 = []
        self.p6500 = []
        self.p7000 = []
        self.t3000 = []
        self.t3500 = []
        self.t4000 = []
        self.t4500 = []
        self.t4800 = []
        self.t5000 = []
        self.t5075 = []
        self.t5500 = []
        self.t5900 = []
        self.t6000 = []
        self.t6100 = []
        self.t6500 = []
        self.t7000 = []
        self.refpt = [] # (rho,p,t,s)
        self.cv = [] # constant high temp cv value for solid eos calcs
        self.cvarr = [] # cv is a function of density for liquid
        self.st0000 = [] # temperatures along the isentropes
        self.st3000 = [] # temperature along rhoarr-st3000 isentrope; vals in svalsarr
        self.st3130 = [] # temperature along rhoarr-st3130 isentrope; vals in svalsarr
        self.st3500 = [] # temperature along rhoarr-st3000 isentrope; vals in svalsarr
        self.st4000 = []
        self.st4500 = []
        self.st4800 = []
        self.st5000 = []
        self.st5075 = []
        self.st5500 = []
        self.st5900 = []
        self.st6000 = []
        self.st6100 = []
        self.st6500 = []
        self.st7000 = []
        self.sp0000 = [] # pressures along the isentropes
        self.sp3000 = [] 
        self.sp4000 = []
        self.sp4500 = []
        self.sp4800 = []
        self.sp5000 = []
        self.sp5075 = []
        self.sp5500 = []
        self.sp5900 = []
        self.sp6000 = []
        self.sp6100 = []
        self.sp6500 = []
        self.sp7000 = []
        self.su0000 = [] # internal energy along the isentropes
        self.su3000 = [] 
        self.su4000 = []
        self.su4500 = []
        self.su4800 = []
        self.su5000 = []
        self.su5075 = []
        self.su5500 = []
        self.su5900 = []
        self.su6000 = []
        self.su6500 = []
        self.su7000 = []
        self.qtzrhotest = []
        self.fsrhotest = []
        self.stempsarr = []
        self.svalsarr = []
        self.uvalsarr = []
        self.isenval_p = []
        self.isenval_t = []
        self.isenval_s = []
        self.isenval_u = []
        self.isotherm_p = []
        self.isotherm_t = []
        self.isotherm_s = []
    def isenval_liquid(self,sval):
        # at the reference density, find the temperature that corresponds to the desired sp. entropy
        # only for the Karki et al. 2007 PVT model (shifted in density)
        karkireftemporiginal=3000
        tstart = np.exp((sval-self.refpt[3])/self.cv+np.log(self.refpt[2]))
        #print('S=',sval,', T=',tstart)
        rhoref = self.refpt[0]
        rhoarr = deepcopy(self.rhoarr)
        garr = deepcopy(self.garr)
        tarr = np.zeros(len(rhoarr))
        parr = np.zeros(len(rhoarr))
        uarr = np.zeros(len(rhoarr))
        t0 = tstart
    
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
        BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-self.vaarr)/1.86)))) # GPa
        parr = 1.e9*(self.p3000+BV*(tarr-karkireftemporiginal)) # karki thermal model reference temp 3000 K; Pa
        # now that we have P, calculate dU along the isentrope dS=0
        # dU = TdS-PdV
        dU = self.cv*(tstart-self.refpt[2]) # reference point for liquid J/kg
        uarr[ind[0]]=dU # U=0 at refpt; this is offset
        for iii in range(ind[1],np.max(ind)+1):
            uarr[iii] = uarr[iii-1] - parr[iii]*(1./(1000*rhoarr[iii])-1./(1000*rhoarr[iii-1]))
        rind = np.flip(np.where(rhoarr < rhoref)[0])
        for iii in rind:
            uarr[iii] = uarr[iii+1] - parr[iii]*(1./(1000*rhoarr[iii])-1./(1000*rhoarr[iii+1]))
        self.isenval_p = deepcopy(parr/1.e9) # GPa
        self.isenval_t = deepcopy(tarr) # K
        self.isenval_s = sval # J/K/kg
        self.isenval_u = deepcopy(uarr) # J/kg
    def isenval_stv(self,sval):
        # stv and sft have different entropy tie points so separate functions
        # only for stv reference 4500 K
        # at the reference density, find the temperature that corresponds to the desired sp. entropy
        tstart = np.exp((sval-self.refpt[3])/self.cv+np.log(self.refpt[2]))
        #print('S=',sval,', T=',tstart)
        #
        rhoref = self.refpt[0]
        rhoarr = deepcopy(self.rhoarr)
        garr = deepcopy(self.garr)
        tarr = np.zeros(len(rhoarr))
        parr = np.zeros(len(rhoarr))
        uarr = np.zeros(len(rhoarr))
        t0 = tstart
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
        dE = self.cv*(tarr-self.refpt[2]) #  reference curve FOR ENTROPY # mks
        parr = self.p4000+ dE*garr*(self.rhoarr*1000) # Pa 
        # now that we have P, calculate dU along the isentrope dS=0
        # dU = TdS-PdV
        dU = self.cv*(tstart-self.refpt[2]) # reference point for stv
        uarr[ind[0]]=dU # U=0 at refpt; this is offset
        for iii in range(ind[1],np.max(ind)+1):
            uarr[iii] = uarr[iii-1] - parr[iii]*(1./(1000*rhoarr[iii])-1./(1000*rhoarr[iii-1]))
        rind = np.flip(np.where(rhoarr < rhoref)[0])
        for iii in rind:
            uarr[iii] = uarr[iii+1] - parr[iii]*(1./(1000*rhoarr[iii])-1./(1000*rhoarr[iii+1]))
        self.isenval_p = deepcopy(parr/1.e9) # GPa
        self.isenval_t = deepcopy(tarr) # K
        self.isenval_s = sval # J/K/kg
        self.isenval_u = deepcopy(uarr) # J/kg
    def isenval_sft(self,sval):
        # stv and sft have different entropy tie points so separate functions
        # only for sft reference 5900 K
        # at the reference density, find the temperature that corresponds to the desired sp. entropy
        tstart = np.exp((sval-self.refpt[3])/self.cv+np.log(self.refpt[2]))
        #print('sft S=',sval,', T=',tstart)
        #
        rhoref = self.refpt[0]
        rhoarr = deepcopy(self.rhoarr)
        garr = deepcopy(self.garr)
        tarr = np.zeros(len(rhoarr))
        uarr = np.zeros(len(rhoarr))
        t0 = tstart
        
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
        dE = self.cv*(tarr-self.refpt[2]) # Temp on ref curve
        parr = self.p5900+ dE*garr*(rhoarr*1000) # Pa  
        
        # now that we have P, calculate dU along the isentrope dS=0
        # dU = TdS-PdV
        dU = self.cv*(tstart-self.refpt[2]) # reference point for stv
        uarr[ind[0]]=dU # U=0 at refpt; this is offset
        for iii in range(ind[1],np.max(ind)+1):
            uarr[iii] = uarr[iii-1] - parr[iii]*(1./(1000*rhoarr[iii])-1./(1000*rhoarr[iii-1]))
        rind = np.flip(np.where(rhoarr < rhoref)[0])
        for iii in rind:
            uarr[iii] = uarr[iii+1] - parr[iii]*(1./(1000*rhoarr[iii])-1./(1000*rhoarr[iii+1]))
        self.isenval_p = deepcopy(parr/1.e9) # GPa
        self.isenval_t = deepcopy(tarr) # K
        self.isenval_s = sval # J/K/kg
        self.isenval_u = deepcopy(uarr) # J/kg
    def isotherm_stv(self,tval):
        # stv and sft have different entropy tie points so separate functions
        dE = self.cv*(tval-4000) # 4000 is the reference curve for stv # mks
        pth = dE*self.garr*(self.rhoarr*1000)/1.e9 # Pa->GPa
        self.isotherm_p = self.p4000+pth
        self.isotherm_t = tval
        self.isotherm_s = np.zeros(len(self.rhoarr))
        # now calculate the entropy along this isotherm
        # for now, do it the same way the solidus is calculated
        # - by interolation on a grid
        tfitarr = [3000,3500,4000,4500,5000,5500,6000,6500]
        pfitarr = np.zeros(len(tfitarr))
        sfitarr = np.zeros(len(tfitarr))
        for iii in range(0,len(tfitarr)):
            ind = np.where(self.stempsarr == tfitarr[iii])[0]
            sfitarr[iii] = self.svalsarr[ind]   
        #print('isotherm calc: ',tval,sfitarr)
        for iii in range(0,len(self.rhoarr)):
            # pressure density point on the set of isentropes pinned to rhoref
            pfitarr[0] = self.sp3000[iii]
            pfitarr[1] = self.sp3500[iii]
            pfitarr[2] = self.sp4000[iii]
            pfitarr[3] = self.sp4500[iii]
            pfitarr[4] = self.sp5000[iii]
            pfitarr[5] = self.sp5500[iii]
            pfitarr[6] = self.sp6000[iii]
            pfitarr[7] = self.sp6500[iii]
            self.isotherm_s[iii]=np.interp(self.isotherm_p[iii],pfitarr,sfitarr)
    def isotherm_sft(self,tval):
        # stv and sft have different entropy tie points so separate functions
        dE = self.cv*(tval-4000) # 4000 is the reference curve for sft # mks
        pth = dE*self.garr*(self.rhoarr*1000)/1.e9 # Pa->GPa
        self.isotherm_p = self.p4000+pth
        self.isotherm_t = tval
        self.isotherm_s = np.zeros(len(self.rhoarr))
        # now calculate the entropy along this isotherm
        # for now, do it the same way the solidus is calculated
        # - by interolation on a grid
        tfitarr = [3000,3500,4000,4500,5000,5500,6000,6500]
        pfitarr = np.zeros(len(tfitarr))
        sfitarr = np.zeros(len(tfitarr))
        for iii in range(0,len(tfitarr)):
            ind = np.where(self.stempsarr == tfitarr[iii])[0]
            sfitarr[iii] = self.svalsarr[ind]   
        #print('isotherm calc: ',tval,sfitarr)
        for iii in range(0,len(self.rhoarr)):
            # pressure density point on the set of isentropes pinned to rhoref
            pfitarr[0] = self.sp3000[iii]
            pfitarr[1] = self.sp3500[iii]
            pfitarr[2] = self.sp4000[iii]
            pfitarr[3] = self.sp4500[iii]
            pfitarr[4] = self.sp5000[iii]
            pfitarr[5] = self.sp5500[iii]
            pfitarr[6] = self.sp6000[iii]
            pfitarr[7] = self.sp6500[iii]
            self.isotherm_s[iii]=np.interp(self.isotherm_p[iii],pfitarr,sfitarr)
    def isotherm_liquid(self,tval):
        sio2mol = 28.085+2*15.999 # g per formula unit
        #dE = self.cv*(tval-3000) # 3000 is the reference curve for liquid # mks
        #pth = dE*self.garr*(self.rhoarr*1000)/1.e9 # Pa->GPa
        #vaarr = sio2mol*1.6605391/self.rhoarr
        T0=3000
        BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-self.vaarr)/1.86)))) # GPa
        #BV = 1.e-3*(10.4 - (9.5 / (1+np.exp((25.2-vaarr)/1.86)))) # GPa
        ptherm =  BV*(tval-T0)
        self.isotherm_p = self.p3000+ptherm
        self.isotherm_t = tval
        self.isotherm_s = np.zeros(len(self.rhoarr))
        # now calculate the entropy along this isotherm
        # for now, do it the same way the solidus is calculated
        # - by interolation on a grid
        tfitarr = [3000,3500,4000,4500,5000,5500,6000,6500]
        pfitarr = np.zeros(len(tfitarr))
        sfitarr = np.zeros(len(tfitarr))
        for iii in range(0,len(tfitarr)):
            ind = np.where(self.stempsarr == tfitarr[iii])[0]
            sfitarr[iii] = self.svalsarr[ind]   
        #print('isotherm calc: ',tval,sfitarr)
        for iii in range(0,len(self.rhoarr)):
            # pressure density point on the set of isentropes pinned to rhoref
            pfitarr[0] = self.sp3000[iii]
            pfitarr[1] = self.sp3500[iii]
            pfitarr[2] = self.sp4000[iii]
            pfitarr[3] = self.sp4500[iii]
            pfitarr[4] = self.sp5000[iii]
            pfitarr[5] = self.sp5500[iii]
            pfitarr[6] = self.sp6000[iii]
            pfitarr[7] = self.sp6500[iii]
            self.isotherm_s[iii]=np.interp(self.isotherm_p[iii],pfitarr,sfitarr)
# End of classes for silica melt curve paper

#==================================================================================================
#### INTERSECTION FUNCTION FROM https://github.com/sukhbinder/intersection/
#### MIT LICENSE
#"""
#Give, two x,y curves this gives intersection points,
#autor: Sukhbinder
#5 April 2017
#Based on: http://uk.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
#"""
def _rect_inter_inner(x1, x2):
    n1 = x1.shape[0]-1
    n2 = x2.shape[0]-1
    X1 = np.c_[x1[:-1], x1[1:]]
    X2 = np.c_[x2[:-1], x2[1:]]
    S1 = np.tile(X1.min(axis=1), (n2, 1)).T
    S2 = np.tile(X2.max(axis=1), (n1, 1))
    S3 = np.tile(X1.max(axis=1), (n2, 1)).T
    S4 = np.tile(X2.min(axis=1), (n1, 1))
    return S1, S2, S3, S4
def _rectangle_intersection_(x1, y1, x2, y2):
    S1, S2, S3, S4 = _rect_inter_inner(x1, x2)
    S5, S6, S7, S8 = _rect_inter_inner(y1, y2)

    C1 = np.less_equal(S1, S2)
    C2 = np.greater_equal(S3, S4)
    C3 = np.less_equal(S5, S6)
    C4 = np.greater_equal(S7, S8)

    ii, jj = np.nonzero(C1 & C2 & C3 & C4)
    return ii, jj
def Intersection(x1, y1, x2, y2):
    """
    INTERSECTIONS Intersections of curves.
       Computes the (x,y) locations where two curves intersect.  The curves
       can be broken with NaNs or have vertical segments.
    usage:
    x,y=Intersection(x1,y1,x2,y2)
        Example:
        a, b = 1, 2
        phi = np.linspace(3, 10, 100)
        x1 = a*phi - b*np.sin(phi)
        y1 = a - b*np.cos(phi)
        x2=phi
        y2=np.sin(phi)+2
        x,y=intersection(x1,y1,x2,y2)
        plt.plot(x1,y1,c='r')
        plt.plot(x2,y2,c='g')
        plt.plot(x,y,'*k')
        plt.show()
    """
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)

    ii, jj = _rectangle_intersection_(x1, y1, x2, y2)
    n = len(ii)

    dxy1 = np.diff(np.c_[x1, y1], axis=0)
    dxy2 = np.diff(np.c_[x2, y2], axis=0)

    T = np.zeros((4, n))
    AA = np.zeros((4, 4, n))
    AA[0:2, 2, :] = -1
    AA[2:4, 3, :] = -1
    AA[0::2, 0, :] = dxy1[ii, :].T
    AA[1::2, 1, :] = dxy2[jj, :].T

    BB = np.zeros((4, n))
    BB[0, :] = -x1[ii].ravel()
    BB[1, :] = -x2[jj].ravel()
    BB[2, :] = -y1[ii].ravel()
    BB[3, :] = -y2[jj].ravel()

    for i in range(n):
        try:
            T[:, i] = np.linalg.solve(AA[:, :, i], BB[:, i])
        except:
            T[:, i] = np.Inf

    in_range = (T[0, :] >= 0) & (T[1, :] >= 0) & (
        T[0, :] <= 1) & (T[1, :] <= 1)

    xy0 = T[2:, in_range]
    xy0 = xy0.T
    return xy0[:, 0], xy0[:, 1]
##==================================================================================================

