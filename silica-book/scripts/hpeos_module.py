## Functions to support High Pressure EOS Database
## include this command in python: import hpeos_module as hp
## S. T. Stewart
## sstewa56@asu.edu
#
# Needs: a way to report fits with covariances; this can be a separate structure
# but it needs a customized format? ok for MongoDB -- a text string that provides the function
# parameters list and covariance list; problems can be ordering of parameters
#
__version__='0.3' # 5/9/2025 - set independent repo version control
#__version__='0.2.2' # 8/9/2023 - set independent repo version control
#__version__='0.2.1' # 7/8/2023 - set independent repo version control
#__version__='0.1.1' # 5/13/2023 adding more data types; meltcurve; phase boundaries; NEED generic EOS/thermo info structure
#__version__='0.1' # 4/27/2023
#
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import deepcopy
from dataclasses import dataclass, field, fields
import json
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)
#
# every entry is a list because empty list is a valid object for regexp
# when different materials have different parameters; empty lists do not throw errors
#
######################################################################################
# Index class for database variables
from enum import IntEnum
class Column(IntEnum):
    """Enumeration for datacube columns for improved readability."""
    rho0 = 0
    srho0 = 1
    p0 = 2
    sp0 = 3
    e0 = 4
    h0 = 5
    t0 = 6
    up = 7
    sup = 8
    us = 9
    sus = 10
    p = 11
    sp = 12
    v = 13
    sv = 14
    t = 15
    st = 16
    h = 17
    sh = 18
    e = 19
    se = 20
    refcode = 21
    emis = 22
    semis = 23
    exp = 24
    xi = 25
    up2 = 26
    ec = 27
    cs = 28
    gamma = 29
######################################################################################
# FUNCTIONS FOR FITTING DATA USING SCIPY ODR
#
import scipy
scipyversion = scipy.__version__
del scipy

from scipy.odr import ODR, Model, RealData
def linear_function(a, x):
    '''Linear function f(x) = a[0] + a[1]*x '''
    # A is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return a[0]+a[1]*x
linear_model = Model(linear_function)

def quadratic_function(a, x):
    '''Quadratic function f(x) = a[0] + a[1]*x + a[2]*x^2 '''
    # A is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return a[2]*x*x + a[1]*x + a[0]
quadratic_model = Model(quadratic_function)

def cubic_function(a, x):
    '''Cubic function f(x) = a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 '''
    # A is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return a[3]*x*x*x + a[2]*x*x + a[1]*x + a[0]
cubic_model = Model(cubic_function)

def fourthpoly_function(a, x):
    '''4th order polynomial function f(x) = a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 + a[4]*x^4 '''
    # A is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return a[4]*x*x*x*x +a[3]*x*x*x + a[2]*x*x + a[1]*x + a[0]
fourthpoly_model = Model(fourthpoly_function)

def simons_function(a, x):
    '''Simons function f(x) = a[0] + a[1]*x^a[2]'''
    # A is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return a[0]+a[1]*np.pow(x,a[2])
simons_model = Model(simons_function)

def simons2_function(a, x):
    '''Simons2 function f(x) = 3073 ((x-13.7)/a[0] - 1)^a[1]'''
    # A is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return 3073*np.pow((x-13.7)/a[0]-1,a[1])
simons2_model = Model(simons2_function)

######################################################################################
class HPeos_Material:
    """Class database materials and associated database entries."""
    def __init__(self):
        # primary database material information
        self.matid         = []     # Database Material ID number = IHED ID for those in IHED
        self.material      = []     # Database Material name
        self.aliases       = []     # Other names used for the same material
        self.formula       = []     # Database Material chemical formula
        self.phase         = []     # list of short name(s) for the material phase/crystal structure
        self.orientation   = []     # xtl orientation
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
# https://realpython.com/python-data-classes/
@dataclass
class HPeos_Fit:
    """ Data class for fits to EOS curves.
        v0.1 sts 20250509
    """
    # primary database material information
    filename:      str = ''     # filename base for writing the EOS fit to file
    outputfitpath: str = './'   # output path for fit files
    matid:         int = 0      # Database Material ID number = IHED ID for those in IHED
    material:      str = ''     # Database Material name
    formula:       str = ''     # Database Material chemical formula
    rho_ref:       float = 0.   # Reference density for starting phase
    rho_ref_std:   float = 0.   # Reference density uncertainty for starting phase
    rho_ref_units: str = ''     # Reference density units
    aliases:       list[str]  = field(default_factory=lambda: [])  # Other names used for the same material
    phase:         list[str]  = field(default_factory=lambda: [])  # list of name(s) for the material phase, descriptions of crystal structure
    orientation:   list[str]  = field(default_factory=lambda: [])  # (initial/reference) xtl orientation
    x_var:         str = ''     # standard variable for independent variable x: e.g., Us, up, T, P
    f_var:         str = ''     # standard variable for dependent variable f(x)
    x_units:       str = ''     # units of independent variable 
    f_units:       str = ''     # units of dependent variable
    x_desc:        str = ''     # descriptive name for independent variable x: e.g. for plot labels, Particle Velocity
    f_desc:        str = ''     # descriptive name  for dependent variable f(x): e.g. for plot labels, Shock Velocity
    curvetype:     str = ''     # 'linear' 'quadratic' 'cubic' 'fourthpoly'
    equation:      str = ''     # f(x) function to define the a[i] variables
    a0:            np.ndarray = field(default_factory=lambda: np.zeros(0))    # initial guess for fit coefficients
    a:             np.ndarray = field(default_factory=lambda: np.zeros(0))    # fit coefficients
    cov:           np.ndarray = field(default_factory=lambda: np.zeros(0))    # covariance matrix for fit. 
    xrange:        np.ndarray = field(default_factory=lambda: np.zeros(0))    # np.asarray([x1,x2]) range of independent variable
    npts:          int = 0      # number of x points fit
    rms:           float = 0.   # rms in f_units
    phase_fit:     list[str]  = field(default_factory=lambda: [])  # name of phase corresponding to the fitted region        
    actions:       list[str]  = field(default_factory=lambda: [])  # append actions: 'create', 'review', 'add', 'remove', 'revise'
                                                                   # action format: [name, action, date, note]
    remlist:       list[str]  = field(default_factory=lambda: [])  # comments on the fit
    hpeos_version:  str = __version__     # keep track of hpeos_module.py version for future field and function changes
    scipy_version:  str = scipyversion    # keep track of scipy version for future field and function changes
    # Fields above are formatted as key: value pairs in the output text file
    # THE FOLLOWING FIELDS MUST COME AT THE END OF THE CLASS DEFINITION
    # The last fields in the class are the data table and reference list
    # These are not formatted as key: value pairs in the output file
    # These will be end of file to read back in
    # These last 6 entries are written at the end of the class output file
    # File indicator for the data table: newline + TABLE
    x:             np.ndarray = field(default_factory=lambda: np.zeros(0))   # values of independent variable x
    x_std:         np.ndarray = field(default_factory=lambda: np.zeros(0))   # std err for independent variable x
    f:             np.ndarray = field(default_factory=lambda: np.zeros(0))   # values of dependent variable f(x)
    f_std:         np.ndarray = field(default_factory=lambda: np.zeros(0))   # std err of dependent variable f(x)
    refid:         list[int]  = field(default_factory=lambda: [])            # reference id number for each point f(x)
    # File indicator for the references: newline + REFERENCES
    reflist:       list[str]  = field(default_factory=lambda: [])            # refid: citation list
    def comparecoefs(self,a2):
        test = 0
        for i in range(len(a2)):
            test += abs(a2[i]-self.a[i])
        if test > 1.e-5:
            print('Coefficients are different = ',test)
            print('***** COEFFICIENTS DO NOT MATCH!! ********')
        else:
            print('Coefficients are the same.')
        print(a2)
        print(self.a)
    def calcf(self,x: np.ndarray) -> np.ndarray:
        '''Returns calculation of the f(x) curve given the independent variable x.'''
        if self.curvetype == 'simons':
            return self.a[0]+self.a[1]*np.pow(x,self.a[2])
        if self.curvetype == 'simons2':
            return 3073*np.pow((x-13.7)/self.a[0],self.a[1])
        if self.curvetype == 'fourthpoly':
            return self.a[0]+self.a[1]*x + self.a[2]*np.pow(x,2) + self.a[3]*np.pow(x,3)+ self.a[4]*np.pow(x,4)
        if self.curvetype == 'cubic':
            return self.a[0]+self.a[1]*x + self.a[2]*np.pow(x,2) + self.a[3]*np.pow(x,3)
        if self.curvetype == 'quadratic':
            return self.a[0]+self.a[1]*x + self.a[2]*np.pow(x,2)
        if self.curvetype == 'linear':
            return self.a[0]+self.a[1]*x
    def fit(self):
        if self.curvetype == 'simons':
            self.equation = 'f(x) = a[0] + a[1]*x^a[2]' 
            if len(self.a0)==0:
                self.a0=[1,1,1]
  #np.concatenate(np.flip(np.polyfit(self.x,self.f,1)),np.asarray([1])) # this equation coefs are flipped from polyfit
            data = RealData(self.x, self.f, sx=self.x_std, sy=self.f_std)
            odr = ODR(data, simons_model, beta0=self.a0)
            output=odr.run()
            self.a = output.beta
            self.cov = output.cov_beta
            fff = self.a[0]+self.a[1]*np.pow(self.x,self.a[2])
            self.rms = np.sqrt(np.mean((fff-self.f)**2))
            self.xrange = np.asarray([min(self.x),max(self.x)])
            self.npts = len(self.x)
            self.remlist.append('Fit with scipy ODR version '+scipyversion)
        if self.curvetype == 'simons2':
            self.equation = 'f(x) = Tref*((x-Pref)/a[0])^a[1]' 
            if len(self.a0)==0:
                self.a0=[1,1]
  #np.concatenate(np.flip(np.polyfit(self.x,self.f,1)),np.asarray([1])) # this equation coefs are flipped from polyfit
            data = RealData(self.x, self.f, sx=self.x_std, sy=self.f_std)
            odr = ODR(data, simons2_model, beta0=self.a0)
            output=odr.run()
            self.a = output.beta
            self.cov = output.cov_beta
            fff = 3073*np.pow((self.x-13.7)/self.a[0],self.a[1])
            self.rms = np.sqrt(np.mean((fff-self.f)**2))
            self.xrange = np.asarray([min(self.x),max(self.x)])
            self.npts = len(self.x)
            self.remlist.append('Fit with scipy ODR version '+scipyversion)
        if self.curvetype == 'fourthpoly':
            self.equation = 'f(x) = a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 + a[4]*x^4' 
            a0 = np.flip(np.polyfit(self.x,self.f,4)) # this equation coefs are flipped from polyfit
            data = RealData(self.x, self.f, sx=self.x_std, sy=self.f_std)
            odr = ODR(data, fourthpoly_model, beta0=a0)
            output=odr.run()
            self.a = output.beta
            self.cov = output.cov_beta
            fff = self.a[0]+self.a[1]*self.x+self.a[2]*np.pow(self.x,2)+self.a[3]*np.pow(self.x,3)+self.a[4]*np.pow(self.x,4)
            self.rms = np.sqrt(np.mean((fff-self.f)**2))
            self.xrange = np.asarray([min(self.x),max(self.x)])
            self.npts = len(self.x)
            self.remlist.append('Fit with scipy ODR version '+scipyversion)
        if self.curvetype == 'cubic':
            self.equation = 'f(x) = a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3' 
            a0 = np.flip(np.polyfit(self.x,self.f,3)) # this equation coefs are flipped from polyfit
            data = RealData(self.x, self.f, sx=self.x_std, sy=self.f_std)
            odr = ODR(data, cubic_model, beta0=a0)
            output=odr.run()
            self.a = output.beta
            self.cov = output.cov_beta
            fff = self.a[0]+self.a[1]*self.x+self.a[2]*np.pow(self.x,2)+self.a[3]*np.pow(self.x,3)
            self.rms = np.sqrt(np.mean((fff-self.f)**2))
            self.xrange = np.asarray([min(self.x),max(self.x)])
            self.npts = len(self.x)
            self.remlist.append('Fit with scipy ODR version '+scipyversion)
        if self.curvetype == 'quadratic':
            self.equation = 'f(x) = a[0] + a[1]*x + a[2]*x^2' 
            a0 = np.flip(np.polyfit(self.x,self.f,2)) # this equation coefs are flipped from polyfit
            data = RealData(self.x, self.f, sx=self.x_std, sy=self.f_std)
            odr = ODR(data, quadratic_model, beta0=a0)
            output=odr.run()
            self.a = output.beta
            self.cov = output.cov_beta
            fff = self.a[0]+self.a[1]*self.x+self.a[2]*np.pow(self.x,2)
            self.rms = np.sqrt(np.mean((fff-self.f)**2))
            self.xrange = np.asarray([min(self.x),max(self.x)])
            self.npts = len(self.x)
            self.remlist.append('Fit with scipy ODR version '+scipyversion)
        if self.curvetype == 'linear':
            self.equation = 'f(x) = a[0] + a[1]*x' 
            a0 = np.polyfit(self.x,self.f,1) # this equation coefs are flipped from polyfit
            data = RealData(self.x, self.f, sx=self.x_std, sy=self.f_std)
            odr = ODR(data, linear_model, beta0=a0)
            output=odr.run()
            self.a = output.beta
            self.cov = output.cov_beta
            fff = self.a[0]+self.a[1]*self.x
            self.rms = np.sqrt(np.mean((fff-self.f)**2))   
            self.xrange = np.asarray([min(self.x),max(self.x)])
            self.npts = len(self.x)
            self.remlist.append('Fit with scipy ODR version '+scipyversion)
    def __iter__(self):
        '''Needed to convert class object into dictionary.'''
        for key in self.__dict__:
            yield key, getattr(self,key)
    def savejson(self):
        # write json file that can be read back into the dataclass
        with open(self.outputfitpath+self.filename+'.txt', 'w') as f:
            temp = json.dumps(dict(self),cls=NumpyEncoder)
            f.write(temp) # writes as one long string
        print('Saved EOS fit in '+self.filename+'.txt')
    def savereadable(self):
        # write with new lines after each field for human readability and csv table of data 
        ddd = dict(self)
        nfields = len({field.name for field in fields(self)})
        with open(self.outputfitpath+self.filename+'-readable.txt', 'w') as f:
            #for i in range(len(ddd.keys())):
            for i in range(nfields-6): # do not write the last 6 fields as key:value pairs
                f.write(str(list(ddd.keys())[i])+': '+str(list(ddd.values())[i])+'\n')
            f.write('\n')
            f.write('TABLE\n')
            #f.write('%16s, %16s, %16s, %16s, %6s\n' % ('x', 'x_std', 'f', 'f_std', 'refid'))
            f.write('%20s, %20s, %20s, %20s, %6s\n' % (self.x_desc, 'stdev', self.f_desc, 'stdev', 'refid'))
            f.write('%20s, %20s, %20s, %20s, %6s\n' % (self.x_units, self.x_units, self.f_units,  self.f_units, '-'))
            for i in range(len(self.x)):
#                f.write(str(self.x[i])+', '+str(self.x_std[i])+', '+str(self.f[i])+', '+str(self.f_std[i])+', '+str(self.refid[i])+'\n')
                f.write('%20.8e, %20.8e, %20.8e, %20.8e, %6i\n' % \
                (self.x[i],self.x_std[i],self.f[i],self.f_std[i],self.refid[i]))
            f.write('\n')
            f.write('REFERENCES\n')
            for item in ddd['reflist']:
                f.write(item+'\n')
        print('Saved EOS fit in readable file '+self.filename+'-readable.txt')
    def loadjson(self):
        # read in json formatted fit record file and then reconstruct a dictionary
        # need to add a function to load the readable version
        with open(self.outputfitpath+self.filename+'.txt', 'r') as f:
            filetext = f.readlines()[0] # read in whole file into text
        datadict = json.loads(filetext) # convert long string into dictionary
        # print(datadict)
        # assign dictionary values into dataclass fields manually
        # simplest brute force to assign correct data types
        fff = 'filename'
        if fff in list(datadict.keys()):  
            self.filename =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'outputfitpath'
        if fff in list(datadict.keys()):  
            self.outputfitpath =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'matid'
        if fff in list(datadict.keys()):  
            self.matid =int(list(datadict.values())[list(datadict.keys()).index(fff)])
        fff = 'material'
        if fff in list(datadict.keys()):  
            self.material =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'formula'
        if fff in list(datadict.keys()):  
            self.formula =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'rho_ref'
        if fff in list(datadict.keys()):  
            self.rho_ref =float(list(datadict.values())[list(datadict.keys()).index(fff)])
        fff = 'rho_ref_std'
        if fff in list(datadict.keys()):  
            self.rho_ref_std =float(list(datadict.values())[list(datadict.keys()).index(fff)])
        fff = 'rho_ref_units'
        if fff in list(datadict.keys()):  
            self.rho_ref_units =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'aliases'
        if fff in list(datadict.keys()):  
            self.aliases =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'phase'
        if fff in list(datadict.keys()):  
            self.phase =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'orientation'
        if fff in list(datadict.keys()):  
            self.orientation =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'x_var'
        if fff in list(datadict.keys()):  
            self.x_var =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'f_var'
        if fff in list(datadict.keys()):  
            self.f_var =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'x_units'
        if fff in list(datadict.keys()):  
            self.x_units =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'f_units'
        if fff in list(datadict.keys()):  
            self.f_units =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'x_desc'
        if fff in list(datadict.keys()):  
            self.x_desc =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'f_desc'
        if fff in list(datadict.keys()):  
            self.f_desc =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'curvetype'
        if fff in list(datadict.keys()):  
            self.curvetype =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'equation'
        if fff in list(datadict.keys()):  
            self.equation =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'a'
        if fff in list(datadict.keys()):  
            self.a =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'cov'
        if fff in list(datadict.keys()):  
            self.cov =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'xrange'
        if fff in list(datadict.keys()):  
            self.xrange =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'npts'
        if fff in list(datadict.keys()):  
            self.npts =int(list(datadict.values())[list(datadict.keys()).index(fff)])
        fff = 'rms'
        if fff in list(datadict.keys()):  
            self.rms =float(list(datadict.values())[list(datadict.keys()).index(fff)])
        fff = 'phase_fit'
        if fff in list(datadict.keys()):  
            self.phase_fit =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'actions'
        if fff in list(datadict.keys()):  
            self.actions =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'remlist'
        if fff in list(datadict.keys()):  
            self.remlist =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'hpeos_version'
        if fff in list(datadict.keys()):  
            self.hpeos_version =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'scipy_version'
        if fff in list(datadict.keys()):  
            self.scipy_version =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'x'
        if fff in list(datadict.keys()):  
            self.x =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'x_std'
        if fff in list(datadict.keys()):  
            self.x_std =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'f'
        if fff in list(datadict.keys()):  
            self.f =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'f_std'
        if fff in list(datadict.keys()):  
            self.f_std =np.asarray(list(datadict.values())[list(datadict.keys()).index(fff)],dtype=float)
        fff = 'refid'
        if fff in list(datadict.keys()):  
            self.refid =list(datadict.values())[list(datadict.keys()).index(fff)]
        fff = 'reflist'
        if fff in list(datadict.keys()):  
            self.reflist =list(datadict.values())[list(datadict.keys()).index(fff)]
        #
@dataclass
class HPeos_ID:
    label: str = ''
    i: int = 0
#
@dataclass
class HPeos_Point:
    """Data Class to hold an EOS data point."""
    up:           float = 0.   # particle velocity
    us:           float = 0.   # shock velocity
    P:            float = 0.   # pressure
    rho:          float = 0.   # density
    S:            float = 0.   # specific entropy
    T:            float = 0.   # temperature
    U:            float = 0.   # specific internal energy
    up_units:  str = ''
    us_units:  str = ''
    P_units:  str = ''
    rho_units:  str = ''
    S_units:  str = ''
    T_units:  str = ''
    U_units:  str = ''
    note: str = ''
    def __iter__(self):
        '''Needed to convert class object into dictionary.'''
        for key in self.__dict__:
            yield key, getattr(self,key)
#
####################################################################
# not a dataclass because most of the possible data fields will not be used
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
        self.filename = []      # Data file name
        self.datapath = []      # Path to data file
        self.matid = []         # Database Material ID number = IHED ID for those in IHED
        self.material = []      # Database Material name
        self.formula = []       # Database Material chemical formula
        self.rho_ref = []       # Database material reference density
        self.rho_ref_std = []   
        self.rho_ref_units = []
        self.phase_ref = []     # short name of phase/crystal structure of reference state
        self.db_keys = ['filename','datapath','matid','material','formula','rho_ref','rho_ref_std','rho_ref_units','phase_ref']
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
        self.orientation = []  # xtl orientation     
        self.init_keys = ['m','rho0','rho0_std','rho0_units','P0','P0_std','P0_units','T0','T0_std','T0_units','phase0','orientation']
        ### these intermediate keys are data fields that may be unused; expect tabular format for included data
        # first shock information (used for decaying shocks; need another category to include time as a variable?)
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
        self.T_source = []        
        self.T_units = []
        self.T_b = []            # alternate temperature; make sure explained in remarks
        self.T_b_std = []          
        self.T_b_units = []
        self.emissivity = [] # emissivity for shock temperature
        self.emissivity_std = [] # emissivity for shock temperature
        self.reflectivity = [] # reflectivity for shock temperature
        self.reflectivity_std = [] # reflectivity for shock temperature error
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
        self.T_rel_source = []    
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
        ##### metadata not in tabular format
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
        # Amodeo set annotations
        self.streak_flag = []
        self.wiggle_flag = []
        self.meta_keys = ['ref','rem','quality','reflist','refidlist','remlist','denotions','groups','collections','actions','datatype','platform','diagnostics','calcinfo','methods','streak_flag','wiggle_flag']
    def calcPrho(self):
        ''' Calculate pressure and density values from particle velocity and shock velocity values. No error propagation. '''
        if len(self.rho0) == len(self.up):
            self.P = self.rho0 * self.us * self.up
            self.rho = 1./( (1./self.rho0) * (1. - self.up/self.us) )
            self.U = 0.5*self.P*1.e9*((1/(1000*self.rho0))-(1/(1000*self.rho))) # J/kg
        else:
            self.P = self.rho_ref * self.us * self.up
            self.rho = 1./( (1./self.rho_ref) * (1. - self.up/self.us) )
            self.U = 0.5*self.P*1.e9*((1/(1000*self.rho_ref))-(1/(1000*self.rho))) # J/kg
    def savejson(self):
        # write json file that can be read back into the dataclass
        with open(self.filename+'.txt', 'w') as f:
            temp = json.dumps(dict(self),cls=NumpyEncoder)
            f.write(temp) # writes as one long string
        print('Saved HPeos_Shock object in '+self.filename+'.txt')
    def savereadable(self):
        # write with new lines after each field for human readability and csv table of data 
        ddd = dict(self)
        nfields = len({field.name for field in fields(self)})
        with open(self.outputfitpath+self.filename+'-readable.txt', 'w') as f:
            #for i in range(len(ddd.keys())):
            for i in range(nfields): # do not write the last 6 fields as key:value pairs
                f.write(str(list(ddd.keys())[i])+': '+str(list(ddd.values())[i])+'\n')
            f.write('\n')
            if 0:
                f.write('TABLE\n')
                #f.write('%16s, %16s, %16s, %16s, %6s\n' % ('x', 'x_std', 'f', 'f_std', 'refid'))
                f.write('%20s, %20s, %20s, %20s, %6s\n' % (self.x_desc, 'stdev', self.f_desc, 'stdev', 'refid'))
                f.write('%20s, %20s, %20s, %20s, %6s\n' % (self.x_units, self.x_units, self.f_units,  self.f_units, '-'))
                for i in range(len(self.x)):
    #                f.write(str(self.x[i])+', '+str(self.x_std[i])+', '+str(self.f[i])+', '+str(self.f_std[i])+', '+str(self.refid[i])+'\n')
                    f.write('%20.8e, %20.8e, %20.8e, %20.8e, %6i\n' % \
                    (self.x[i],self.x_std[i],self.f[i],self.f_std[i],self.refid[i]))
                f.write('\n')
                f.write('REFERENCES\n')
                for item in ddd['reflist']:
                    f.write(item+'\n')
        print('Saved HPeos_Shock object in readable file '+self.filename+'-readable.txt')                
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
        self.p8000 = []
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
        self.t8000 = []
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
        self.st8000 = []
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
        self.sp8000 = []
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
        self.su8000 = []
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


class Dict2Class(object):
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])
    def __iter__(self):
        '''Print the class.'''
        for key in self.__dict__:
            yield key, getattr(self,key)

# End of classes for silica melt curve paper
#####################################
