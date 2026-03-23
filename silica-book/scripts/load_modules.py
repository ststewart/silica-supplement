# import the same modules for every notebook

import os
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import pickle
# smoothing function
from scipy.signal import savgol_filter
from scipy import linalg

# Get date and time values in ISO 8601 format
# importing datetime and time module
import datetime
import zoneinfo

# used to document curve fits to data
import json
from dataclasses import dataclass,field, fields

# optional for main Jupyterbook
#from pprint import pprint
#import pandas as pd # pandas is used in impedance match subfolders only

# for customizing colors
from cycler import cycler

# customized python classes and functions for EOS data
import hpeos_module as hp # holds shock Hugoniot python class
import hpeos_mc_module as hpmc # holds shock Hugoniot python classes for montecarlo calculations


# Calling the isoformat() function over the current time
# include time zone:
#DateTime_in_ISOFormat = datetime.datetime.now(zoneinfo.ZoneInfo('UTC')).isoformat('T', 'seconds') #

# if changes are made to hpeos_module; use the following to reload the module without restarting the kernel
#from importlib import reload #import function "reload"
#reload(hp)
#import hpeos_module as hp # holds shock Hugoniot python class

# if need to load json EOS fit file at the top level need this class:
#class NumpyEncoder(json.JSONEncoder):
#    def default(self, obj):
#        if isinstance(obj, np.ndarray):
#            return obj.tolist()
#        return super().default(obj)
