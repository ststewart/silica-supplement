# fitting
from scipy.odr import ODR, Model, RealData

def linear_function(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]
linear_model = Model(linear_function)

def quadratic_function(B, x):
    '''Quadratic function y = a*x*x + b*x + c'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x*x + B[1]*x + B[2]
quadratic_model = Model(quadratic_function)

def cubic_function(B, x):
    '''Cubic function y =  a*x*x*x + b*x*x + c*x + d'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x*x*x + B[1]*x*x + B[2]*x + B[3]
cubic_model = Model(cubic_function)
