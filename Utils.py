
import numpy as np

def get_circ_area(diam):
    return np.pi*(diam/2)**2

def get_circ_diam(area):
    return 2*((area/np.pi)**0.5)