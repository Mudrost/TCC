import matplotlib.pyplot as plot
import numpy as np
import Utils as utils
from statistics import mean
from scipy import optimize
import random




def get_reynolds(V, L, h):
    isa = utils.get_real_atmosphere(h)
    rho = isa['rho']
    mu = isa['mu']
    return rho*V*L/mu

def get_friction(reynolds, surface_roughness, rocket_ref_length):
    Re_crit = 51*(surface_roughness/rocket_ref_length)**(-1.039)
    a = 1.0/(((1.5*np.log(reynolds))-5.6)**2)
    b = 0.032*(surface_roughness/rocket_ref_length)**0.2
    return np.maximum(a,b)
    if reynolds < 1e4:
        return 1.48e-2
    elif reynolds >= 1e4 and reynolds < Re_crit:
        return 1.0/(((1.5*np.log(reynolds))-5.6)**2)
    else:
        return 0.032*(surface_roughness/rocket_ref_length)**0.2

def get_friction_compressed(friction_coefficient, mach):                       
    return friction_coefficient*(1-0.1*mach**2)

def get_scaled_friction_coef(friction_coefficient_compressed, param, wet_area, ref_area):
    return friction_coefficient_compressed*(1+ param)*wet_area/ref_area

def get_von_karman_nosecone_friction_drag(length, base_diam, rocket_length, rocket_ref_diam ,mach = 0.3, h = 0, surface_roughness = 20e-6):
    # Arrasto de fricção de uma coifa de Von Karman, com base no seu comprimento e diâmetro.
    # Integração numérica para área molhada:
    n = 5000
    dx = length/n
    x = np.linspace(0, length, n)
    theta = np.arccos(1 - 2*x/length)
    y = base_diam/2 * np.sqrt(theta - np.sin(2*theta)/2)/np.sqrt(np.pi)
    wet_area = 0
    for i,k in enumerate(x, start = 0):
        if i == 0:
            continue
        wet_area += np.pi*(y[i-1]+y[i])*np.sqrt((y[i-1]-y[i])**2 + dx**2)
    #
    a = utils.get_real_atmosphere(h)['a']
    v = mach*a
    Re = get_reynolds(v, rocket_length, h)
    fineness_ratio = rocket_length/(base_diam/2)

    coef1 = get_friction(Re, surface_roughness, rocket_length)
    coef2 = get_friction_compressed(coef1, mach)
    param = 1/(2*fineness_ratio)
    coef3 = get_scaled_friction_coef(coef2, param, wet_area, utils.get_circ_area(rocket_ref_diam))
    print(coef3)

get_von_karman_nosecone_friction_drag(0.1, 0.068, 0.67, 0.068, 0.3, surface_roughness=500e-6)

def get_transition_friction_drag(length, diam_1, diam_2, rocket_length, rocket_ref_diam ,mach = 0.3, h = 0, surface_roughness = 20e-6):
    # Arrasto de fricção para uma transição, com base nos dois diâmetros.
    # Área molhada de um tronco cônico:
    wet_area = np.pi * (diam_1/2 + diam_2/2)*np.sqrt((diam_1/2 - diam_2/2)**2 + length**2)
    #
    a = utils.get_real_atmosphere(h)['a']
    v = mach*a
    Re = get_reynolds(v, rocket_length, h)
    fineness_ratio = rocket_length/(rocket_ref_diam/2)
    coef1 = get_friction(Re, surface_roughness, rocket_length)
    coef2 = get_friction_compressed(coef1, mach)
    param = 1/(2*fineness_ratio)
    coef3 = get_scaled_friction_coef(coef2, param, wet_area, utils.get_circ_area(rocket_ref_diam))
    print(coef3)

def get_fins_friction_drag(chord_root, chord_tip, span, thickness, rocket_length ,rocket_ref_diam, mach = 0.3, no_fins = 4, h = 0, surface_roughness = 20e-6):
    # Arrasto de fricção para o conjunto de aletas
    # Área molhada de uma aleta trapezoidal
    taper = chord_tip/chord_root
    mac = 2/3 *chord_root*(taper**2 + taper + 1)/(taper + 1)
    wet_area = no_fins * 2 * mac * span
    
    a = utils.get_real_atmosphere(h)['a']
    v = mach*a
    Re = get_reynolds(v, rocket_ref_diam, h)

    coef1 = get_friction(Re, surface_roughness, rocket_length)
    coef2 = get_friction_compressed(coef1, mach)
    param = 2*thickness/mac
    coef3 = get_scaled_friction_coef(coef2, param, wet_area, utils.get_circ_area(rocket_ref_diam))
    print(coef3)


get_transition_friction_drag(0.05, 0.042, 0.068, 0.67, 0.068, 0.3, surface_roughness=500e-6)

get_fins_friction_drag(0.07, 0.06, 0.05, 0.003, 0.69, 0.068,  surface_roughness=500e-6)