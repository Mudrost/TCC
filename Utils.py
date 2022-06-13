
import numpy as np

def get_circ_area(diam):
    return np.pi*(diam/2)**2

def get_circ_diam(area):
    return 2*((area/np.pi)**0.5)

def get_real_atmosphere(h, deltaT = 0):
    H_0 = 0
    p_0 = 101325
    rho_0= 1.225
    T_0=  288.15
    g_0 = -9.80665
    a0 = 340.294
    R = 287.058
    L = -0.0065
    gamma = 1.4
    mu_0 = 1.783e-5
    T_tropopausa = T_0 + 11000*L

    # Calcular a temperatura real:
    T_ISA = 0
    if h <= 11000:
       T_ISA = T_0 + L*h
    elif 11000 < h <= 20000:
       T_ISA = T_tropopausa
    T_real = T_ISA + deltaT


    # Calcular pressão real: igual à variação da pressão ISA (U.S. Standard Atmosphere 1976)
    if h <= 11000:
        p_p0 =(T_ISA/T_0)**(g_0/(L*R))
        p_real = p_p0 * p_0
    elif h > 11000 and h <= 20000:
        p_p0 = ((T_tropopausa)/(T_0))**(g_0/(L*R))
        p_real = p_p0* p_0 * np.exp((11000-h)/6342)

    # Calcular densidade real:
    rho_real = p_real/T_real/R
    
    # Calcular velocidade local do som:
    a_real = np.sqrt(gamma*R*T_real)

    # Calcular visocidade real:
    S = 113 # para o ar seco
    mu_mu0 = (T_real/T_0)**1.5 * ((T_0 + S)/(T_real + S))
    mu_real = mu_mu0 * mu_0
    return {'T': T_real, 'rho': rho_real, 'p': p_real, 'mu': mu_real, 'a': a_real}
