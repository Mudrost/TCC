# Apêndice A - Propulsion.py

import matplotlib.pyplot as plot
import numpy as np
import Utils as utils
from statistics import mean
from scipy import optimize
import random

def get_burn_rate(pressure):
    a = 0
    b = 0
    if pressure < 103000:
        return 0.0026
    elif pressure >= 103000 and pressure < 807000:
        a = 10.71
        b = 0.625
    elif pressure >= 807000 and pressure < 1500000:
        a = 8.763
        b = -0.314
    elif pressure >= 1500000 and pressure < 3790000:
        a = 7.852
        b = -0.013
    elif pressure >= 3790000 and pressure < 7030000:
        a = 3.907
        b = 0.535
    elif pressure >= 7030000 and pressure < 10670000:
        a = 9.653
        b = 0.064
    else:
        #raise Exception("Pressão na câmara acima do intervalo conhecido")
        # Para evitar problemas com a otimização, manter mesmas constantes, mas verificar se a pressão superou 1550 psi
        a = 9.653
        b = 0.064
    return (a*(pressure/1000000)**b)/1000

def get_bates_burning_area(inner_diam, outer_diam, length,no_grains = 1):
    return no_grains*np.pi*(0.5*(outer_diam**2 - inner_diam**2) + length*inner_diam)

def get_kn(throat_area, surface_area):
    return surface_area/throat_area

def get_chamber_pressure(Kn, propelant_density, burn_rate, c_star):
    return Kn*propelant_density*burn_rate*c_star

def get_optimal_expansion_ratio(chamber_pressure, k, ambient_pressure):
    if chamber_pressure < 1e-8:
         return 0
    return 1/(((k+1)/2)**(1/(k-1)) *(ambient_pressure/chamber_pressure)**(1/k)*np.sqrt((k+1)/(k-1) *(1 - (ambient_pressure/chamber_pressure)**((k-1)/k))))

def get_optimal_exit_area(chamber_pressure, k, ambient_pressure):
    return get_optimal_expansion_ratio(chamber_pressure, k, ambient_pressure)*A_t

def get_exit_pressure(throat_area, exit_area, k, chamber_pressure):
    return optimize.fsolve(lambda x: (throat_area/exit_area) - (((k+1)/2)**(1/(k-1)) *(x/chamber_pressure)**(1/k)*np.sqrt(((k+1)/(k-1)) *(1 - (x/chamber_pressure)**((k-1)/k)))), 0)[0]

def get_coef_thrust(k, exit_pressure, chamber_pressure, exit_area, ambient_pressure, throat_area):
    return np.sqrt((2*k**2)/(k-1) * (2/(k+1))**((k+1)/(k-1))*(1 - (exit_pressure/chamber_pressure)**((k-1)/k))) + (exit_pressure - ambient_pressure)*exit_area/chamber_pressure/throat_area

def get_thrust(thrust_coef, throat_area, chamber_pressure):
    return thrust_coef*throat_area*chamber_pressure

def get_mass_flux(burn_rate, surface_area, d, rho):
    return burn_rate*surface_area*rho/(utils.get_circ_area(d))


# Funções para alterar a geometria dos grãos ou área da garganta
def set_outer_diam(new_D):
    global D
    D = new_D

def set_inner_diam(new_d):
    global d
    d = new_d

def set_length(new_L):
    global L
    L = new_L

def set_no_grains(new_N):
    global N
    N = new_N

def set_throat_diam(new_dstar):
    global A_t
    A_t = np.pi*(new_dstar/2)**2



# Parâmetros do motor
D = 0.03738                                       # Diâmetro externo (m)
d = 0.015                                         # Diâmetro interno (m)
L = 0.05                                          # Comprimento do grão (m)
A_t= np.pi*(0.01/2)**2                            # Área da garganta (m^2)
N = 6                                            #  Número de grãos 


# Parâmetros do propelente
r0 = 0.0026                                     # Taxa de queima em 1 atm (m/s)
c_star_knsb = 885                               # Velocidade característica (m/s)
rho_knsb_ideal = 1750                           # Densidade do propelente (kg/m^3)
k = 1.1361                                      # Coeficiente de expansão adiabática 
knsb = {'r0' : r0, 'c_star' : c_star_knsb, 'rho' : rho_knsb_ideal, 'k' : k}

# Outros parâmetros
ambient_pressure = 101325                       # Pressão a nível do mar (Pa)
propelant_mass = N*np.pi *(D**2 - d**2)*L * knsb['rho']/4


def simulate_nozzle(params):
    global ambient_pressure
    global N
    d, L, d_t = params
    A_t = np.pi*(d_t/2)**2
    cur_diam = d                                                    # Diâmetro atual
    cur_length = L                                                  # Comprimento atual
    cur_pressure = ambient_pressure                                 # Pressão atual
    it = 0                                                          # Iteração
    D = 0.03664
    burning_area = [get_bates_burning_area(cur_diam, D, L, N)]      # Dados de superfície exposta
    kn = [burning_area[0]/A_t]                                      # Dados de Kn
    chamber_pressure = [ambient_pressure]                           # Dados de pressão na câmara
    timestep = 0.02                                                 # Passo de simulação
    mass_flow = [0]
    # Loop de queima
    while (True):
        r = get_burn_rate(cur_pressure)
        cur_diam = cur_diam + 2*r*timestep
        cur_length = cur_length - 2*r*timestep
        if (D - cur_diam) < 1e-8 or cur_length < 1e-8:
            break
        burning_area.append(get_bates_burning_area(cur_diam, D, cur_length, N))
        cur_pressure = get_chamber_pressure(get_kn(A_t,burning_area[-1]), knsb['rho'], r, knsb['c_star'])
        chamber_pressure.append(cur_pressure)
        kn.append(burning_area[-1]/A_t)
        mass_flow.append(get_mass_flux(r, burning_area[-1], cur_diam, knsb['rho']))
        it = it + 1

    # Listas de plot
    time =     np.linspace(0, it*timestep, it + 1)
    pressure = np.linspace(0, 1e7, 50)

    # Cálculo da pressão média e razão de expansão ideal
    chamber_avg = mean(chamber_pressure)
    ae_at = get_optimal_expansion_ratio(chamber_avg, knsb['k'], ambient_pressure)
    ae = ae_at * A_t

    # Cálculo da pressão de saída, coeficiente de empuxo e empuxo ao longo da simulação
    pe = []
    thrust_coef = []
    thrust = []
    for i in chamber_pressure:
        pe.append(get_exit_pressure(A_t, ae, knsb['k'],i))
        thrust_coef.append(get_coef_thrust(knsb['k'], pe[-1], i, ae, ambient_pressure, A_t))
        thrust.append(get_thrust(thrust_coef[-1], A_t, i))
   
    burn_time = timestep*it
    total_impulse = mean(thrust)*timestep*it


    print("Tempo de queima: ", burn_time,"s")
    print("Razão de expansão ideal: ",ae_at)
    print("Diâmetro de saída: ", utils.get_circ_diam(ae)*1e3,"mm")
    print("Coeficiente de empuxo ideal: ", get_coef_thrust(knsb['k'], ambient_pressure, chamber_avg, ae, ambient_pressure, A_t))
    print("Impulso total: ", total_impulse,"Ns" )
    print("Fluxo de massa máx.: ", max(mass_flow),"Ns" )
    print("--------")
    return burning_area, chamber_pressure, kn, pe, thrust_coef, thrust, burn_time, total_impulse, mass_flow, time, pressure


## Funções de otimização
def optimize_total_impulse(params):
    _, _, _, _, _, _, _, total_impulse,_, _, _= simulate_nozzle(params)
    return abs(total_impulse - 620)

def restraint_mean_pressure(params):
    _, chamber_pressure, _, _,_,_,_,_,_,_,_ = simulate_nozzle(params)
    return 1.5e6 - abs(8.03e6 - mean(chamber_pressure))

def restraint_peak_pressure(params):
    _, chamber_pressure, _, _,_,_,_,_,_,_,_ = simulate_nozzle(params)
    return 8.03e6 - max(chamber_pressure)

def restraint_kn_variation(params):
    _, _, kn, _,_,_,_,_,_,_,_ = simulate_nozzle(params)
    return 0.2 - (max(kn) - min(kn))/(max(kn))

def restraint_port_ratio(params):
    d, L, d_t = params
    port_area = np.pi*d**2/4
    throat_area = np.pi*d_t**2/4
    return (port_area/throat_area - 2)

def restraint_mass_flux(params):
    _, _, _, _,_,_,_,_,mass_flux,_,_ = simulate_nozzle(params)
    return 1300 - max(mass_flux)


# Otimiza um motor com um impulso total como objetivo e restrições
def optimize_for_impulse():
    constraint = [{'type' : 'ineq', 'fun' : restraint_mean_pressure},                 # Restrição de design #1 pressão média prox de 8.55e6 Pa
                  {'type' : 'ineq', 'fun' : restraint_port_ratio},                    # Restrição de design #2 Port/Throat > 2 
                  {'type' : 'ineq', 'fun' : restraint_peak_pressure},                 # Pressão máx. < 11.37e6 Pa
                  {'type' : 'ineq', 'fun' : restraint_kn_variation},                  # Restrição de Kn no máx 15%
                  {'type' : 'ineq', 'fun' : restraint_mass_flux}]                     # Restrição de design #3 fluxo de massa máx. 2 lb/sqin/s
    bounds = ((0.01, 0.03664), (0.1, 0.225), (0.005, 0.015))
    it = 0
    while (True):
        initial_guess = [random.uniform(0.01, 0.03664), random.uniform (0.1, 0.225), random.uniform(0.005, 0.015)]
        a = optimize.minimize(optimize_total_impulse, initial_guess, bounds=bounds, constraints=constraint ) 
        if a.success and a.fun < 50:
            print(initial_guess)
            a
            break
        elif it > 1000:
            global N
            N = N+1
        it = it + 1
        print(it)
    
optimize_for_impulse()


# Taxa de queima vs pressão (KNSB)
f= plot.figure()
ax = f.gca()
#plot.title('Burn rate vs pressure')
plot.plot(pressure, list(map(get_burn_rate, pressure)))
plot.xlabel("Pressure (Pa)")
plot.ylabel("Burn rate (m/s)")
plot.grid()
 

# Kn vs tempo
f= plot.figure()
ax = f.gca()
plot.plot(time, kn)
ax.set_ylim(ymin=0)
plot.xlabel("Tempo (s)")
plot.ylabel("$K_n$")
plot.title("$K_n$ vs time")
plot.grid()

# Pressão vs tempo
f= plot.figure()
ax = f.gca()
plot.plot(time, chamber_pressure)
ax.set_ylim(ymin=0)
plot.xlabel("Tempo (s)")
plot.ylabel("Pressure (Pa)")
plot.title("Pressure vs time")
plot.grid()


# Razão de expansão ideal vs tempo
f= plot.figure()
ax = f.gca()
plot.plot(time, list(map(lambda n: get_optimal_expansion_ratio(n,k =k, ambient_pressure = ambient_pressure), chamber_pressure)))
ax.set_ylim(ymin=0)
plot.xlabel("Tempo (s)")
plot.ylabel("Ae/A*")
#plot.title("Pressure vs time")
plot.grid()

f= plot.figure()
ax = f.gca()
plot.plot(time, thrust_coef)
ax.set_ylim(ymin=0)
plot.xlabel("Tempo (s)")
plot.ylabel("Ae/A*")
#plot.title("Pressure vs time")
plot.grid()


f= plot.figure()
ax = f.gca()
plot.plot(time, thrust)
ax.set_ylim(ymin=0)
plot.xlabel("Tempo (s)")
plot.ylabel("Ae/A*")
#plot.title("Pressure vs time")
plot.grid()
plot.show()