import matplotlib.pyplot as plot
import numpy as np
import Utils as utils
from statistics import mean
from scipy import optimize
import random

def get_burn_rate(pressure, is_optimizing = False):
    a = 0
    b = 0
    if pressure <= 103000:
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
    elif is_optimizing:
        a = 9.653
        b = 0.064
    else:
        raise Exception("Pressão na câmara acima do intervalo conhecido")
    return (a*(pressure/1000000)**b)/1000

def get_bates_burning_area(inner_diam, outer_diam, length,no_grains = 1):
    return no_grains*np.pi*(0.5*(outer_diam**2 - inner_diam**2) + length*inner_diam)

def get_kn(throat_area, surface_area):
    return surface_area/throat_area

def get_chamber_pressure(Kn, propelant_density, burn_rate, c_star):
    return Kn*propelant_density*burn_rate*c_star

def get_optimal_expansion_ratio(chamber_pressure, k, ambient_pressure):
    if chamber_pressure < 1e-8 or chamber_pressure == ambient_pressure:
         return 0
    return 1/(((k+1)/2)**(1/(k-1)) *(ambient_pressure/chamber_pressure)**(1/k)*np.sqrt((k+1)/(k-1) *(1 - (ambient_pressure/chamber_pressure)**((k-1)/k))))


def get_optimal_exit_area(chamber_pressure, k, ambient_pressure):
    return get_optimal_expansion_ratio(chamber_pressure, k, ambient_pressure)*A_t

def get_exit_pressure(throat_area, exit_area, k, chamber_pressure):
    if exit_area < 1e-8:
        return 0
    return optimize.fsolve(lambda x: (throat_area/exit_area) - (((k+1)/2)**(1/(k-1)) *(x/chamber_pressure)**(1/k)*np.sqrt(((k+1)/(k-1)) *(1 - (x/chamber_pressure)**((k-1)/k)))), 0)[0]

def get_coef_thrust(k, exit_pressure, chamber_pressure, exit_area, ambient_pressure, throat_area):
    return np.sqrt((2*k**2)/(k-1) * (2/(k+1))**((k+1)/(k-1))*(1 - (exit_pressure/chamber_pressure)**((k-1)/k))) + (exit_pressure - ambient_pressure)*exit_area/chamber_pressure/throat_area

def get_thrust(thrust_coef, throat_area, chamber_pressure):
    return thrust_coef*throat_area*chamber_pressure

def get_mass_flux(burn_rate, surface_area, d, rho):
    return burn_rate*surface_area*rho/(utils.get_circ_area(d))

def get_propellant_mass(D, d, L, N, rho):
    return N*np.pi *(D**2 - d**2)*L *rho/4


# Parâmetros do motor
N = 4                                                               # Número de grãos 
D = 0.03664                                                         # Diâmetro externo no grão, fixo em projeto

# Parâmetros do propelente
knsb = {'r0' : 0.0026, 'c_star' : 885, 'rho' : 1750, 'k' : 1.1361}

# Outros parâmetros
ambient_pressure = 101325                       # Pressão a nível do mar (Pa)

def simulate_motor(params, is_optimizing = False):
    cur_diam, cur_length, d_t = params
    prop_mass = get_propellant_mass(D, cur_diam, cur_length, N, knsb['rho'])
    A_t = np.pi*(d_t/2)**2                                                      # Calcular área da garganta partindo do diâmetro
    it = 0                                                                      # Iteração
    burning_area = [get_bates_burning_area(cur_diam, D, cur_length, N)]         # Dados de superfície exposta
    kn = [burning_area[0]/A_t]                                                  # Dados de Kn
    chamber_pressure = [ambient_pressure]                                       # Dados de pressão na câmara
    timestep = 0.001                                                            # Passo de simulação
    mass_flow = [0]
    burn_rate = [knsb['r0']]
    # Loop de queima
    while (True):
        cur_diam += 2*burn_rate[-1]*timestep
        cur_length -= 2*burn_rate[-1]*timestep
        if (D - cur_diam) < 1e-8 or cur_length < 1e-8:
            break
        burning_area.append(get_bates_burning_area(cur_diam, D, cur_length, N))
        kn.append(get_kn(A_t,burning_area[-1]))
        chamber_pressure.append(get_chamber_pressure(kn[-1], knsb['rho'], burn_rate[-1], knsb['c_star']))
        mass_flow.append(get_mass_flux(burn_rate[-1], burning_area[-1], cur_diam, knsb['rho']))
        burn_rate.append(get_burn_rate(chamber_pressure[-1], is_optimizing))
        it = it + 1

    # Listas de plot
    time =     np.linspace(0, it*timestep, it + 1)
    pressure = np.linspace(103000, 1e7, 500)

    # Cálculo da pressão média, razão de expansão ideal e área de saída
    chamber_avg = mean(chamber_pressure)
    ae_at = get_optimal_expansion_ratio(chamber_avg, knsb['k'], ambient_pressure)
    ae = ae_at * A_t

    # Cálculo da pressão de saída, coeficiente de empuxo e empuxo ao longo da simulação
    pe = []
    thrust_coef = []
    thrust = []
    for i,k in enumerate(chamber_pressure):
        if i == 0:
            # Condições de contorno
            pe.append(ambient_pressure)
            thrust_coef.append(0)
            thrust.append(0)
            continue
        pe.append(get_exit_pressure(A_t, ae, knsb['k'],k))
        thrust_coef.append(get_coef_thrust(knsb['k'], pe[-1], k, ae, ambient_pressure, A_t))
        thrust.append(get_thrust(thrust_coef[-1], A_t, k))
   
    burn_time = timestep*it
    total_impulse = mean(thrust)*timestep*it

    print("RESULTADOS DA SIMULAÇÃO: ")
    print("Tempo de queima: ", burn_time,"s")
    print("Razão de expansão ideal: ",ae_at)
    print("Diâmetro de saída: ", utils.get_circ_diam(ae)*1e3,"mm")
    print("Diâmetro da garganta: ", d_t*1e3, "mm")
    print("Coeficiente de empuxo ideal: ", get_coef_thrust(knsb['k'], ambient_pressure, chamber_avg, ae, ambient_pressure, A_t))
    print("Impulso total: ", total_impulse,"Ns" )
    print("Fluxo de massa máx.: ", max(mass_flow),"kg/(s*m^2)" )
    print("MEOP: ", np.max(chamber_pressure)/1e6, "MPa")
    print("Pressão média: ", chamber_avg/1e6, "MPa")
    print("Massa de propelente: ", prop_mass, "g")
    print("------------------------------------------------")
    return burning_area, chamber_pressure, kn, pe, thrust_coef, thrust, burn_time, total_impulse, mass_flow, time, pressure, burn_rate


## Funções de otimização
def optimize_total_impulse(params):
    _, _, _, _, _, _, _, total_impulse,_, _, _, _= simulate_motor(params, True)
    impulse_goal = 434
    return abs(total_impulse - impulse_goal)

def restraint_mean_pressure(params):
    _, chamber_pressure, _, _,_,_,_,_,_,_,_, _ = simulate_motor(params, True)
    mean_chamber_pressure_goal = 8.03e6
    max_deviation = 1e6
    return max_deviation - abs(mean_chamber_pressure_goal - mean(chamber_pressure))

def restraint_peak_pressure(params):
    _, chamber_pressure, _, _,_,_,_,_,_,_,_,_ = simulate_motor(params, True)
    max_chamber_pressure_goal = 8.03e6
    return max_chamber_pressure_goal - max(chamber_pressure)

def restraint_kn_variation(params):
    _, _, kn, _,_,_,_,_,_,_,_,_ = simulate_motor(params, True)
    max_kn_variation_goal = 0.15
    return max_kn_variation_goal - (max(kn) - min(kn))/(max(kn))

def restraint_port_ratio(params):
    d, _, d_t, = params
    port_area = np.pi*d**2/4
    throat_area = np.pi*d_t**2/4
    min_port_throat_ratio_goal = 2
    return (port_area/throat_area - min_port_throat_ratio_goal)

def restraint_mass_flux(params):
    _, _, _, _,_,_,_,_,mass_flux,_,_,_ = simulate_motor(params,True)
    max_mass_flux_goal = 1406
    return max_mass_flux_goal - max(mass_flux)


# Otimiza um motor com um impulso total como objetivo e restrições
def optimize_for_impulse():
    constraint = [{'type' : 'ineq', 'fun' : restraint_mean_pressure},                 # Restrição de design #1 pressão média prox de 8.03e6 Pa
                  {'type' : 'ineq', 'fun' : restraint_port_ratio},                    # Restrição de design #2 Port/Throat > 2 
                  {'type' : 'ineq', 'fun' : restraint_peak_pressure},                 # Pressão máx. < 8.03e6  Pa
                  {'type' : 'ineq', 'fun' : restraint_kn_variation},                  # Restrição de Kn no máx 15%
                  {'type' : 'ineq', 'fun' : restraint_mass_flux}]                     # Restrição de design #3 fluxo de massa máx. 2 lb/sqin/s
    bounds = ((0.005, 0.03664*0.6), (0.05, 0.225), (0.003, 0.015))                    # Limites dimensionais
    it = 0
    while (True):
        initial_guess = [random.uniform(0.005, 0.03664*0.6), random.uniform (0.1, 0.225), random.uniform(0.005, 0.015)]
        a = optimize.minimize(optimize_total_impulse, initial_guess, bounds=bounds, constraints=constraint ) 
        if a.success and a.fun < 10:
            return a.x, it
        elif it > 500:
            print("Nenhuma configuração adequada foi encontrada em ", it, "iterações.")
            break
        it += 1


data, it = optimize_for_impulse()

# Plotar curvas do METEOR, OpenMotor juntas:
def compare_curves():
    METEOR_DATA =[]
    OPEN_MOTOR_DATA = []
    TCC_DATA = []

    METEOR = open('meteor.txt', 'r')
    OPEN_MOTOR = open('openMotor.txt', 'r')
    TCC = open('motor.txt', 'r')

    for line in METEOR:
        METEOR_DATA.append(np.array(line.split(), dtype=np.float32))
    for line in OPEN_MOTOR:
        OPEN_MOTOR_DATA.append(np.array(line.split(), dtype=np.float32))
    for line in TCC:
        TCC_DATA.append(np.array(line.split(), dtype=np.float32))
    f= plot.figure()
    ax = f.gca()
    plot.plot(list(list(zip(*METEOR_DATA))[0]), list(list(zip(*METEOR_DATA))[1]), '--g',label = 'METEOR')
    plot.plot(list(list(zip(*OPEN_MOTOR_DATA))[0]), list(list(zip(*OPEN_MOTOR_DATA))[1]), '-.r', label = 'openMotor')
    plot.plot(list(list(zip(*TCC_DATA))[0]), list(list(zip(*TCC_DATA))[1]), label = 'Apêndice A')
    ax.set_ylim(ymin=0)
    plot.xlabel("Tempo (s)")
    plot.ylabel("Empuxo (N)")
    plot.legend(loc = 'best')
    plot.grid()
    plot.show()

if __name__ == "__main__":
    burning_area, chamber_pressure, kn, pe, thrust_coef, thrust, burn_time, total_impulse, mass_flow, time, pressure, burn_rate = simulate_motor([0.01859, 0.05463, 0.007560])

    # Taxa de queima vs pressão (KNSB)
    f= plot.figure()
    ax = f.gca()
    plot.title('Taxa de queima x pressão')
    plot.plot(pressure, list(map(get_burn_rate, pressure)))
    plot.xlabel("Pressure (Pa)")
    plot.ylabel("Burn rate (m/s)")
    plot.grid()
    

    # Kn vs tempo
    f= plot.figure()
    ax = f.gca()
    plot.title('Kn x tempo')
    plot.plot(time, kn)
    ax.set_ylim(ymin=0)
    plot.xlabel("Tempo (s)")
    plot.ylabel("$K_n$")
    plot.grid()


    f= plot.figure()
    ax = f.gca()
    plot.plot(time, burn_rate)
    ax.set_ylim(ymin=0)
    plot.xlabel("Tempo (s)")
    plot.ylabel("Taxa de queima (m/s)")
    plot.title("Taxa de queima x tempo")
    plot.grid()

    f= plot.figure()
    ax = f.gca()
    plot.plot(time, chamber_pressure)
    ax.set_ylim(ymin=0)
    plot.xlabel("Tempo (s)")
    plot.ylabel("Pressão (Pa)")
    plot.title("Pressão vs tempo")
    plot.grid()

    f= plot.figure()
    ax = f.gca()
    plot.plot(time, thrust)
    ax.set_ylim(ymin=0)
    plot.xlabel("Tempo (s)")
    plot.ylabel("Empuxo (N)")
    plot.title("Empuxo x tempo ")
    plot.grid()
    plot.show()
    # Salvar curva de empuxo-tempo
    file = open('motor.txt', 'w')
    for i, time in enumerate(time):
        file.write("\t%s\t%s\n" % (time, thrust[i]))
    file.close()
    
    compare_curves()


