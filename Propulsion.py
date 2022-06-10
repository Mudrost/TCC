# Apêndice A - Propulsion.py

import matplotlib.pyplot as plot
import numpy as np
import Utils as utils
from statistics import mean

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
        raise Exception("Pressão na câmara acima do intervalo conhecido")
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





# Parâmetros do motor
D = 0.038                                       # Diâmetro externo (m)
d = 0.015                                       # Diâmetro interno (m)
L = 0.05                                       # Comprimento do grão (m)
A_t= np.pi*(0.008/2)**2                         # Área da garganta (m^2)
N = 4                                           # Número de grãos 

# Parâmetros do propelente
r0 = 0.0026                                     # Taxa de queima em 1 atm (m/s)
c_star_knsb = 885                               # Velocidade característica (m/s)
rho_knsb_ideal = 1750                           # Densidade do propelente (kg/m^3)
k = 1.1361                                      # Coeficiente de expansão adiabática 
# Outros parâmetros
ambient_pressure = 101325                       # Pressão a nível do mar (Pa)


# Variáveis internas
cur_diam = d                                                    # Diâmetro atual
cur_length = L                                                  # Comprimento atual
cur_pressure = ambient_pressure                                 # Pressão atual
it = 0                                                          # Iteração
burning_area = [get_bates_burning_area(cur_diam, D, L, N)]      # Dados de superfície exposta
kn = [burning_area[0]/A_t]                                      # Dados de Kn
chamber_pressure = [ambient_pressure]                           # Dados de pressão na câmara
timestep = 0.001                                                # Passo de simulação

# Loop
while (True):
    r = get_burn_rate(cur_pressure)
    cur_diam = cur_diam + 2*r*timestep
    cur_length = cur_length - 2*r*timestep
    if (D - cur_diam) < 1e-8 or cur_length < 1e-8:
        break
    burning_area.append(get_bates_burning_area(cur_diam, D, cur_length, N))
    cur_pressure = get_chamber_pressure(get_kn(A_t,burning_area[-1]), rho_knsb_ideal, r, c_star_knsb)
    chamber_pressure.append(cur_pressure)
    kn.append(burning_area[-1]/A_t)
    it = it + 1

# Condição de contorno:
kn[-1] = 0

# Listas de plot
time =     np.linspace(0, it*timestep, it + 1)
pressure = np.linspace(0, 1e7, 50)

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
plot.show()

chamber_avg = mean(chamber_pressure)
ae_at = get_optimal_expansion_ratio(chamber_avg, k, ambient_pressure)
print(utils.get_circ_diam(ae_at*A_t))