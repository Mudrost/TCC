import Aerodynamics as aero
import Propulsion as prop
import Utils as utils
import matplotlib.pyplot as plot
import numpy as np

# Parâmetros do foguete
m_empty = 0.9176 + .350 + .200 + .300 + .200 - .299
m_propellant = .299
m_full = m_empty + m_propellant
d_ref = 0.07
area_ref = utils.get_circ_area(d_ref)
burn_time = 0.79

# Outros parâmetros
g = 9.81
current_height = 0
current_mass = m_full
average_mdot = m_propellant/burn_time
burning_area, chamber_pressure, kn, pe, thrust_coef, thrust, burn_time, total_impulse, mass_flow, time, pressure, burn_rate = prop.simulate_motor([0.01859, 0.05463, 0.007560])

valores_Cd = [aero.get_rocket_cd(0)]
valores_V = [0]
valores_S = [0]

def get_thrust_at_time(t):
    if t == 0 or t>=max(time):
        return 0
    else:
        # Interpolação linear
        max_t = time[time > t].min()
        min_t = time[time < t].max()
        d = time.tolist()
        return thrust[d.index(min_t)] + (t - min_t)*(thrust[d.index(max_t)] - thrust[d.index(min_t)])



def drag_dv_dt(t, v):
    if v < 0:
        return 0
    thrust = get_thrust_at_time(t)
    rho = utils.get_real_atmosphere(valores_S[-1])['rho']
    a = utils.get_real_atmosphere(valores_S[-1])['a']
    mach = v/a
    drag_coef = aero.get_rocket_cd(mach)
    if t <= burn_time:
        accel = -g + thrust/current_mass -0.5*v**2*drag_coef*area_ref*rho/current_mass
        return accel
    elif t > burn_time:
        return -g -0.5*v**2*drag_coef*area_ref*rho/current_mass

def rungeKutta(t0, V0, h):
    V = V0
    k1 = h * drag_dv_dt(t0, V)
    k2 = h * drag_dv_dt((t0 + h/2), (V + k1/2))
    k3 = h * drag_dv_dt((t0 + h/2), (V + k2/2))
    k4 = h * drag_dv_dt(t0 + h, V + k3)
    V = V + (k1 + 2*k2 + 2*k3 + k4)/6
    if V < 0:
        return 0
    valores_Cd.append(aero.get_rocket_cd(V/utils.get_real_atmosphere(valores_S[-1])['a']))
    return V


def solve_trajectory():
    global h, current_height, current_mass
    stopFlag = False
    i = 0
    timestep = 0.01
    while stopFlag is False:
        if valores_V[-1] < 0.01 and i*timestep > burn_time:
            stopFlag = True
        valores_V.append(rungeKutta(i*timestep, valores_V[-1],timestep))
        valores_S.append(valores_S[-1]+ valores_V[-1]*timestep)
        if i*timestep < burn_time:
            current_mass -= average_mdot*timestep
        i=i+1
    valores_t = np.linspace(0, i*timestep, i+1)
    return valores_t
    

# Plotar curvas do OpenRocket juntas:
def compare_curves():
    valores_t = solve_trajectory()
    OPEN_ROCKET_DATA = []
    OPEN_ROCKET = open('openrocket_sim.csv', 'r')

    for line in OPEN_ROCKET:
        OPEN_ROCKET_DATA.append(np.array(line.split(), dtype=np.float32))

    f= plot.figure()
    ax = f.gca()
    plot.plot(list(list(zip(*OPEN_ROCKET_DATA))[0]), list(list(zip(*OPEN_ROCKET_DATA))[1]), '--g',label = 'OpenRocket')
    plot.plot(valores_t, valores_S, '-.r', label = 'Apêndice D')
    ax.set_ylim(ymin=0)
    ax.set_xlim(xmax=max(valores_t), xmin = 0)
    plot.xlabel("Tempo (s)")
    plot.ylabel("Altitude (m)")
   
    plot.grid()
    xmax = valores_t[np.argmax(valores_S)]
    ymax = max(valores_S)
    text= "t={:.3f}, h={:.3f}".format(xmax, ymax)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(0.9,0.7), **kw)

    xmax = list(list(zip(*OPEN_ROCKET_DATA))[0])[np.argmax(list(list(zip(*OPEN_ROCKET_DATA))[1]))]
    ymax = max(list(list(zip(*OPEN_ROCKET_DATA))[1]))
    text= "t={:.3f}, h={:.3f}".format(xmax, ymax)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=5,angleB=-60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(0.65,0.95), **kw)
    plot.legend(loc = 'best')
    plot.show()

    f= plot.figure()
    ax = f.gca()
    plot.plot(list(list(zip(*OPEN_ROCKET_DATA))[0]), list(list(zip(*OPEN_ROCKET_DATA))[2]), '--g',label = 'OpenRocket')
    plot.plot(valores_t, valores_V, '-.r', label = 'Apêndice D')
    plot.xlabel("Tempo (s)")
    plot.ylabel("Velocidade vertical (m/s)")
    xmax = valores_t[np.argmax(valores_V)]
    ymax = max(valores_V)
    text= "t={:.3f}, h={:.3f}".format(xmax, ymax)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=170,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(0.65,0.9), **kw)

    ax.set_ylim(ymin=0)
    ax.set_xlim(xmax=max(valores_t), xmin = 0)
    plot.legend(loc = 'best')
    plot.grid()
    plot.show()

compare_curves()