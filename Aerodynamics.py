import numpy as np
import Utils as utils

##### Definição das funções de arrasto
def get_reynolds(V, L, h):
    isa = utils.get_real_atmosphere(h)
    rho = isa['rho']
    mu = isa['mu']
    return rho*V*L/mu

def get_friction(reynolds, surface_roughness, rocket_ref_length):
    Re_crit = 51*(surface_roughness/rocket_ref_length)**(-1.039)
    #a = 1.0/(((1.5*np.log(reynolds))-5.6)**2)
    #b = 0.032*(surface_roughness/rocket_ref_length)**0.2
    #return np.maximum(a,b)
    if reynolds < 1e4:
        return 1.48e-2
    elif reynolds >= 1e4 and reynolds < Re_crit:
        return 1.0/(((1.5*np.log(reynolds))-5.6)**2)
    else:
        return 0.032*(surface_roughness/rocket_ref_length)**0.2

def get_friction_compressed(friction_coefficient, mach):                       
    return friction_coefficient*(1-0.1*mach**2)

def get_scaled_friction_coef(friction_coefficient_compressed, param, wet_area, ref_area):
    return friction_coefficient_compressed*((1+ param)*wet_area)/ref_area

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
    return coef3


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
    return coef3

def get_fins_friction_drag(chord_root, chord_tip, span, thickness, rocket_length ,rocket_ref_diam, mach = 0.3, no_fins = 1, h = 0, surface_roughness = 20e-6):
    # Arrasto de fricção para o conjunto de aletas
    # Área molhada de uma aleta trapezoidal
    taper = chord_tip/chord_root
    mac = 2/3 *chord_root*(taper**2 + taper + 1)/(taper + 1)
    mac = chord_tip/2 + chord_root/2
    wet_area = no_fins * 2 * mac * span
    
    a = utils.get_real_atmosphere(h)['a']
    v = mach*a
    Re = get_reynolds(v, rocket_length, h)

    coef1 = get_friction(Re, surface_roughness, rocket_length)
    coef2 = get_friction_compressed(coef1, mach)
    param = 2*thickness/mac
    coef3 = get_scaled_friction_coef(coef2, param, wet_area, utils.get_circ_area(rocket_ref_diam))
    return coef3

def get_body_friction_drag(length, diam, rocket_length, rocket_ref_diam ,mach = 0.3, h = 0, surface_roughness = 20e-6):
    # Arrasto de fricção para uma transição, com base nos dois diâmetros.
    # Área molhada de um tronco cônico:
    wet_area = 2 * length *np.pi * diam/2
    #
    a = utils.get_real_atmosphere(h)['a']
    v = mach*a
    Re = get_reynolds(v, rocket_length, h)
    fineness_ratio = rocket_length/(rocket_ref_diam/2)
    coef1 = get_friction(Re, surface_roughness, rocket_length)
    coef2 = get_friction_compressed(coef1, mach)
    param = 1/(2*fineness_ratio)
    coef3 = get_scaled_friction_coef(coef2, param, wet_area, utils.get_circ_area(rocket_ref_diam))
    return coef3


def get_base_cd(mach):
    if mach < 1:
        return 0.12 + 0.13*mach**2
    else:
        return 0.25/mach
def get_stag_cd(mach):
    return .85*(1 + mach**2/4 + mach**4/40)

def get_von_karman_pressure_drag(length, base_diam, rocket_ref_diam):
    # Para a coifa de Von Karman, temos coef. de arrasto de pressão nulo no começo da região transônica
    # Portanto, basta usar a eq. para o regime subsônico
    dx = 0.00001
    theta_1 = np.arccos(1 - 2*length/length)
    y_1 = base_diam/2 * np.sqrt(theta_1 - np.sin(2*theta_1)/2)/np.sqrt(np.pi)

    theta_2 = np.arccos(1 - 2*(length-dx)/length)
    y_2 = base_diam/2 * np.sqrt(theta_2 - np.sin(2*theta_2)/2)/np.sqrt(np.pi)
    phi = np.arctan((y_1-y_2)/dx)
    unscaled_coef = 0.8*np.sin(phi)**2
    scaled_coef = unscaled_coef*utils.get_circ_area(base_diam)/utils.get_circ_area(rocket_ref_diam)
    return scaled_coef
  
def get_boattail_pressure_drag(length, diam_1, diam_2,rocket_ref_diam, mach = 0.3):
    aft_diam = np.minimum(diam_1, diam_2)
    fore_diam = np.maximum(diam_1,diam_2)
    gamma = length/np.abs(diam_1 - diam_2)
    param = 0
    if gamma <= 1:
        param = 1
    elif gamma > 1 and gamma < 3:
        param = (3-gamma)/2
    boattail_area = np.pi*(fore_diam/2)**2 - np.pi*(aft_diam/2)**2
    cd = get_base_cd(mach)*boattail_area*param/utils.get_circ_area(rocket_ref_diam)
    return cd

def get_end_base_drag(end_diam, rocket_ref_diam, mach = 0.3):
    return get_base_cd(mach)*np.pi*end_diam**2/4/utils.get_circ_area(rocket_ref_diam)

def get_fins_pressure_drag(span, thickness,sweep_angle, rocket_ref_diam ,mach = 0.3, no_fins = 1, type = 'square'):
    cd_le = 0
    cd_te = 0
    cd_le_perp = 0
    if type == 'airfoil' or type == 'rounded':
        cd_le_perp = (1-mach**2)**(-0.417) - 1
        cd_le = cd_le_perp*np.cos(np.deg2rad(sweep_angle))**2
        cd_te = 0
        if type == 'rounded':
            cd_te = get_base_cd(mach)/2
    else:
        cd_le_perp = get_stag_cd(mach)
        cd_le = cd_le_perp*np.cos(np.deg2rad(sweep_angle))**2
        cd_te = get_base_cd(mach)
    return (cd_le + cd_te)*no_fins*span*thickness/utils.get_circ_area(rocket_ref_diam)

def get_launch_lug_drag(length, diam_1, diam_2, rocket_ref_diam, mach = 0.3):
    ld = length/np.maximum(diam_1,diam_2)
    cd = np.maximum(1.3 - ld, 1)*get_stag_cd(mach)
    area = np.pi*(np.maximum(diam_1, diam_2)/2)**2 - (np.pi*(np.minimum(diam_1, diam_2)/2)**2 * np.maximum(1 - ld, 0))
    return cd * area/utils.get_circ_area(rocket_ref_diam)


####### Fim de definições de funções


## Função a ser chamada na simulação:

def get_rocket_cd(mach, h = 0):
    # Parâmetros do foguete
    # Coifa
    nose_length = 0.15
    nose_diam = 0.07
    nose_roughness = 60e-6

    # Fairing (corpo da carga-paga)
    fairing_length = 0.16
    fairing_diam = 0.07
    fairing_roughness = 60e-6

    # Transição Fairing->Corpo
    transition_length = 0.08
    transition_diam_1 = 0.07
    transition_diam_2 = 0.043
    transition_roughness = 60e-6

    # Corpo
    body_length = 0.40
    body_diam = 0.043
    body_roughness = 60e-6

    # Aletas
    fin_count = 3
    fin_chord_root = 0.11
    fin_chord_tip = 0.03
    fin_sweep_angle = 41.6
    fin_type = 'airfoil'
    fin_roughness = 60e-6
    fin_span = 0.045
    fin_thickness = 0.0044

    # Encaixes da guia de lançamento
    launchlug_count = 2
    launchlug_diam_1 = 0.007
    launchlug_diam_2 = 0.005
    launchlug_length = 0.025


    # Foguete inteiro
    rocket_ref_diam = np.max([nose_diam, fairing_diam, transition_diam_1, transition_diam_2, body_diam])
    rocket_length = nose_length + fairing_length + transition_length + body_length


    ## Arrasto da coifa
    nose_pressure_cd = get_von_karman_pressure_drag(nose_length, nose_diam, rocket_ref_diam)
    nose_friction_cd = get_von_karman_nosecone_friction_drag(nose_length, nose_diam, rocket_length, rocket_ref_diam, mach, h, nose_roughness)
   

    ## Arrasto do fairing
    fairing_friction_cd = get_body_friction_drag(fairing_length, fairing_diam, rocket_length, rocket_ref_diam, mach, h, fairing_roughness)


    # Arrasto da transição fairing -> corpo
    transition_pressure_cd = get_boattail_pressure_drag(transition_length, transition_diam_1, transition_diam_2, rocket_ref_diam, mach)
    transition_friction_cd = get_transition_friction_drag(transition_length, transition_diam_1, transition_diam_2, rocket_length, rocket_ref_diam, mach, h , transition_roughness)


    # Arrasto do corpo
    body_friction_cd = get_body_friction_drag(body_length, body_diam, rocket_length, rocket_ref_diam, mach, h, body_roughness)
    body_base_cd = get_end_base_drag(body_diam, rocket_ref_diam, mach)
 
    # Arrasto das aletas
    fins_friction_cd = get_fins_friction_drag(fin_chord_root, fin_chord_tip, fin_span, fin_thickness, rocket_length, rocket_ref_diam, mach, fin_count, h, fin_roughness)
    fins_pressure_cd = get_fins_pressure_drag(fin_span, fin_thickness, fin_sweep_angle, rocket_ref_diam, mach, fin_count, fin_type)

    # Arrasto dos encaixes
    launchlug_pressure_cd = get_launch_lug_drag(launchlug_length, launchlug_diam_1, launchlug_diam_2, rocket_ref_diam, mach)
    
    total_friction_cd = nose_friction_cd + fairing_friction_cd + transition_friction_cd + body_friction_cd + fins_friction_cd
    total_pressure_cd = nose_pressure_cd + transition_pressure_cd + fins_pressure_cd + launchlug_count * launchlug_pressure_cd
    total_base_cd = body_base_cd
    total_cd = total_friction_cd + total_pressure_cd + total_base_cd

    
    return total_cd




