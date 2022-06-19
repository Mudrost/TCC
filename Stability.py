import numpy as np

# Parâmetros da coifa
m_nosecone = .0152                              # Massa da coifa (kg)
l_nosecone = 0.15                               # Comprimento da coifa (m)
d_nosecone = 0.07                                # Diâmetro da base da coifa (m)
xcg_nosecone = (2/3)*l_nosecone                   # Distância do CG da coifa em relação à ponta da coifa (m) 
cn_n = 2                                        # Coef. para a coifa
x_n = l_nosecone*(1 - 0.500)                      # Distância relativa, coifa Von Karman

# Parâmetros da carenagem
m_fairing = .0642                           
l_fairing = 0.16                              
d_fairing = 0.07
xcg_fairing = 0.5*l_fairing                      
xcg_fairing_ref = l_nosecone +xcg_fairing      
r_fairing_end = d_fairing/2
# Parâmetros da transição
m_transition = 0.0262
l_transition = 0.08 
d_fore_transition = 0.07  
d_aft_transition = 0.043
x_p = l_nosecone + l_fairing                                                                                                                          
xcg_transition = 0.03676
xcg_transition_ref = l_nosecone + l_fairing + xcg_transition  
cn_t = 2*((d_aft_transition/d_nosecone)**2 - (d_fore_transition/d_nosecone)**2 )
x_t = x_p + l_transition/3 * (1 + (1 - d_fore_transition/d_aft_transition)/(1 - (d_fore_transition/d_aft_transition)**2))

# Parâmetros do corpo principal
m_body = .191                            
l_body = 0.4                              
d_body = 0.043
xcg_body = 0.5*l_body                       
xcg_body_ref = l_nosecone + l_fairing + l_transition + xcg_body        
r_body_end = d_body/2

# Parâmetros do motor
m_motor = 0.595
l_motor = 0.29
d_motor = 0.38
xcg_motor = 0.5*l_motor
xcg_motor_ref = l_nosecone + l_fairing + l_transition + l_body - l_motor + xcg_motor

# Parâmetros da aleta
n_fins = 3                                        
m_fin = 0.00557                                        
s_fin = 0.045                                           
lf_fin = 0.045                                           
cr_fin = 0.11                                           
ct_fin = 0.03                                            
xcg_fin = 0.5*cr_fin
xcg_fin_ref = l_nosecone + l_fairing + l_transition + l_body - cr_fin + xcg_fin
x_r = (cr_fin - ct_fin)/2                                           
x_b = l_nosecone + l_fairing + l_transition + l_body - cr_fin                                          
cn_f = (1 + (r_body_end/(r_body_end + s_fin)))*( (4*n_fins*(s_fin/d_nosecone)**2)/(1 + np.sqrt(1 + (2*lf_fin/(cr_fin+ct_fin))**2)))
x_f = x_b + (x_r/3) *((cr_fin + 2*ct_fin)/(cr_fin + ct_fin)) + (1/6) * ((cr_fin + ct_fin) - (cr_fin*ct_fin)/(cr_fin + ct_fin))

## COMPONENTES INTERNOS
# Parâmetros dos eletronicos
m_bay = 0.3
l_bay = 0.025
xcg_bay = l_bay/2
xcg_bay_ref = xcg_bay + 0.033

# Paraquedas
m_parachute = 0.038
l_parachute = 0.075
xcg_parachute = l_parachute/2
xcg_parachute_ref = xcg_parachute + 0.07

# CanSat (externo)
m_cansat_external = 0.2
l_cansat_external = 0.025
xcg_cansat_external = l_cansat_external/2
xcg_cansat_external_ref = l_nosecone + xcg_cansat_external + 0.013

# CanSat (interno)
m_cansat_internal = 0.35
l_cansat_internal = 0.115
xcg_cansat_internal = l_cansat_internal/2
xcg_cansat_internal_ref = l_nosecone + xcg_cansat_internal+ 0.045

# Carga 
m_charge = 0.2
l_charge = 0.055
xcg_charge = l_charge/2
xcg_charge_ref = l_nosecone + l_fairing + xcg_charge + 0.01

# Centro de gravidade
xcg = (m_nosecone * xcg_nosecone + m_body * xcg_body_ref + m_transition * xcg_transition_ref + m_motor * xcg_motor_ref +\
     n_fins*m_fin * xcg_fin_ref + m_bay*xcg_bay_ref + m_parachute*xcg_parachute_ref + m_cansat_internal*xcg_cansat_internal_ref + \
        m_cansat_external*xcg_cansat_external_ref + m_charge*xcg_charge_ref)/ \
            (m_nosecone + m_body + m_transition + m_motor + n_fins*m_fin + m_bay + m_parachute + m_cansat_internal + m_cansat_external + m_charge)

xcp = (cn_t * x_t + cn_f*x_f + cn_n*x_n)/(cn_t + cn_f + cn_n)

static_margin = (xcp-xcg)/d_fairing


print("Xcg em relação à ponta da coifa: ", xcg, "m")
print("Xcp em relação à ponta da coifa: ", xcp, "m")
print("Margem de estabilidade estática: ", static_margin, " cal")