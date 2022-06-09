# Apêndice B - Stability.py

import numpy as np


# Parâmetros da coifa
m_nosecone = .200                               # Massa da coifa (kg)
l_nosecone = 0.20                               # Comprimento da coifa (m)
d_nosecone = 0.1                                # Diâmetro da base da coifa (m)
xcg_nosecone = 0.6*l_nosecone                   # Distância do CG da coifa em relação à ponta da coifa (m) 
cn_n = 2                                        # Coef. para a coifa
xn = l_nosecone(1 - 0.500)                      # Distância relativa, coifa Von Karman


# Parâmetros do tubo
m_body = .200                               # Massa do corpo (kg)
l_body = 0.20                               # Comprimento do corpo (m)
d_body = 0.2
xcg_body = 0.6*l_body                       # Distância do CG do corpo em relação à base da coifa (m) 
xcg_body_ref = l_nosecone + xcg_body        # Distância do CG do corpo em relação à linha de referencia na ponta da coifa (m)
r_body_end = d_body/2


# Parâmetros de transição (se houver)
m_transition = 0.1
l_transition = 0.1 
d_front_transition = 0.1   
d_rear_transition = 0.8
x_p = 0.06                                                                                                                              # Distância entre a ponta da coifa e o começo da transição
xcg_transition = 0.5*l_transition
xcg_transition_ref = l_nosecone + 0.5  
cn_t = 2*((d_rear_transition/d_nosecone)**2 - (d_front_transition/d_nosecone)**2 )
x_t = x_p + l_transition/3 * (1 + (1 - d_front_transition/d_rear_transition)/(1 - (d_front_transition/d_rear_transition)**2))

# Parâmetros do motor
m_motor = 1.5
l_motor = 0.3
d_motor = 0.1
xcg_motor = 0.55*l_motor
xcg_motor_ref = l_nosecone + l_body - l_motor + xcg_motor


# Parâmetros da aleta
n_fins = 4                                          # Número de aletas
m_fin = 0.2                                         # Massa da aleta
s_fin = 0.5                                         # Envergadura da aleta
lf_fin = 0.7                                        # Comprimento da linha de corda média
cr_fin = 0.6                                        # Corda na raiz
ct_fin = 0.4                                        # Corda na ponta
xcg_fin = 0.6*cr_fin
xcg_fin_ref = l_nosecone + l_body - cr_fin + xcg_fin
x_r = 0.4                                           # Distância paralela ao corpo, entre o começo do bordo de ataque e o final do bordo de ataque 
x_b = 0.8                                           # Distância entre a ponta da coifa e o começo do bordo de ataque da aleta
cn_f = (1 + (r_body_end/(r_body_end + s_fin)))*( (4*n_fins*(s_fin/d_nosecone)**2)/(1 + np.sqrt(1 + (2*lf_fin/(cr_fin+ct_fin))**2)))
x_f = x_b + x_r/3 *((cr_fin + 2*ct_fin)/(cr_fin + ct_fin)) + 1/6 * ((cr_fin + ct_fin) - (cr_fin*ct_fin)/(cr_fin + ct_fin))

# Centro de gravidade
xcg = m_nosecone * xcg_nosecone + m_body * xcg_body_ref + m_transition * xcg_transition_ref + m_motor * xcg_motor_ref + m_fin * xcg_