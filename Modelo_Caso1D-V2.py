from math import *
from numpy import *
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp


########################Modelo atmosférico ##########################################################################

def P_air(z):
    if (0<= z <= 11000):
        P = 101290 * (T_air(z)/288.08)**(5.256)
        return P
    if (11000< z <=25000):
        P = 22650 * exp(1.73 - 0.000157*z)
        return P
    if (25000< z ):
        P = 2488 * (T_air(z)/ 216.65)**(-11.388)
        return P
    else:
        return print('altitude error')

def T_air(z):
    if (0<= z <= 11000):
        T = 15.05 + 273.1 - 0.00649*z
        return T 
    if (11000< z <=25000):
        T = -56.46 + 273.1
        return T 
    if (25000< z ):
        T = -131.21 + 273.1 + 0.00299*z 
        return T 
    else:
        return print('Temperature error')    


def rho_air(z):
    rho = P_air(z) / (286.9*T_air(z))
    return rho
##################### Constantes y Parámetros ##############################################################
#rho_he = 0.167

R_gas = 8.3144   # [m^3·Pa / K mol]

g0 = 9.81
R_e = 6371 * 1000
Cv_air = 0.5


M_mol_he = 4.0026 # [g / mol]
m_film = 30  # [Kg]
m_parachute = 270 # [Kg]
m_c = 5 # [Kg]
m_PL = 300 # [Kg]

def Vol_b(z,m_he):
    Vol = ((m_he*1000) *R_gas*T_air(z)) / (M_mol_he*P_air(z))
    return Vol

def Dm_b(z,m_he):
    Dm = 1.383 * (Vol_b(z,m_he))**(1/3)
    return Dm

def S_b(z,m_he):
    S = pi * (Dm_b(z,m_he)/2)**2
    return S

def rho_he(z):
    rho_he = (P_air(z) * M_mol_he) / (R_gas * T_air(z))
    return rho_he

################## Fuerzas Aerodinamicas ##################################################################


def L_b(z,m_he):
    L = rho_air(z)*Vol_b(z,m_he)*g_z(z)
    return L


def D_b(z,Vz,m_he):
    D = 0.5*rho_air(z)*(Vz*Vz) * 0.5*S_b(z,m_he)
    return D
########################## Gravedad  y Pesos #######################################################################

def g_z (z):
    g = g0* (R_e/(R_e + z))**2
    return g


def m_T(z,m_he):
    m_virt = Cv_air*rho_air(z)*Vol_b(z,m_he)
    m_Total = m_film + m_he + m_c + m_PL + m_virt + m_parachute
    return m_Total

def W_T(z,m_he):
    W = m_T(z,m_he) * g_z(z)
    return W


########################## Solver ########################################################################
t_0 = 0
t_f = 12800
m0_he = 330 # [Kg]
t_Interval = linspace(0,t_f)
y_01 = [0.0,0.0,m0_he]
y_02 = [0.0,0.0,325]
y_03 = [0.0,0.0,320]
def f1D(t,y):
    
    z = y[0]
    Vz = y[1]
    m_he = y[2]

    C1 = 0.6
    D_val = 0.2  #[m]
    A1 = pi*(D_val/2)**2

    if y[1] > 0:      #Vertical Motion  
        z = y[0]
        Vz = y[1]
        m_he = y[2]
        return [Vz, (L_b(z,m_he) - W_T(z,m_he) - D_b(z,Vz,m_he)) / m_T(z,m_he), - C1*A1*(2*0.01*rho_he(z))**0.5]
    if y[1] <= 0:       # Stable Floatin
        z = y[0] 
        Vz  = 0
        m_he = y[2]
        return [Vz, (L_b(z,m_he) - W_T(z,m_he) - D_b(z,Vz,m_he)) / m_T(z,m_he), 0]



sol1 = solve_ivp( fun= f1D, t_span=(t_0,t_f), y0=y_01, method= 'RK45', t_eval=t_Interval)

sol2 = solve_ivp( fun= f1D, t_span=(t_0,t_f), y0=y_02, method= 'RK45', t_eval=t_Interval)

sol3 = solve_ivp( fun= f1D, t_span=(t_0,t_f), y0=y_03, method= 'RK45', t_eval=t_Interval)

############################## Representación de resultados #####################################

plt.subplot(2,3,1)
plt.plot(sol1.t,sol1.y[0],'b',sol2.t,sol2.y[0],'r',sol3.t,sol3.y[0],'g')
plt.xlabel('time [s]')
plt.ylabel('altitud [m]')
plt.legend(['He = 330 [Kg]', 'He = 325 [Kg]', 'He = 320 [Kg]'],loc='upper left')
plt.grid()

plt.subplot(2,3,2)
plt.plot(sol1.y[0],sol1.y[1],'b',)
plt.xlabel('altitud [m]')
plt.ylabel('velocity [m/s]')
plt.legend(['He = 330 [Kg]', 'He = 325 [Kg]', 'He = 320 [Kg]'],loc='upper right')
plt.grid()

plt.subplot(2,3,3)
plt.plot(sol1.t,sol1.y[1],'b',)
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend(['He = 330 [Kg]', 'He = 325 [Kg]', 'He = 320 [Kg]'],loc='upper right')
plt.grid()

plt.subplot(2,3,4)
plt.plot(sol1.t,sol1.y[2],'b',sol2.t,sol2.y[2],'r',sol3.t,sol3.y[2],'g')
plt.xlabel('time [s]')
plt.ylabel('Helio [Kg]')
plt.legend(['He = 330 [Kg]', 'He = 325 [Kg]', 'He = 320 [Kg]'],loc='upper right')
plt.grid()

plt.show()

