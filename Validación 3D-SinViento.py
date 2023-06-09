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
    if (z<0):
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
    if (z<0):
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
m_film = 10  # [Kg]
m_parachute = 0 # [Kg]
m_c = 6 # [Kg]
m_PL = 170 # [Kg]

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
m0_he = 109 # [Kg]
t_Interval = linspace(0,t_f)
y_01 = [0,0,0,0,0,0,m0_he]
y_02 = [0,0,0,0,0,0,325]
y_03 = [0,0,0,0,0,0,320]

    ######## Componentes vector Viento en m/s ########
V_wx = 0
V_wy = 0
V_wz = 0

alpha = 0
phi = 0

def f3D(t,y):
    
    x_e = y[0]
    y_e = y[1]
    z_e = y[2]

    Vx = y[3]
    Vy = y[4]
    Vz = y[5]

    m_he = y[6]

    C1 = 0.6
    D_val = 0.2  #[m]
    A1 = pi*(D_val/2)**2
    
    F1 = Vx*sin(alpha)*cos(phi) + V_wx
    F2 = Vy*sin(alpha)*sin(phi) + V_wz
    F3 = Vz*cos(alpha) + V_wz
    F4 = D_b(z_e,Vx,m_he)*cos(phi)*sin(alpha)/m_T(z_e,m_he)
    F5 = D_b(z_e,Vy,m_he)*sin(phi)*sin(alpha)/m_T(z_e,m_he)
    F6 = (L_b(z_e,m_he)-W_T(z_e,m_he)-D_b(z_e,Vz,m_he)*cos(alpha))/m_T(z_e,m_he)
    F7 = - C1*A1*(2*0.00129*rho_he(z_e))**0.5
    
    return [F1,F2,F3,F4,F5,F6,F7]



sol1 = solve_ivp( fun= f3D, t_span=(t_0,t_f), y0=y_01, t_eval=t_Interval)

#sol2 = solve_ivp( fun= f3D, t_span=(t_0,t_f), y0=y_02, method= 'RK45', t_eval=t_Interval)

#sol3 = solve_ivp( fun= f3D, t_span=(t_0,t_f), y0=y_03, method= 'RK45', t_eval=t_Interval)

############################## Representación de resultados #####################################

plt.subplot(2,3,1)
plt.plot(sol1.t,sol1.y[2],'b')
#plt.plot(sol1.t,sol1.y[2],'b',sol2.t,sol2.y[2],'r',sol3.t,sol3.y[2],'g')
plt.xlabel('time [s]')
plt.ylabel('altitud [m]')
#plt.legend(['caso 1', 'caso 2', 'caso 3'],loc='upper left')
plt.grid()

plt.subplot(2,3,2)
plt.plot(sol1.y[2],sol1.y[5],'b')
#plt.plot(sol1.y[2],sol1.y[5],'b',sol2.y[2],sol2.y[5],'r',sol3.y[2],sol3.y[5],'g')
plt.xlabel('altitud [m]')
plt.ylabel('velocity [m/s]')
#plt.legend(['caso 1', 'caso 2', 'caso 3'],loc='upper left')
plt.grid()

plt.subplot(2,3,3)
plt.plot(sol1.t,sol1.y[5],'b')
#plt.plot(sol1.t,sol1.y[5],'b',sol2.t,sol2.y[5],'r',sol3.t,sol3.y[5],'g')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
#plt.legend(['caso 1', 'caso 2', 'caso 3'],loc='upper left')
plt.grid()

plt.subplot(2,3,4)
plt.plot(sol1.t,sol1.y[6],'b')
#plt.plot(sol1.t,sol1.y[6],'b',sol2.t,sol2.y[6],'r',sol3.t,sol3.y[6],'g')
plt.xlabel('time [s]')
plt.ylabel('Helio [Kg]')
#plt.legend(['caso 1', 'caso 2', 'caso 3'],loc='upper right')
plt.grid()

fig = plt.figure()

# Agregamos un plano 3D
ax1 = fig.add_subplot(111,projection='3d')
plt.title('Trayectoria del Globo')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [m]')

# Datos en array bi-dimensional
x = array([sol1.y[0]])
y = array([sol1.y[1]])
z = array([sol1.y[2]])

# plot_wireframe nos permite agregar los datos x, y, z. Por ello 3D
# Es necesario que los datos esten contenidos en un array bi-dimensional
ax1.plot_wireframe(x, y, z)

plt.show()
