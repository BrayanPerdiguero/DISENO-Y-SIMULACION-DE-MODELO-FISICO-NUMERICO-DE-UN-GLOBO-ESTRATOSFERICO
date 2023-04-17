from math import *
from numpy import *
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import axes3d
import random 
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
m_film = 30  # [Kg]
m_parachute = 270 # [Kg]
m_c = 5 # [Kg]
m_PL = 265 # [Kg]

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
z_e0 = 0.1
m0_he = 330 # [Kg]
t_Interval = linspace(0,t_f)
y_01 = [0,0,z_e0,0,0,0,m0_he]

                ######## Modelo de Viento Atmosférico ########

z_0 = 0.1
C_D = 0.0075
k = 0.4
M_10 = 4.7

U_fricc = sqrt(C_D * M_10**2)

def viento(z):
    V_wx = random.uniform(-((U_fricc / k) * log(z/z_0)), (U_fricc / k) * log(z/z_0))
    V_wy = random.uniform(-((U_fricc / k) * log(z/z_0)), (U_fricc / k) * log(z/z_0))
    V_wz = random.uniform(-5, 5)
    return [V_wx, V_wy, V_wz]


def f3D(t,y):
    

###### Variables lado izquierdo, deribadas respecto al tiempo #####
    x_e = y[0]
    y_e = y[1]
    z_e = y[2]

    Vx = y[3]
    Vy = y[4]
    Vz = y[5]

    m_he = y[6]
########### Valores y dirección del viento como inputs en funcion de la altura #####
   
    V_wind = viento(z_e) 

    alpha = random.uniform(-pi/2, pi/2)                     #[rad]
    phi = random.uniform(-2*pi, 2*pi)                      #[rad]

######## Caracterización de la valvula de escape ########
    C1 = 0.6
    D_val = 0.2                    #[m]
    A1 = pi*(D_val/2)**2


########## Expresiones del modelo fisico ##################

    F1 = Vx*sin(alpha)*cos(phi) + V_wind[0]                        #[m/s]
    F2 = Vy*sin(alpha)*sin(phi) + V_wind[1]                       #[m/s]
    F3 = Vz*cos(alpha) + V_wind[2]                                #[m/s]
    F4 = D_b(z_e,Vx,m_he)*cos(phi)*sin(alpha)/m_T(z_e,m_he)
    F5 = D_b(z_e,Vy,m_he)*sin(phi)*sin(alpha)/m_T(z_e,m_he)
    F6 = (L_b(z_e,m_he)-W_T(z_e,m_he)-D_b(z_e,Vz,m_he)*cos(alpha))/m_T(z_e,m_he)
    F7 = - C1*A1*(2*0.01*rho_he(z_e))**0.5
    
    return [F1,F2,F3,F4,F5,F6,F7]



sol1 = solve_ivp( fun= f3D, t_span=(t_0,t_f), y0=y_01, t_eval=t_Interval)

############################## Representación de resultados #####################################

plt.subplot(2,3,1)
plt.plot(sol1.t,sol1.y[2],'b',)
plt.xlabel('time [s]')
plt.ylabel('altitud [m]')
plt.grid()

plt.subplot(2,3,2)
plt.plot(sol1.y[2],sol1.y[5],'b',)
plt.xlabel('altitud [m]')
plt.ylabel('velocity [m/s]')
plt.grid()

plt.subplot(2,3,3)
plt.plot(sol1.t,sol1.y[5],'b',)
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.grid()

plt.subplot(2,3,4)
plt.plot(sol1.t,sol1.y[6],'b',)
plt.xlabel('time [s]')
plt.ylabel('masa_sistema [Kg]')
plt.grid()

plt.subplot(2,3,5)
plt.plot(sol1.t,sol1.y[0],'b',)
plt.xlabel('time [s]')
plt.ylabel(' X [m]')
plt.grid()

plt.subplot(2,3,6)
plt.plot(sol1.t,sol1.y[1],'b',)
plt.xlabel('time [s]')
plt.ylabel(' Y [m]')
plt.grid()

plt.show()



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


# Mostramos el gráfico
plt.show()