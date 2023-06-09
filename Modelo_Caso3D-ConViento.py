from math import *
from numpy import *
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import axes3d
import random 
import pandas as pd
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
z_e0 = 1
m0_he = 330 # [Kg]
t_Interval = linspace(0, t_f, num = 100)
y_01 = [0,0,z_e0,0,0,0,m0_he]

                ######## Modelo de Viento Atmosférico ########

z_0 = 0.1
C_D = 0.0075
k = 0.4
M_10 = 4.7

U_fricc = sqrt(C_D * M_10**2)


#def viento(z):
    
#    V_wx = random.uniform(-((U_fricc / k) * log(z/z_0)), (U_fricc / k) * log(z/z_0))
#    V_wy = random.uniform(-((U_fricc / k) * log(z/z_0)), (U_fricc / k) * log(z/z_0))
#    V_wz = random.uniform(-5, 5)
#    return [V_wx, V_wy, V_wz]
    

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
    
    alpha = random.uniform(-pi/2, pi/2)                     #[rad]
    phi = random.uniform(-2*pi, 2*pi)                       #[rad]
    if (random.randint(0,1) == True): 
        V_wx = random.uniform(-((U_fricc / k) * log(z_e/z_0)), (U_fricc / k) * log(z_e/z_0))
        V_wy = random.uniform(-((U_fricc / k) * log(z_e/z_0)), (U_fricc / k) * log(z_e/z_0))
        V_wz = random.uniform(-5, 5)  
        V_wind = [V_wx, V_wy, V_wz]
                             
    else:
        V_wx = 0
        V_wy = 0
        V_wz = 0
        V_wind = [V_wx, V_wy, V_wz]

  
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
plt.plot(sol1.t /3600,sol1.y[2] /1000,'b',)
plt.xlabel('time [Horas]')
plt.ylabel('altitud [Km]')
plt.grid()

plt.subplot(2,3,2)
plt.plot(sol1.y[2] /1000,sol1.y[5],'b',)
plt.xlabel('altitud [Km]')
plt.ylabel('velocity [m/s]')
plt.grid()

plt.subplot(2,3,3)
plt.plot(sol1.t /3600,sol1.y[5],'b',)
plt.xlabel('time [Horas]')
plt.ylabel('velocity [m/s]')
plt.grid()

plt.subplot(2,3,4)
plt.plot(sol1.t /3600,sol1.y[6],'b',)
plt.xlabel('time [Horas]')
plt.ylabel('masa_sistema [Kg]')
plt.grid()

plt.subplot(2,3,5)
plt.plot(sol1.t /3600,sol1.y[0],'b',)
plt.xlabel('time [Horas]')
plt.ylabel(' X [m]')
plt.grid()

plt.subplot(2,3,6)
plt.plot(sol1.t /3600,sol1.y[1],'b',)
plt.xlabel('time [Horas]')
plt.ylabel(' Y [m]')
plt.grid()

plt.show()

plt.subplot(1,3,1)
plt.plot(sol1.y[0],sol1.y[2] /1000,'b',)
plt.xlabel('X [m]')
plt.ylabel(' Z [Km]')
plt.grid()

plt.subplot(1,3,2)
plt.plot(sol1.y[1],sol1.y[2] /1000,'b',)
plt.xlabel('Y [m]')
plt.ylabel(' Z [Km]')
plt.grid()

plt.subplot(1,3,3)
plt.plot(sol1.y[0],sol1.y[1],'b',)
plt.xlabel('X [m]')
plt.ylabel(' Y [m]')
plt.grid()

plt.show()

fig = plt.figure()

# Agregamos un plano 3D
ax1 = fig.add_subplot(111,projection='3d')
plt.title('Trayectoria del Globo')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_zlabel('z [Km]')

# Datos en array bi-dimensional
x = array([sol1.y[0]])
y = array([sol1.y[1]])
z = array([sol1.y[2]/1000])

# plot_wireframe nos permite agregar los datos x, y, z. Por ello 3D
# Es necesario que los datos esten contenidos en un array bi-dimensional
ax1.plot_wireframe(x, y, z)


# Mostramos el gráfico
plt.show()

#############################   Pandas   #####################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

fig, ax_nstd = plt.subplots(figsize=(6, 6))

ax_nstd.axvline(c='grey', lw=1)
ax_nstd.axhline(c='grey', lw=1)

############### X e Y ####################
ax_nstd.scatter(sol1.y[0], sol1.y[1])

confidence_ellipse(sol1.y[0], sol1.y[1], ax_nstd, n_std=1, label=r'$1\sigma$', edgecolor='firebrick')
confidence_ellipse(sol1.y[0], sol1.y[1], ax_nstd, n_std=2, label=r'$2\sigma$', edgecolor='fuchsia', linestyle='--')
confidence_ellipse(sol1.y[0], sol1.y[1], ax_nstd, n_std=3, label=r'$3\sigma$', edgecolor='blue', linestyle=':')

ax_nstd.set_title('Elipse de Error')
ax_nstd.legend()
plt.xlabel('X [m]')
plt.ylabel('Y [m]')

plt.show()