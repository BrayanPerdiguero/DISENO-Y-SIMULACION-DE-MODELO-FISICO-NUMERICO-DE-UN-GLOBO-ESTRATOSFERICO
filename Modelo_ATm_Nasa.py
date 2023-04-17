from math import *
from numpy import *
from matplotlib import pyplot as plt

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



h = range(0,50000)


plt.subplot(3,1,1)
plt.plot([rho_air(i) for i in h],h)
plt.xlabel('Densidad [Kg/m^3]')
plt.ylabel('h [m]')
plt.grid()


plt.subplot(3,1,2)
plt.plot([T_air(i) for i in h],h)
plt.xlabel('Temperatura [K]')
plt.ylabel('h [m]')
plt.grid()


plt.subplot(3,1,3)
plt.plot([P_air(i) for i in h],h)
plt.xlabel('Presion [Pa]')
plt.ylabel('h [m]')
plt.grid()


plt.show()



