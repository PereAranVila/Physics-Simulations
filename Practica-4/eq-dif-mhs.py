
import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Osciladors Harmonics (MHS). Resoluccio de l'equaccio diferencial per un
# pendol simple.

def equadif_mhs(y,t):
    
    '''
    y[0] = velocitat
    y[1] = angle
    '''
    g = 9.8
    L = 1
    return -g * y[1] / L,  y[0]

t = np.linspace(0, 2*np.pi, 100)
ci = np.array([1, 0])
resultat = spi.odeint(equadif_mhs, ci, t)
vel_ang = resultat[:,0]
angle = resultat[:,1]

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, vel_ang)
plt.title('Pèndol Simple')
plt.xlabel('t (s)')
plt.ylabel('w (rad/s)')
plt.subplot(2,1,2)
plt.plot(t,angle)
plt.ylabel(r"$\theta (rad)$")
plt.xlabel('t (s)')
plt.savefig('mhs-pendol-simple.png')
plt.show()