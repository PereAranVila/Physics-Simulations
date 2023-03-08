
import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Oscilador Esmorteit Forçat amb A=0 i b=0

def equacio_dif(y,t):
    
    '''
    y[0] = velocitat
    y[1] = angle
    '''
    
    L = 1
    g = 9.8
    
    
    return -g*np.sin(y[1]) /L, y[0]

t = np.linspace(0, 4*np.pi, 1000)
ci = np.array([1,0])

resultat_0 = spi.odeint(equacio_dif, ci, t)
vel_ang_0 = resultat_0[:,0]
angle_0 = resultat_0[:,1]

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, vel_ang_0)
plt.title('b = 0 , A = 0')
plt.xlabel('t (s)')
plt.ylabel('w (rad/s)')
plt.subplot(2,1,2)
plt.plot(t,angle_0)
plt.xlabel('t(s)')
plt.ylabel(r'$\theta (rad)$')
plt.tight_layout()
plt.savefig('osilador-esmorteit-forçat(A=0-i-b=0).png')
plt.show()
    
###############################################################################
# Oscilador Esmorteit Forçat amb b!=0 i A=0

def equacio_dif_a_0(y,t):
    
    '''
    y[0] = velocitat
    y[1] = angle
    '''
    L = 1
    m = 1
    b = 1
    g = 9.8
    
    return (-g*np.sin(y[1])/L) + (-b*y[0]/m*(L^2)), y[0]

resultat_1 = spi.odeint(equacio_dif_a_0, ci, t)
vel_ang_1 = resultat_1[:,0]
angle_1 = resultat_1[:,1]

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, vel_ang_1)
plt.title('Esmorteït')
plt.xlabel('t (s)')
plt.ylabel('w (rad/s)')
plt.subplot(2,1,2)
plt.plot(t,angle_1)
plt.xlabel('t(s)')
plt.ylabel(r'$\theta$ (rad)')
plt.ylim(-0.05,0.2)
plt.tight_layout()
plt.savefig('oscilador-esmorteit-forçat(b!=0-i-A=0)')
plt.show()

###############################################################################
# Oscilador Esmorteit Forçat amb b!=0 i A!=0

def equacio_dif_tercer(y,t):
    
    '''
    y[0] = velocitat
    y[1] = angle
    '''

    L = 1
    m = 1
    b = 0.5
    A = 1.35
    g = 1
    omega = 0.666
    
    return (-g*np.sin(y[1])/L) + (-b*y[0] + A*np.cos(omega*t))/m*(L^2), y[0]

t = np.linspace(0, 10*np.pi, 1000)

resultat_2 = spi.odeint(equacio_dif_tercer, ci, t)
vel_ang_2 = resultat_2[:,0]
angle_2 = resultat_2[:,1]

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, vel_ang_2)
plt.title('Forçat')
plt.xlabel('t (s)')
plt.ylabel('w (rad/s)')
plt.subplot(2,1,2)
plt.plot(t,angle_2)
plt.xlabel('t(s)')
plt.ylabel(r'$\theta$ (rad)')
plt.tight_layout()
plt.savefig('oscilador-esmorteit-forçat-(b!=0-i-A=0)')
plt.show()
