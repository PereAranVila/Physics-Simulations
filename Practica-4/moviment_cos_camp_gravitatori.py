
import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Moviment d'un Cos en un Camp Gravitatori. En el nostre cas farem les
# simulacions per determinar les posicions (x,y) i velocitats (vx,vy) de la
# Terra i el cometa Halleys

def eq_dif_cos_camp_gravitatori(y,t):
    
    '''
    y[0] = velocitat en el eix x
    y[1] = posicio en el eix x
    y[2] = velocitat en el eix y
    y[3] = posicio en el eix y
    '''
    
    # paràmetres del sistema físic en SI
    G = 6.67e-11
    M = 1.9891e30
    
    R = (y[1]**2+y[3]**2)**0.5
    
    return -(G*M*y[1])/(R**3), y[0], -(G*M*y[3])/(R**3), y[2]

# Moviment pel cas de la Terra al voltant del Sol
ci_terra = [0,1.52098e11,2.93e4,0]
t = np.linspace(0,365.242198*24*3600,500)

resultat_terra = spi.odeint(eq_dif_cos_camp_gravitatori, ci_terra, t)

v_x_terra = resultat_terra[:,0]
posicio_x_terra = resultat_terra[:,1]
v_y_terra = resultat_terra[:,2]
posicio_y_terra = resultat_terra[:,3]

modul_velocitat = (v_x_terra**2 + v_y_terra**2)**0.5

plt.figure()
plt.suptitle('Terra-Sol')
plt.subplot(2,2,1)
plt.plot(t, posicio_x_terra)
plt.xlabel('t(s)')
plt.ylabel('x(m)')
plt.subplot(2,2,2)
plt.plot(t, posicio_y_terra)
plt.xlabel('t(s)')
plt.ylabel('y(m)')
plt.subplot(2,2,3)
plt.plot(t, modul_velocitat)
plt.xlabel('t(s)')
plt.ylabel('|v(m/s)|')
plt.subplot(2,2,4)
plt.plot(posicio_x_terra, posicio_y_terra)
plt.xlim(-0.5,0.5)
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.axis('equal')
plt.tight_layout()
plt.savefig('simulacio-Terra-Sol.png')
plt.show()

# Moviment del cometa Halley al voltant del Sol
ci_halley = [0,5.24e12, 1000, 0]
t = np.linspace(0,76*365.242198*24*3600, 300)

resultat_halley = spi.odeint(eq_dif_cos_camp_gravitatori, ci_halley, t)

v_x_halley = resultat_halley[:,0]
posicio_x_halley = resultat_halley[:,1]
v_y_halley = resultat_halley[:,2]
posicio_y_halley = resultat_halley[:,3]

velocitat_global = (v_x_halley**2 + v_y_halley**2)**0.5

plt.figure()
plt.suptitle('Halley-Sol')
plt.subplot(2,2,1)
plt.plot(t, posicio_x_halley)
plt.xlabel('t(s)')
plt.ylabel('x(m)')
plt.subplot(2,2,2)
plt.plot(t, posicio_y_halley)
plt.xlabel('t(s)')
plt.ylabel('y(m)')
plt.subplot(2,2,3)
plt.plot(t, velocitat_global)
plt.xlabel('t(s)')
plt.ylabel('|v(m/s)|')
plt.subplot(2,2,4)
plt.plot(posicio_x_halley, posicio_y_halley)
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.axis('equal')
plt.tight_layout()
plt.savefig('simulacio-Halley-Sol.png')
plt.show()
