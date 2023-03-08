
import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt 


def eq_dif_acoblats(x,t):
    
    '''
    x[0] = velocitat de x_1
    x[1] = posicio de x_1
    x[2] = velocitat de x_2
    x[3] = posicio de x_2
    '''
    # parametres del sistema
    m = 1
    k_1 = 10
    k_2 = 0.5
    
    return(-((k_1+k_2)*x[1]/m) +(k_2*x[3])/m ,x[0], -((k_1+k_2)*x[3])/m + (k_1*x[1])/m, x[2] )

ci = [0, 1, 0, 0]  # condicions inicials

t = np.linspace(0, 40, 10000)
resultat = spi.odeint(eq_dif_acoblats, ci, t)
v_1 = resultat[:,0]
x_1 = resultat[:,1]
v_2 = resultat[:,2]
x_2 = resultat[:,3]


# plot de les posicions x1 i x2 en funcio del temps i en una altre grafica 
# les velocitats v1 i v2 en funcio del temps
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,x_1,label='x1')
plt.plot(t,x_2,label='x2')
plt.grid()
plt.ylim(-5,5)
plt.xlim(0,40)
plt.xlabel('Temps (s)')
plt.ylabel('Posició objecte ')
plt.legend()
plt.subplot(2,1,2)
plt.plot(t,v_1,label='v1')
plt.plot(t, v_2,label='v2')
plt.grid()
plt.xlim(0,40)
plt.ylim(-15,15)
plt.xlabel('Temps (s)')
plt.ylabel('Velocitat objecte ')
plt.savefig('osciladors-acoblats.png')
plt.legend()
plt.show()