
import scipy.integrate as spi
import scipy.special as sps
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Resoluccio de l'equaccio diferencial de Schrodinguer per a l'Oscilador
# Harmònic. 

def equacio_dif_osc_quantic(y,x):
    
    ''' y[0] = velocitat
        y[1] = posicio '''
    
    global n
    B = 1
    
    return -((2*n+1)-B*x**2)*y[1], y[0]

valors_n = [0,2,4]
condicions_incials = [[0, 0.751126],[0, -0.531125],[0, 0.459969]]

plt.figure()
plt.subplot(2,1,1)
for j in range(0,3):
    n = valors_n[j]
    ci = condicions_incials[j]
    
    # per x entre [0,-5]
    x_neg = np.linspace(0,-5,1000)
    x_neg_fliped = np.flip(x_neg) 
    resultat_neg = spi.odeint(equacio_dif_osc_quantic, ci, x_neg)
    posicio_neg = resultat_neg[:,1]
    posicio_neg_fliped = np.flip(posicio_neg)

    # per x entre[0,5] i n = 0
    x_pos = np.linspace(0,5,1000)
    resultat_pos = spi.odeint(equacio_dif_osc_quantic, ci, x_pos)
    posicio_pos = resultat_pos[:,1]

    # concateno per tal de graficar les arrays x i la posicio
    x = np.concatenate((x_neg_fliped, x_pos ))
    posicio = np.concatenate((posicio_neg_fliped, posicio_pos))

    plt.plot(x, posicio,label='n='+str(n))
    
plt.xlim(-5,5)
plt.ylim(-0.6,0.8)
plt.grid()
plt.xlabel('x')
plt.ylabel('$\psi$')
plt.title("Numèrica")
plt.legend()

###############################################################################
# Mostrem la soluccio analítica (la qual interve els polinomis d'Hermite) de
# l'equaccio de Schrodinguer per a un oscilador harmonic


def analitic_schrodinguer(x,n):
    '''Aquesta funció te per entrada el valor de les x i n és un paràmtre
    que intervé en la solució de l'equació'''
    
    H_n = sps.eval_hermite(n,x)  # poliomi d'Hermite
    
    return ((((2**n)*np.math.factorial(n)*(np.pi)**0.5)**-0.5) * 
            np.exp((-(x**2) / 2))*H_n)

plt.subplot(2,1,2)
for j in range(0,3):
    n = valors_n[j]  # agafo el valor de la n
    x = np.linspace(-5,5,2000)
    plt.plot(x, analitic_schrodinguer(x, n),label='n='+str(n))

plt.xlim(-5,5)
plt.ylim(-0.6,0.8)
plt.grid()
plt.xlabel('x')
plt.ylabel('$\psi$')
plt.title("Analítica")
plt.legend()
plt.tight_layout()
plt.savefig('eq-Schrodinguer-oscilador-harmonic.png')
plt.show()

    