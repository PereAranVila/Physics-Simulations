
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate
import numpy.random as random


###############################################################################
# Càlcul del valor de pi fent series numèriques i gràfica del error comés al
# càlcul d'aquest en funció de la "k" del sumatori.
 
k_llista = []  
error_llista = [] 
k = 0
valor_calculat_pi = 0
error = 1

while error > 0.0000000000001:
    valor_calculat_pi = (valor_calculat_pi + 
            +((1/16)**k)*((4/(8*k+1)-(2/(8*k+4))-(1/(8*k+5))-(1/(8*k+6)))))
    error = np.abs(np.pi-valor_calculat_pi)
    k_llista.append(k)
    error_llista.append(error)
    k = k +1

print('Valor de pi exacte: ', np.pi)
print('Valor de pi calculat: ', valor_calculat_pi)
print('Amb un error de: ', error)


plt.figure()
line_1 = plt.plot(k_llista, error_llista, label = 'error')
plt.legend(line_1, 'error')
plt.semilogy()
plt.grid()
plt.legend()
plt.title('Error de pi')
plt.xlabel('Numero Iteracions')
plt.ylabel('Error')
plt.savefig('error-pi-series.png')
plt.show()


# Guardem els valors generats en format .xlsx.
xx = np.array(error_llista)
yy = np.array(k_llista)

array_2d = np.transpose(np.vstack((xx,yy)))
df = pd.DataFrame(array_2d)
df.to_excel('calculs-series.xlsx', header = False, index = False)

###############################################################################
# Càlcul de pi mitjançant integrals (analítica == quad(), numèrica == simpson())

# Integrals Analítiques (quad())
def gaussiana(x):
    return np.exp(-x**2)

def racional(x):
    return (1/(1+x**2))

integral_quad_gauss = integrate.quad(gaussiana, -np.Inf, np.Inf)
print('Càlcul de pi, integral de Gauss: ', integral_quad_gauss)
integral_quad_racional = integrate.quad(racional, -np.Inf, np.Inf)
print('Calcul de pi, integral racional', integral_quad_racional)

# Integrals numèriques (simpson())
x_1 = np.linspace(-100000, 100000, 1000000)
y_1 = gaussiana(x_1)
y_2 = racional(x_1)
integral_simpson_gauss = integrate.simpson(y_1,x_1)
print('Integral gaussiana simpson',integral_simpson_gauss)
integral_simpson_racional = integrate.simpson(y_2, x_1)
print('Integral racional simpson',integral_simpson_racional)

###############################################################################
# Calcul de pi com a punts (quocients entre l'àrea del cercle de radi=1 i un 
# quadrat de L = 1)

num_iteracions_llista = []
error_llista_area = []
error = 1               # incialitzem amb un valor aleatori per el "while"
contador_dintre = 0 
contador_total = 0 

while error > 0.000001:
    x = random.rand()    # x pertany [0,1]
    y = random.rand()
    r = (x**2+y**2)**0.5
    if r <= 1:
        contador_dintre = contador_dintre +1
        contador_total = contador_total +1
    if r > 1:
        contador_total = contador_total +1 
        
    pi_calculat = 4*(contador_dintre/contador_total)
    error = np.abs(pi_calculat-np.pi)
    num_iteracions_llista.append(contador_total)
    error_llista_area.append(error)

print('Valor de pi calculat amb les àreas', pi_calculat)
plt.figure()
plt.subplot(2,1,1)
plt.plot(num_iteracions_llista, error_llista_area)
plt.yscale('log')
plt.grid()
plt.xlabel('Num iteracions')
plt.ylabel('Error absolut')
plt.title('Càlcul de pi per Àreas')
plt.subplot(2,1,2)
plt.plot(k_llista, error_llista)
plt.yscale('log')
plt.grid()
plt.xlabel('Num Iteracions')
plt.ylabel('Error absolut')
plt.title('Error utiitzant series')
plt.tight_layout()
plt.show()
plt.savefig('error-pi-areas-vs-series.png')