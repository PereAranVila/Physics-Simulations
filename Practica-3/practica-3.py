
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.stats as stats
import pandas as pd

###############################################################################
# APARTAT A) Representacio de l'equaccio de Debye de les permitivitats
# dielectriques de certs materials com l'aigua en funcio de la frequencia per
# i per a certes temperatures determinades 

def epsilon (f,T):
    '''Calcula la formula 3.1 del guió de pràctiques, f és la freqüencia
    i T la temperatura'''
    
    # costants necessaries
    a_1 = 87.9
    b_1 = 0.404
    c_1 = 9.59e-4
    d_1 = 1.33e-6
    a_2 = 80.7
    b_2 = 4.42e-3
    c_2 = 1.37e-13
    d_2 = 651
    T_0 = 133
    
    # expresions auxiliars
    e_0 = (a_1 -b_1*T +c_1*(T^2) -d_1*(T^3))
    e_inf = (e_0 -a_2*np.exp(-b_2*T))
    tau = c_2 * np.exp(d_2/(T+T_0))
    
    return e_inf + ((e_0 - e_inf) / (1-(complex(0, 2*np.pi*tau))*f))

temperatures = [0,20,40,60,80,100]
frequencies = np.logspace(8, 12, num = 1000, base= 10)  

# plot de les 12 funcions
plt.figure()
for t in temperatures:
    plt.plot(frequencies, np.real(epsilon(frequencies, t)))
    plt.plot(frequencies, np.imag(epsilon(frequencies, t)))
    
plt.semilogx()
plt.xlabel('Freqüències (Hz)')
plt.ylabel("$\epsilon_r$ i $\epsilon_i$ ")    
plt.grid()
plt.savefig('eq-Debye.png')
plt.show()


###############################################################################
# APARTAT B) Calcul de les frequencies per les quals la permitivitat
# dielectrica imaginaria es maxima

for t in temperatures:
    imag_epsilon = np.imag(epsilon(frequencies, t))
    arg_max_imag_epsilon = np.argmax(imag_epsilon)
    maxim = frequencies[arg_max_imag_epsilon]
    print('Màxim de la part imaginària per a T=',t,' a la freqüencia:', maxim)

###############################################################################
# APARTAT C) Calcul de la frequencia per la qual la corba corresponent a
# la permitivitat dielectrica real es creua amb la imaginaria. (aixo es 
# calcula per a cada corba de diferent temperatura)

def function(f,T):
    '''Defineixo aquesta funció ja per poder determinar les interseccions 
    entre la part imaginaria i la part real de epsilon(f,T)'''
    
    a_1 = 87.9
    b_1 = 0.404
    c_1 = 9.59e-4
    d_1 = 1.33e-6
    a_2 = 80.7
    b_2 = 4.42e-3
    c_2 = 1.37e-13
    d_2 = 651
    T_0 = 133
    
    e_0 = (a_1 -b_1*T +c_1*(T^2) -d_1*(T^3))
    e_inf = (e_0 -a_2*np.exp(-b_2*T))
    tau = c_2 * np.exp(d_2/(T+T_0))
    
    return (e_inf*(1+(2*np.pi*f*tau)**2)/(e_0-e_inf)) +1 -2*np.pi*f*tau

data_f = []  # variable auxiliar per emmegatzemar les f on es tallen
for T in temperatures:
    frequencia_tall = opt.fsolve(function, 10**10, args=(T) )
    data_f.append(float(frequencia_tall))

# creo un DataFrame per mostra les dades en forma de taula com se'm demana    
dicc_data = {'Freqüencies de tall': data_f, 'Temperatures': temperatures}
df = pd.DataFrame(dicc_data)

print()
print(df)

###############################################################################
#APARTAT D) Frequencies de tall en funcio de la temperatura 
contador = 0
data_y = []   #variable auxilar on emmagatzem els valors de epsilon

for frequencia in data_f:
    y = np.real(epsilon(frequencia, temperatures[contador]))
    data_y.append(y)
    contador = contador +1 

result_lin_reg = stats.linregress(data_f, data_y)  
pendent = result_lin_reg[0]           # pendent
ordenada_origen = result_lin_reg[1]   # ordenada a l'origen
coef_correlacio = result_lin_reg[2]   # coeficient de correlacio
print()
print('Pendent de la recta de regressió: ', pendent)
print("Ordenada a l'origen de la recta de regressió:" , ordenada_origen)
print('Coeficient de Correlació:', coef_correlacio)

# representacio de les frequencies de tall en funcio de la temperatura
plt.figure()
plt.plot(data_f, temperatures)
plt.grid()
plt.xlabel('Freqüencies de tall (Hz)')
plt.ylabel('Temperatures (°C) ')
plt.savefig('T(f_tall).png')
plt.show()