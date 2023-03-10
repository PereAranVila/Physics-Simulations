
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Index de refraccio en funcio de la longitud d'ona n(lambda) pel 
# Ti02, MgFe2 i BK7

def index_refraccio_TiO2(longitud_ona):
    '''longituds ona entre 0.43-1.53 micrometres'''
    
    n = ((5.913 + (0.2411/(longitud_ona**2 - 0.0803))))**0.5
    return n

def index_refraccio_MgF2(longitud_ona):
    '''lognituds d'ona entre 0.2-7 micrometres'''
    
    n=((1+0.48755108/(1-(0.04338408/longitud_ona)**2)
        +0.39875031/(1-(0.09461442/longitud_ona)**2)+
        2.3120353/(1-(23.793604/longitud_ona)**2))**0.5)
    return n

def index_refraccio_BK7(longitud_ona):
    '''Es vàlida l'expressió per logntiuds d'ona entre 0.3-2.5 mirometres'''
    
    n=((1+1.03961212/(1-0.00600069867/longitud_ona**2)
        +0.231792344/(1-0.0200179144/longitud_ona**2)
        +1.01046945/(1-103.560653/longitud_ona**2))**0.5)
    return n

# calculs dels n(\lambda)
longituds_ona_tio2 = np.linspace(0.43, 1.53, 1000)
indexs_refrac_tio2 = index_refraccio_TiO2(longituds_ona_tio2)
longituds_ona_mgf2 = np.linspace(0.2, 7, 1000)
indexs_refra_mgf2 = index_refraccio_MgF2(longituds_ona_mgf2)
longituds_ona_bk7 = np.linspace(0.3, 2.5, 1000)
indexs_refrac_bk7 = index_refraccio_BK7(longituds_ona_bk7)

# representacio n(lambda)
figure = plt.figure(num = 0)
plt.subplot(2,2,1)
plt.plot(longituds_ona_tio2, indexs_refrac_tio2)
plt.title('$TiO_{2}$')
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('n')
plt.subplot(2,2,2)
plt.plot(longituds_ona_mgf2, indexs_refra_mgf2)
plt.title('$MgF_{2}$')
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('n')
plt.subplot(2,2,3)
plt.plot(longituds_ona_bk7, indexs_refrac_bk7)
plt.title('$BK_{7}$')
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('n')
plt.tight_layout()
plt.savefig('n-TiO2-MgFe2-BK7.png')
plt.show()

# representaccio conjunta de n(lambda)
plt.figure()
plt.plot(longituds_ona_tio2, indexs_refrac_tio2,label='$TiO_{2}$')
plt.plot(longituds_ona_mgf2, indexs_refra_mgf2,label='$MgF_{2}$')
plt.plot(longituds_ona_bk7, indexs_refrac_bk7,label='$BK_{7}$')
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('n')
plt.legend()
plt.savefig('n-conjunt.png')
plt.show()


###############################################################################

'''Funció que crea la matriu M'''

def matrix_M(n, d, longitud_ona):

    '''Matriu M per incidencia normal'''
    
    M = np.array([[np.cos((2*np.pi*n*d)/longitud_ona), 
                (complex(0+1j)*np.sin((2*np.pi*n*d)/longitud_ona))],
                [complex(0+1j)*n*np.sin((2*np.pi*n*d)/longitud_ona),
                np.cos((2*np.pi*n*d)/longitud_ona)]])
    
    return M


########################
# Simulacio de dues capes de TiO2 i MgF2 de gruix lambda/4

longitud_referencia = 0.550e-6  # longitud d'ona = 550 nm 
n_tio2_ref = index_refraccio_TiO2(longitud_referencia)
n_mgf2_ref = index_refraccio_MgF2(longitud_referencia)
# gruixos de les dues capes
d_tio2 = longitud_referencia/(4*n_tio2_ref)
d_mgf2 = longitud_referencia/(4*n_mgf2_ref)

# matriu de transferencia total
longituds_ona = np.linspace(0.43, 1.53, 1000)
indexs_refrac_tio2 = index_refraccio_TiO2(longituds_ona)
indexs_refrac_mgf2 = index_refraccio_MgF2(longituds_ona)
n_0 = 1
n_s = 1.52
R_array = []
T_array = []
for j in range (0,1000):
    n_tio2 = indexs_refrac_tio2[j]
    n_mgf2 = indexs_refrac_mgf2[j]
    matrix_tio2 = matrix_M(n_tio2, d_tio2, longitud_referencia)
    matrix_mgf2 = matrix_M(n_mgf2, d_tio2, longitud_referencia)


    matrix_transf_total = matrix_tio2 @ matrix_mgf2

    # valors B i C
    matrix_n_0_s = np.array([[n_0], [n_s]])
    matrix_B_C = matrix_transf_total @ matrix_n_0_s
    B = matrix_B_C[0][0]
    C = matrix_B_C[1][0]

    R = np.abs((n_0*B-C)/(n_0*B+C))**2
    T = 1-R
    R_array.append(R)
    T_array.append(T)

# representem T(lambda) i R(lambda)
plt.figure()
plt.suptitle('$TiO_{2}$ i $MgFe_{2}$ amb gruix $\lambda$ /4')
plt.subplot(2,1,1)
plt.plot(longituds_ona, R_array)
plt.xlabel('$\lambda (\mu m)$')
plt.xlim(0.4, 1.6)
plt.ylim(0.16, 0.20)
plt.ylabel('R')
plt.subplot(2,1,2)
plt.plot(longituds_ona, T_array)
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('T')
plt.xlim(0.4, 1.6)
plt.savefig('TiO2-MgFe2.png')
plt.show()


########################
# Sistema de 5 capes de TiO2 i MgF2 sobre vidre BK7

longituds_ona = np.linspace(0.43, 1.53, 1000)
indexs_refrac_tio2 = index_refraccio_TiO2(longituds_ona)
indexs_refrac_mgf2 = index_refraccio_MgF2(longituds_ona)
n_0 = 1
n_s = 1.52
R_array = []
T_array = []

for j in range(0,1000):
    n_tio2 = indexs_refrac_tio2[j]
    n_mgf2 = indexs_refrac_mgf2[j]
    matrix_tio2 = matrix_M(n_tio2, d_tio2, longitud_referencia)
    matrix_mgf2 = matrix_M(n_mgf2, d_tio2, longitud_referencia)
    
    matrix_transf_total = (matrix_tio2 @ matrix_mgf2 @ matrix_tio2 @ 
                matrix_mgf2 @ matrix_tio2 @ matrix_mgf2 @ matrix_tio2 @ 
                matrix_mgf2 @ matrix_tio2 @ matrix_mgf2)
    
    
    # valors B i C
    matrix_n_0_s = np.array([[n_0], [n_s]])
    matrix_B_C = matrix_transf_total @ matrix_n_0_s
    B = matrix_B_C[0][0]
    C = matrix_B_C[1][0]

    R = np.abs((n_0*B-C)/(n_0*B+C))**2
    T = 1-R
    R_array.append(R)
    T_array.append(T)

# representacio T(lambda) i R(lambda)
plt.figure()
plt.suptitle('$TiO_{2}$ i $MgFe_{2}$ 5 parell de capes amb gruix '+
         '$\lambda$ /4')
plt.subplot(2,1,1)
plt.plot(longituds_ona, R_array)
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('R')
plt.subplot(2,1,2)
plt.plot(longituds_ona, T_array)
plt.xlabel('$\lambda (\mu m)$')
plt.ylabel('T')
plt.savefig('5-parells-TiO2-MgFe2.png')
plt.show()

