
import sys
import numpy as np
import matplotlib.pyplot as plt
 
###############################################################################
# CONJUNT DE MANDELBROT

# apartat 1 i 2, crear la matriu de valors de c
cx = np.linspace(-2,1,1000)
cy = np.linspace(-1.5,1.5,1000)
x, y = np.meshgrid(cx, cy)   # coordenades en x i en y

# la matriu conjunta de valors de c
c = x+1j*y

print("Valor maxim del sistema: ",sys.float_info.max)  # valor màxim
fita_superior = sys.float_info.max

# calcul conjunt de Manderlot
z_i = 0
vc_i = 0
num_simulacions = 0

plt.figure()
while num_simulacions <= 35:
       
    # Conjunt de Manderlot
    z = z_i**2 + c
    control_M = np.abs(z) < 1e10    
    z = z*control_M + (1-control_M)*1e10
    cm = np.abs(z) < 2
    
    z_i = z 
    
    # Velocitat de Convergencia (VC)
    vc = 1*cm + vc_i
    vc_i = vc
    
    if num_simulacions == 1:
        plt.subplot(2,2,1)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')

        
    if num_simulacions == 3:
        plt.subplot(2,2,2)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')

        
    if num_simulacions == 6:
        plt.subplot(2,2,3)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')

        
    if num_simulacions == 10:
        plt.subplot(2,2,4)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')

        
    num_simulacions = num_simulacions + 1

plt.suptitle('CM Mandelbrot')
plt.savefig('Mandelbrot-cm-variable.png')
plt.show()

plt.figure()
plt.imshow(cm, cmap = 'gray')
plt.axis('off')
plt.title('CM final Mandelbrot')
plt.savefig('Mandelbrot-cm-final.png')
plt.show()

plt.figure()
plt.imshow(vc, cmap = 'jet')
plt.axis('off')
plt.title('VM Mandelbrot')
plt.savefig('Mandelbrot-vm-final.png')
plt.show()

###############################################################################
# CONJUNT DE JULIA

# apartat 1 i 2, crear la matriu de valors de c
cx = np.linspace(-2,2,1000)
cy = np.linspace(-2,2,1000)
x, y = np.meshgrid(cx, cy)   # coordenades en x i en y

# la matriu conjunta de valors de c
c = x+1j*y

fita_superior = sys.float_info.max

# calcul conjunt de Julia
z_i = c
vc_i = 0
num_simulacions = 0

plt.figure()
while num_simulacions <= 100:
       
    # Conjunt de Julia
    z = z_i**2 -0.7269 + 0.1889*(1j)
    control_M = np.abs(z) < 1e10    
    z = z*control_M + (1-control_M)*1e10
    cm = np.abs(z) < 2
    
    z_i = z 
    
    # Velocitat de Convergencia (VC)
    vc = 1*cm + vc_i
    vc_i = vc
    
    if num_simulacions == 1:
        plt.subplot(2,2,1)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')
        
    if num_simulacions == 3:
        plt.subplot(2,2,2)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')
        
    if num_simulacions == 6:
        plt.subplot(2,2,3)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')
        
    if num_simulacions == 10:
        plt.subplot(2,2,4)
        plt.imshow(cm, cmap = 'gray')
        plt.axis('off')
        
    num_simulacions = num_simulacions + 1

plt.suptitle('CM Julia')
plt.tight_layout()
plt.savefig('Julia-cm-variable.png')
plt.show()

plt.figure()
plt.imshow(cm, cmap = 'gray')
plt.title('CM final Julia')
plt.axis('off')
plt.show()

plt.figure()
plt.imshow(vc, cmap = 'jet')
plt.title('VC Julia')
plt.axis('off')
plt.show()