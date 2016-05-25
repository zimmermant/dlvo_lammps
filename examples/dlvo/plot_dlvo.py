# Plotting DLVO force from pair_dlvo.cpp file with real constants
import numpy as np
import matplotlib.pyplot as plt

x_max=5e-6
points=10000


espr = 80 
epso =8.85E-12 
radius = 3e-6 
psi = 7e-2	 
psisqrd=psi**2
debyeinv = 104022291 
hamaker =1.3e-20 
X=np.linspace(0,x_max,points)
h=X-radius

f_dlvo = espr*epso*radius*psisqrd*debyeinv/(np.exp(-debyeinv*h)+1)*0.5 + radius*hamaker/(12*h**2)

#print f_dlvo
evdwl = espr*epso*radius*psisqrd*0.5*np.log(1+np.exp(-debyeinv*h))-radius*hamaker/(12*h)

plt.semilogy(X,f_dlvo)
plt.show()

plt.plot(X,evdwl)
plt.show()
