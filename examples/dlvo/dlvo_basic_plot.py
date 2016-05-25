# Plotting DLVO force from pair_dlvo.cpp file with real constants
import numpy as np
import matplotlib.pyplot as plt

x_max=1.2
points=10000


espr = 1 
epso =1
radius =1 
psi =1	 
psisqrd=psi**2
debyeinv = 1
hamaker =1
X=np.linspace(0,x_max,points)
h=X-radius

f_dlvo = espr*epso*radius*psisqrd*debyeinv/(np.exp(-debyeinv*h)+1)*0.5 + radius*hamaker/(12*h**2)

#print f_dlvo
evdwl = espr*epso*radius*psisqrd*0.5*np.log(1+np.exp(-debyeinv*h))-radius*hamaker/(12*h)

plt.semilogy(X,f_dlvo)
plt.show()

plt.semilogy(X,evdwl)
plt.show()
