#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
import matplotlib.pyplot as plt
import numpy as np
import numpy 
from matplotlib.colors import ListedColormap
from matplotlib.colors import BoundaryNorm
import lecture as lec
import matplotlib

#configuramos la fuente de matplotlib
font={'family':'serif',
	'size':5}
matplotlib.rc('font',**font)


#leemos los parametros
lec.rutine.parametros()

#lista de colores, maximo 10
COLORS=['blue','orange','red','green','yellow',  'gray','pink', 'brown', 'white', 'black']


#generamos una lista de colores segun la cantidad de resistividades
n_rho=np.int_(lec.rutine.n_rho)
color=range(n_rho)
color=COLORS[0:n_rho]
c=' '
c=ListedColormap(color)

#ajustamos los limites entre cada color
bounds=range(n_rho+1)
aux=-1.5
for i in range(n_rho+1):
	aux=aux+1.0
	bounds[i]=aux
norm = BoundaryNorm(bounds, c.N)


#creamos una imagen segun la cantidad de datos
n_data=np.int_(lec.rutine.n_data)
im=range(n_data)
fig, ax = plt.subplots(1)

#ajustar el grosor de los bordes del grafico
ax.spines['top'].set_linewidth(0.25)
ax.spines['right'].set_linewidth(0.25)
ax.spines['bottom'].set_linewidth(0.25)
ax.spines['left'].set_linewidth(0.25)
plt.yticks(lec.rutine.z)

for i in range(n_data):
	value=np.int_(lec.rutine.dat[i,0])
	ext=(lec.rutine.dat[i,1],lec.rutine.dat[i,2],lec.rutine.dat[i,3],lec.rutine.dat[i,4])
	im[i] = ax.imshow([[value]], extent=ext, picker=False, cmap=c, norm=norm, interpolation='nearest')

#ajustamos los limites del plot
limits=lec.rutine.limits
ax.axis([limits[0], limits[1], limits[3], limits[2]])

#generamos una imagen de dimension cero para crear la paleta de colores
val=range(n_rho)
img=ax.imshow([val], extent=(0, 0, 0, 0), picker=False, cmap=c, norm=norm, interpolation='nearest')
#fig.canvas.mpl_connect('pick_event', onpick4) 
cb=plt.colorbar(img, cmap=c, norm=norm, ticks=val, orientation='horizontal')
cb.outline.set_linewidth(0.5)

#escribimos las resitividades en las etiquetas
labels=range(n_rho)
for i in range(n_rho):
	labels[i]=str(lec.rutine.rho[i])
cb.set_ticklabels(labels)

fig=plt.gcf()
fig.set_size_inches(8.0,3.0)
plt.text(0.0,-0.5,'MODELO DE RESISTIVIDAD [Ohm-m]', fontdict={'size':7})
plt.savefig('modelo_resistividad.jpg',dpi=200)

