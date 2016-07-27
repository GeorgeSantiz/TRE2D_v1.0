#!/usr/bin/env python
import numpy                        #modulo para manejo de arrays

#modulo para el trazo de graficas
#import matplotlib    
import matplotlib.tri as tri                
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap
#from matplotlib.ticker import MultipleLocator
#import matplotlib.figure as figure
#import matplotlib.font_manager as font_manager
import matplotlib.backends.backend_tkagg #sino se importa hay errores al compilar el *.py
import lecture as l
#from matplotlib.colors import BoundaryNorm

l.rutine.lres()
font={'family':'serif',
	'size':6}
matplotlib.rc('font',**font)

triang = tri.Triangulation(l.rutine.xr, l.rutine.yr)

# pcolor plot.
fig,ax=plt.subplots(1)
plt.subplots_adjust(left=0.1, right=0.75, bottom=None, top=0.80)
#plt.gca().set_aspect(1.0)

refiner = tri.UniformTriRefiner(triang)
tri_refi, z_test_refi = refiner.refine_field(l.rutine.ar, subdiv=4)


plt.tripcolor(tri_refi, z_test_refi, cmap=plt.get_cmap('jet'))#, norm=matplotlib.colors.LogNorm())
#plt.tripcolor(triang, l.rutine.ar)

plt.colorbar(orientation='horizontal', label='Ohm-m', extend='neither')
plt.clim(min(l.rutine.ar),max(l.rutine.ar))
#tic=numpy.arange(numpy.log10(min(r.res)), 0.3*10, )
#cb.set_ticks([numpy.log10(min(r.res)),numpy.log10(max(r.res)*0.999),numpy.log10(max(r.res))] )
#cb.set_ticklabels([min(r.res),max(r.res)*0.999,max(r.res)] )

plt.xlim(l.rutine.limits[0],l.rutine.limits[1])
plt.ylim(l.rutine.limits[3],l.rutine.limits[2])

plt.title('Pseudoseccion de Resistividad Aparente')
fig.set_size_inches(8.0,3.0)
fig=plt.gcf()
plt.savefig('resistividad.jpg',dpi=200)


