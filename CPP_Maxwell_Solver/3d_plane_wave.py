import numpy as np
import h5py
from mayavi import mlab
# filename = 'Ezdata.h5'
# f = h5py.File(filename, 'r')


x,y,z=np.mgrid[0:99,0:99,0:100]
dx=10**(-3)
data=np.exp(-(((x*dx+y*dx+z*dx)/(3**0.5)+00*dx-100*dx/(3**0.5))**2)/(10**(-3)));
# data=np.array(f['t10'])
field=mlab.pipeline.scalar_field(x, y, z,(data))
mlab.pipeline.volume(field)#, contours=10, transparent=True)
# # mlab.savefig(filename='fig'+str(i)+'.png')
mlab.show()
