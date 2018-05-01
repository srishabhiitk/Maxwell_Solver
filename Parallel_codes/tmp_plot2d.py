import numpy as np
import h5py
import matplotlib.pyplot as plt


f=h5py.File('Ey_data_2d.h5','r')

i=140
key='t='+str(i)
data=f[key]
data=data[:,0,:]
data=data**2
data=data[450:550,0:100]
plt.imshow(data,animated=True,cmap='coolwarm')
plt.clim([0,0.05])
plt.axis('off')
plt.savefig("pml1.eps")
plt.show()
