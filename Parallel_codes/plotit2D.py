import numpy as np
import h5py
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 10, 10


fig = plt.figure()
imz=[]

f=h5py.File('Ey_data_2d.h5','r')
# rf=h5py.File('wg2.h5','r')
# crystal=rf['data']
# crystal=crystal[:,0,:]
# crystal=np.gradient(crystal)
# crystal=np.sqrt(crystal[0]**2+crystal[1]**2)

for i in np.arange(100,2000,5):
    key='t='+str(i)
    data=f[key]
    data=data[:,0,:]
    # data=data+crystal*0.5
    data=data**2
    # data=data[450:550,0:100]
    frame=plt.imshow(data,animated=True,cmap='coolwarm')
    # plt.clim([0,1])
    plt.axis('off')
    plt.tight_layout()
    imz.append([frame])


ani = animation.ArtistAnimation(fig, imz, interval=50, blit=True,repeat_delay=1000)
plt.show()
