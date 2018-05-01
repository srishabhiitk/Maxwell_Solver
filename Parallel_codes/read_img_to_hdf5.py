import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import misc


img=misc.imread('wg2.png')
img=img[:,:,0]
img=img/255
img=img*1.0
img=1-img
img2=np.zeros([1002,1,1002])
img2[:,0,:]=img


f=h5py.File('wg2.h5','w')
f.create_dataset('data',data=img2)
f.close();

plt.imshow(img,cmap='gray')
plt.show()
