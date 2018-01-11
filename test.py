import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation





gg=np.load('out_2245.npy')
gg[500-50,500-50:500+50+1]=np.amax(gg)/2
gg[500+50,500-50:500+50+1]=np.amax(gg)/2
gg[500-50:500+50+1,500-50]=np.amax(gg)/2
gg[500-50:500+50+1,500+50]=np.amax(gg)/2
plt.imshow(gg,cmap='coolwarm')

plt.colorbar()

plt.savefig('out3.eps')
plt.clf()


gg=np.load('out_1745.npy')
gg[500-50,500-50:500+50+1]=np.amax(gg)/2
gg[500+50,500-50:500+50+1]=np.amax(gg)/2
gg[500-50:500+50+1,500-50]=np.amax(gg)/2
gg[500-50:500+50+1,500+50]=np.amax(gg)/2
plt.imshow(gg,cmap='coolwarm')

plt.colorbar()

plt.savefig('out2.eps')
plt.clf()

gg=np.load('out_1000.npy')
gg[500-50,500-50:500+50+1]=np.amax(gg)/2
gg[500+50,500-50:500+50+1]=np.amax(gg)/2
gg[500-50:500+50+1,500-50]=np.amax(gg)/2
gg[500-50:500+50+1,500+50]=np.amax(gg)/2
plt.imshow(gg,cmap='coolwarm')

plt.colorbar()

plt.savefig('out1.eps')
plt.clf()
