import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



num_timesteps=1200
n=1000

epsilon=8.85*(10**(-12))
mew=4*np.pi*(10**(-7))

Ez=np.zeros([n,n])
Hy=np.zeros([n-1,n])
Hx=np.zeros([n,n-1])

c=1/np.sqrt(epsilon*mew)
S=1/(2**0.5)
dx=1.0*(10**(-6))
dt=S*dx/c
lambd=20.0*dx
omega=2*np.pi*c/lambd

Ca=Da=1.0
Cb=dt/(epsilon*dx)
Db=dt/(mew*dx)









fig = plt.figure()
imz=[]
for timesteps in range(num_timesteps):
	if (timesteps==0):
		Ez[int(n/2),int(n/2)]=1

	Hy=Da*Hy+Db*(Ez[1:n,0:n]-Ez[0:n-1,0:n])
	Hx=Da*Hx-Db*(Ez[0:n,1:n]-Ez[0:n,0:n-1])
	Ez[1:n-1,1:n-1]=Ca*Ez[1:n-1,1:n-1]+Cb*(Hy[1:n-1,1:n-1]-Hy[0:n-2,1:n-1]+Hx[1:n-1,0:n-2]-Hx[1:n-1,1:n-1])

	Ez[int(n/2),int(n/2)]=1
	print (timesteps)

	if(np.mod(timesteps,5)==0):
		frame=plt.imshow(Ez,animated=True)
		imz.append([frame])





ani = animation.ArtistAnimation(fig, imz, interval=50, blit=True,repeat_delay=1000)
plt.show()
