import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import png



num_timesteps=2500
n=1000

epsilon=8.85*(10**(-12))
mew=4*np.pi*(10**(-7))

Ez=np.zeros([n,n])
Hy=np.zeros([n-1,n])
Hx=np.zeros([n,n-1])

x,y=np.meshgrid(np.arange(0,n),np.arange(0,n))
# x=x[0]
# x=np.array(x)

c=1/np.sqrt(epsilon*mew)
S=1/(2**0.5)
dx=1.0*(10**(-6))
dt=S*dx/c
lambd=20.0*dx
omega=2*np.pi*c/lambd

Ca=Da=1.0
Cb=np.ones_like(Ez)
Cb=Cb*dt/(epsilon*dx)
Db=dt/(mew*dx)

x=x*(dx)
y=y*dx






Cb[int(n/2)-50:int(n/2)+50,int(n/2)-50:int(n/2)+50]=Cb[int(n/2)-50:int(n/2)+50,int(n/2)-50:int(n/2)+50]/10
fig = plt.figure()
imz=[]
for timesteps in range(num_timesteps):
	# if (timesteps==0):
	# 	Ez[int(n/2),int(n/2)]=1

	# Ez[int(n/2)-50:int(n/2)+50,int(n/2)-50:int(n/2)+50]=0
	Ez[0,0:n]=np.exp(-(((x[0,0:n]+y[0,0:n])/np.sqrt(2)+300*dx-timesteps*c*dt)**2)/0.000000001)
	Ez[n-1,0:n]=np.exp(-(((x[n-1,0:n]+y[n-1,0:n])/np.sqrt(2)+300*dx-timesteps*c*dt)**2)/0.000000001)
	Ez[1:n-1,0]=np.exp(-(((x[1:n-1,0]+y[1:n-1,0])/np.sqrt(2)+300*dx-timesteps*c*dt)**2)/0.000000001)
	Ez[1:n-1,n-1]=np.exp(-(((x[1:n-1,n-1]+y[1:n-1,n-1])/np.sqrt(2)+300*dx-timesteps*c*dt)**2)/0.000000001)


	Hy=Da*Hy+Db*(Ez[1:n,0:n]-Ez[0:n-1,0:n])
	Hx=Da*Hx-Db*(Ez[0:n,1:n]-Ez[0:n,0:n-1])
	Ez[1:n-1,1:n-1]=Ca*Ez[1:n-1,1:n-1]+Cb[1:n-1,1:n-1]*(Hy[1:n-1,1:n-1]-Hy[0:n-2,1:n-1]+Hx[1:n-1,0:n-2]-Hx[1:n-1,1:n-1])

	# Ez[int(n/2)-50:int(n/2)+50,int(n/2)-50:int(n/2)+50]=0

	# Ez[int(n/2),int(n/2)]=1
	print (timesteps, np.amax(Ez))

	if(np.mod(timesteps,5)==0):
		# Ez[int(n/2)-50:int(n/2)+50,int(n/2)-50:int(n/2)+50]=1
		gg=Ez**2
		gg[500-50,500-50:500+50+1]=0.5
		gg[500+50,500-50:500+50+1]=0.5
		gg[500-50:500+50+1,500-50]=0.5
		gg[500-50:500+50+1,500+50]=0.5
		frame=plt.imshow(gg,animated=True,cmap='coolwarm')
		plt.clim([0,1])
		plt.axis('off')
		imz.append([frame])
	# if (timesteps==1000 or timesteps==1745 or timesteps==2245):
	# 	filegg='out_'+str(timesteps)
	# 	np.save(filegg,Ez**2)
		# plt.imshow(Ez**2,cmap='GnBu')
		# plt.savefig(filegg)






ani = animation.ArtistAnimation(fig, imz, interval=50, blit=True,repeat_delay=1000)
plt.show()
