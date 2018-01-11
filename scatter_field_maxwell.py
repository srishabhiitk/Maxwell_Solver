import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



def Hx_aux_interp(i0,j0,i_start,i_end,j_start,j_end):
	i=np.arange(i_start,i_end)
	j=np.arange(j_start,j_end)
	# d=(i-i0)/np.sqrt(2)+(j-j0)/np.sqrt(2)
	d=(i-i0)/2.0+(j-j0+0.5)/2.0
	d=d-0.5
	d_int=d.astype(int)
	d=d-d_int
	# d=1-d
	return ((1-d)*H_aux[n+d_int+1]+d*H_aux[n+d_int])/np.sqrt(2)

def Hy_aux_interp(i0,j0,i_start,i_end,j_start,j_end):
	i=np.arange(i_start,i_end)
	j=np.arange(j_start,j_end)
	# d=(i-i0)/np.sqrt(2)+(j-j0)/np.sqrt(2)
	d=(i-i0+0.5)/2.0+(j-j0)/2.0
	d=d-0.5
	d_int=d.astype(int)
	d=d-d_int
	# d=1-d
	return -((1-d)*H_aux[n+d_int+1]+d*H_aux[n+d_int])/np.sqrt(2)

def E_aux_interp(i0,j0,i_start,i_end,j_start,j_end):
	i=np.arange(i_start,i_end)
	j=np.arange(j_start,j_end)
	# d=(i-i0)/np.sqrt(2)+(j-j0)/np.sqrt(2)
	d=(i-i0)/2.0+(j-j0)/2.0
	# d=d+0.5
	d_int=d.astype(int)
	d=d-d_int
	# d=1-d
	return (1-d)*E_aux[n+d_int+1]+d*E_aux[n+d_int]











num_timesteps=1200
n=1000
m=400 #Total field region mXm

epsilon0=8.85*(10**(-12))
# epsilon=epsilon0*np.ones([n,n])
mew0=8.85*(10**(-12))
# mew=mew0*np.ones([n,n])
# epsilon[400:600,400:600]=epsilon[400:600,400:600]*2
# mew[400:600,400:600]=mew[400:600,400:600]*2

Ez=np.zeros([n,n])
Hy=np.zeros([n-1,n])
Hx=np.zeros([n,n-1])
E_aux=np.zeros(2*n)
H_aux=np.zeros(2*n-1)

c=1/np.sqrt(epsilon0*mew0)
S=1/(2**0.5)
dx=1.0*(10**(-6))
dt=S*dx/c
lambd=20.0*dx
N_lamdb=lambd/dx
vp_ratio=np.sqrt(2)*np.arcsin(np.sin(np.pi*S/N_lamdb)/S/np.sqrt(2))/np.arcsin(np.sin(np.pi*S/N_lamdb)/S)
vp_ratio=(1/vp_ratio)
omega=2*np.pi*c/lambd

Ca=Da=1.0
Cb=dt/(epsilon0*dx)
Db=dt/(mew0*dx)
Cb_aux=dt/(epsilon0*dx*np.sqrt(2))*vp_ratio
Db_aux=dt/(mew0*dx*np.sqrt(2))*vp_ratio








p=int((n+m)/2)
q=int((n-m)/2)
fig = plt.figure()
imz=[]
for timesteps in range(num_timesteps+1):
	if (timesteps==0):
		# Ez[int(n/2),int(n/2)]=1
		# E_aux[n]=np.sin(0.01*omega*timesteps*dt)
		# E_aux[n]=1
		E_aux[n]=np.exp(-((timesteps-300)**2)/400)

	H_aux=H_aux+Db_aux*(E_aux[0:2*n-1]-E_aux[1:2*n])
	E_aux[1:2*n-1]=E_aux[1:2*n-1]+Cb_aux*(H_aux[0:2*n-2]-H_aux[1:2*n-1])



	Ez[450:550,450:550]=0

	# Hy=Da*Hy+Db*(Ez[1:n,0:n]-Ez[0:n-1,0:n])
	Hy[0:n-1,0:q]=Da*Hy[0:n-1,0:q]+Db*(Ez[1:n,0:q]-Ez[0:n-1,0:q]) #Scatter top
	Hy[0:n-1,p:n]=Da*Hy[0:n-1,p:n]+Db*(Ez[1:n,p:n]-Ez[0:n-1,p:n]) #Scatter bottom
	Hy[0:q-1,q:p]=Da*Hy[0:q-1,q:p]+Db*(Ez[1:q,q:p]-Ez[0:q-1,q:p]) #Scatter left
	Hy[q-1,q:p]=Da*Hy[q-1,q:p]+Db*(Ez[q,q:p]-Ez[q-1,q:p])-Db*E_aux_interp(q,q,q,q+1,q,p) #Special left
	Hy[q:p-1,q:p]=Da*Hy[q:p-1,q:p]+Db*(Ez[q+1:p,q:p]-Ez[q:p-1,q:p]) # Total
	Hy[p-1,q:p]=Da*Hy[p-1,q:p]+Db*(Ez[p,q:p]-Ez[p-1,q:p])+Db*E_aux_interp(q,q,p-1,p,q,p) #Special right
	Hy[p:n-1,q:p]=Da*Hy[p:n-1,q:p]+Db*(Ez[p+1:n,q:p]-Ez[p:n-1,q:p]) #Scatter right

	# Hx=Da*Hx-Db*(Ez[0:n,1:n]-Ez[0:n,0:n-1])
	Hx[0:q,0:n-1]=Da*Hx[0:q,0:n-1]-Db*(Ez[0:q,1:n]-Ez[0:q,0:n-1]) #Scatter left
	Hx[p:n,0:n-1]=Da*Hx[p:n,0:n-1]-Db*(Ez[p:n,1:n]-Ez[p:n,0:n-1]) #Scatter right
	Hx[q:p,0:q-1]=Da*Hx[q:p,0:q-1]-Db*(Ez[q:p,1:q]-Ez[q:p,0:q-1]) #Scatter top
	Hx[q:p,q-1]=Da*Hx[q:p,q-1]-Db*(Ez[q:p,q]-Ez[q:p,q-1])+Db*E_aux_interp(q,q,q,p,q,q+1) #Special top
	Hx[q:p,q:p-1]=Da*Hx[q:p,q:p-1]-Db*(Ez[q:p,q+1:p]-Ez[q:p,q:p-1]) # Total
	Hx[q:p,p-1]=Da*Hx[q:p,p-1]-Db*(Ez[q:p,p]-Ez[q:p,p-1])-Db*E_aux_interp(q,q,q,p,p-1,p) #Special Bottom
	Hx[q:p,p:n-1]=Da*Hx[q:p,p:n-1]-Db*(Ez[q:p,p+1:n]-Ez[q:p,p:n-1]) #Scatter bottom

	# Ez[1:n-1,1:n-1]=Ca*Ez[1:n-1,1:n-1]+Cb*(Hy[1:n-1,1:n-1]-Hy[0:n-2,1:n-1]+Hx[1:n-1,0:n-2]-Hx[1:n-1,1:n-1])
	Ez[1:n-1,1:q]=Ca*Ez[1:n-1,1:q]+Cb*(Hy[1:n-1,1:q]-Hy[0:n-2,1:q]+Hx[1:n-1,0:q-1]-Hx[1:n-1,1:q]) #Scatter top
	Ez[1:n-1,p:n-1]=Ca*Ez[1:n-1,p:n-1]+Cb*(Hy[1:n-1,p:n-1]-Hy[0:n-2,p:n-1]+Hx[1:n-1,p-1:n-2]-Hx[1:n-1,p:n-1]) #Scatter bottom
	Ez[1:q,q:p]=Ca*Ez[1:q,q:p]+Cb*(Hy[1:q,q:p]-Hy[0:q-1,q:p]+Hx[1:q,q-1:p-1]-Hx[1:q,q:p]) #Scatter left
	Ez[p:n-1,q:p]=Ca*Ez[p:n-1,q:p]+Cb*(Hy[p:n-1,q:p]-Hy[p-1:n-2,q:p]+Hx[p:n-1,q-1:p-1]-Hx[p:n-1,q:p]) #Scatter right
	Ez[q+1:p-1,q+1:p-1]=Ca*Ez[q+1:p-1,q+1:p-1]+Cb*(Hy[q+1:p-1,q+1:p-1]-Hy[q:p-2,q+1:p-1]+Hx[q+1:p-1,q:p-2]-Hx[q+1:p-1,q+1:p-1]) #Total
	Ez[q+1:p-1,q]=Ca*Ez[q+1:p-1,q]+Cb*(Hy[q+1:p-1,q]-Hy[q:p-2,q]+Hx[q+1:p-1,q-1]-Hx[q+1:p-1,q])+Cb*Hx_aux_interp(q,q,q+1,p-1,q-1,q) #Special top edge
	Ez[q+1:p-1,p-1]=Ca*Ez[q+1:p-1,p-1]+Cb*(Hy[q+1:p-1,p-1]-Hy[q:p-2,p-1]+Hx[q+1:p-1,p-2]-Hx[q+1:p-1,p-1])-Cb*Hx_aux_interp(q,q,q+1,p-1,p-1,p) #Special bottom edge
	Ez[q,q+1:p-1]=Ca*Ez[q,q+1:p-1]+Cb*(Hy[q,q+1:p-1]-Hy[q-1,q+1:p-1]+Hx[q,q:p-2]-Hx[q,q+1:p-1])-Cb*Hy_aux_interp(q,q,q-1,q,q+1,p-1) #Special left edge
	Ez[p-1,q+1:p-1]=Ca*Ez[p-1,q+1:p-1]+Cb*(Hy[p-1,q+1:p-1]-Hy[p-2,q+1:p-1]+Hx[p-1,q:p-2]-Hx[p-1,q+1:p-1])+Cb*Hy_aux_interp(q,q,p-1,p,q+1,p-1) #Special right edge
	Ez[q,q]=Ca*Ez[q,q]+Cb*(Hy[q,q]-Hy[q-1,q]+Hx[q,q-1]-Hx[q,q])-Cb*Hy_aux_interp(q,q,q-1,q,q,q+1)+Cb*Hx_aux_interp(q,q,q,q+1,q-1,q) #Special top-left
	Ez[p-1,q]=Ca*Ez[p-1,q]+Cb*(Hy[p-1,q]-Hy[p-2,q]+Hx[p-1,q-1]-Hx[p-1,q])+Cb*Hy_aux_interp(q,q,p-1,p,q,q+1)+Cb*Hx_aux_interp(q,q,p-1,p,q-1,q) #Special top-right
	Ez[q,p-1]=Ca*Ez[q,p-1]+Cb*(Hy[q,p-1]-Hy[q-1,p-1]+Hx[q,p-2]-Hx[q,p-1])-Cb*Hy_aux_interp(q,q,q-1,q,p-1,p)-Cb*Hx_aux_interp(q,q,q,q+1,p-1,p) #Special bottom-left
	Ez[p-1,p-1]=Ca*Ez[p-1,p-1]+Cb*(Hy[p-1,p-1]-Hy[p-2,p-1]+Hx[p-1,p-2]-Hx[p-1,p-1])+Cb*Hy_aux_interp(q,q,p-1,p,p-1,p)-Cb*Hx_aux_interp(q,q,p-1,p,p-1,p) #Special bottom-right


	Ez[450:550,450:550]=0

	# Ez[int(n/2),int(n/2)]=1
	# E_aux[n]=np.sin(0.01*omega*timesteps*dt)
	E_aux[n]=np.exp(-((timesteps-300)**2)/400)
	# E_aux[n]=1
	print (timesteps)

	if(np.mod(timesteps,5)==0):
		Ez[450:550,450:550]=1
		frame=plt.imshow(np.abs(Ez)**2,animated=True,cmap='gray')
		plt.clim([0,1])
		imz.append([frame])





ani = animation.ArtistAnimation(fig, imz, interval=50, blit=True,repeat_delay=1000)
plt.show()
# plt.plot(E_aux)
# plt.show()
