import numpy as np
import matplotlib.pyplot as plt

xdim=200
ydim=200


time_tot=100

xsource=100
ysource=100


S=1/(2**0.5)


epsilon0=(1/(36*np.pi))*(10**-(9))
mu0=4*np.pi*(10**(-7))
c=3*(10**8)



delta=10**(-6)
deltat=S*delta/c

Ez=np.zeros([xdim,ydim])
Hy=np.zeros([xdim,ydim])
Hx=np.zeros([xdim,ydim])


epsilon=epsilon0
mu=mu0


sigma=4*(10**(-4))
sigma_star=4*(10**(-4))


A=((mu-0.5*deltat*sigma_star)/(mu+0.5*deltat*sigma_star))
B=(deltat/delta)/(mu+0.5*deltat*sigma_star)


C=((epsilon-0.5*deltat*sigma)/(epsilon+0.5*deltat*sigma))
D=(deltat/delta)/(epsilon+0.5*deltat*sigma)


for n in np.arange(1,time_tot+1):



    n1=0
    n2=xdim-1
    n11=0
    n21=ydim-1




    Hy[n1:n2,n11:n21]=A*Hy[n1:n2,n11:n21]+B*(Ez[n1+1:n2+1,n11:n21]-Ez[n1:n2,n11:n21])
    Hx[n1:n2,n11:n21]=A*Hx[n1:n2,n11:n21]-B*(Ez[n1:n2,n11+1:n21+1]-Ez[n1:n2,n11:n21])

    Ez[n1+1:n2+1,n11+1:n21+1]=C*Ez[n1+1:n2+1,n11+1:n21+1]+D*(Hy[n1+1:n2+1,n11+1:n21+1]-Hy[n1:n2,n11+1:n21+1]-Hx[n1+1:n2+1,n11+1:n21+1]+Hx[n1+1:n2+1,n11:n21])






    Ez[0:xdim-1,0]=0
    Ez[0:xdim-1,ydim-1]=0
    Ez[0,0:ydim-1]=0
    Ez[xdim-1,0:ydim-1]=0

    Ez[xsource,ysource]=1

    print(n)


plt.matshow(Ez)
plt.colorbar()
plt.show()
