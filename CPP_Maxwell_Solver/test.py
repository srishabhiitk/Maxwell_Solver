import numpy as np
import h5py
from mayavi import mlab
from tvtk.util.ctf import ColorTransferFunction
from tvtk.util.ctf import PiecewiseFunction
filename = 'Ezdata2.h5'
f = h5py.File(filename, 'r')
filename = 'Ezdata.h5'
f_p = h5py.File(filename, 'r') #just plane wave

ctf = ColorTransferFunction()
ctf.add_rgb_point(1, 1, 0, 0)
# ctf.add_rgb_point(1, 0, 0, 1)


# Changing the otf:
otf = PiecewiseFunction()
# for i in range(101):
#     otf.add_point(i/100, np.sqrt(np.sqrt(i/100)))
# otf2 = PiecewiseFunction()
for i in range(101):
    otf.add_point(i/100.0, (i/100.0)**(0.7))

for i in [30,65,85]:
    mlab.figure(size=(800,799),bgcolor = (1,1,1))
    x,y,z=np.mgrid[0:99,0:99,0:100]
    data=np.array(f[f.keys()[i]]) #scattered wave
    data2=np.array(f_p[f_p.keys()[i]]) #plane wave
    data2[:,:,:]=0;
    # data2=np.abs(data2-data)
    # m=np.amax(data)
    m=np.amax(data)
    data[0,0,:]=m;
    data[0,-1,:]=m;
    data[-1,0,:]=m;
    data[-1,-1,:]=m;
    data[0,:,0]=m;
    data[0,:,-1]=m;
    data[-1,:,0]=m;
    data[-1,:,-1]=m;
    data[:,0,0]=m;
    data[:,0,-1]=m;
    data[:,-1,0]=m;
    data[:,-1,-1]=m;
    data2[45:56,45,45]=0.1;
    data2[45,45:56,45]=0.1;
    data2[45,45,45:56]=0.1;
    data2[45:56,55,55]=0.1;
    data2[55,45:56,55]=0.1;
    data2[55,55,45:56]=0.1;
    data2[45:56,45,55]=0.1;
    data2[45:56,55,45]=0.1;
    data2[45,45:56,55]=0.1;
    data2[55,45:56,45]=0.1;
    data2[45,55,45:56]=0.1;
    data2[55,45,45:56]=0.1;
    field=mlab.pipeline.scalar_field(x, y, z,(data))
    field2=mlab.pipeline.scalar_field(x, y, z,data2)
    vol=mlab.pipeline.volume(field)#, contours=10, transparent=True)
    vol._otf = otf
    vol._volume_property.set_scalar_opacity(otf)


    vol2=mlab.pipeline.volume(field2, vmax=0)
    vol2._volume_property.set_color(ctf)
    # # vol2._ctf = ctf
    # # vol2.update_ctf = True
    vol2._otf = otf
    vol2._volume_property.set_scalar_opacity(otf)

    mlab.view(azimuth=-60, elevation=80, distance=350, focalpoint=(49,49,49))
    # mlab.move(right=30)
    mlab.savefig(filename='./anim/fig%03d.png'%i)
    # mlab.show()
    mlab.close()


# mlab.show()
#ffmpeg -r 10 -i fig%03d.png -c:v libx264 -vf fps=10 -pix_fmt yuv420p out.mp4
#ffmpeg -f concat -i list.txt -c copy output.mp4
