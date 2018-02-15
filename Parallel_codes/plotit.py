import numpy as np
import h5py
from mayavi import mlab
from tvtk.util.ctf import ColorTransferFunction
from tvtk.util.ctf import PiecewiseFunction
filename = 'Ez_data.h5'
f = h5py.File(filename, 'r')
# filename = 'Ezdata.h5'
# f_p = h5py.File(filename, 'r') #just plane wave


# mlab.options.offscreen = True
# List all groups
# print("Keys: %s" % f.keys())
# a_group_key = f.keys()[0]

# Get the data
#
# for i in range(200):
#     if (np.mod(i,10)==0):
#         data = np.array(f['t'+str(i)])
#         mlab.contour3d(data, contours=10, transparent=True)
#         mlab.savefig(filename='fig'+str(i)+'.png')



ctf = ColorTransferFunction()
# ctf.add_rgb_point(0, 0, 0, 1)
# ctf.add_rgb_point(200, 0, 0, 1)
# ctf.add_rgb_point(499, 1, 0, 0)
# ctf.add_rgb_point(501, 1, 0, 0)
ctf.add_rgb_point(1, 1, 0, 0)


# Changing the otf:
otf = PiecewiseFunction()
# otf.add_point(502, 1)
for i in range(101):
    otf.add_point(i/100.0, (np.sqrt(i/100.0)))
otf2 = PiecewiseFunction()
for i in range(101):
    otf.add_point(i/100, np.sqrt(np.sqrt(np.sqrt(np.sqrt(i/100)))))

for i in range(len(f.keys())):
    mlab.figure(size=(800,799),bgcolor = (1,1,1))
    x,y,z=np.mgrid[0:102,0:102,0:101]
    data=np.array(f[f.keys()[i]]) #scattered wave
    # data2=np.array(f_p[f_p.keys()[i]]) #plane wave
    # data[:,:,:]=0
    m=np.amax(data)
    # m=0.5;
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
    # m1=0.1;
    # data2[45:56,45,45]=m1;
    # data2[45,45:56,45]=m1;
    # data2[45,45,45:56]=m1;
    # data2[45:56,55,55]=m1;
    # data2[55,45:56,55]=m1;
    # data2[55,55,45:56]=m1;
    # data2[45:56,45,55]=m1;
    # data2[45:56,55,45]=m1;
    # data2[45,45:56,55]=m1;
    # data2[55,45:56,45]=m1;
    # data2[45,55,45:56]=m1;
    # data2[55,45,45:56]=m1;
    field=mlab.pipeline.scalar_field(x, y, z,(data))
    # field2=mlab.pipeline.scalar_field(x, y, z,data2)
    vol=mlab.pipeline.volume(field)#, contours=10, transparent=True)
    # vol._ctf = ctf
    # vol.update_ctf = True
    vol._otf = otf
    vol._volume_property.set_scalar_opacity(otf)


    # vol2=mlab.pipeline.volume(field2, vmax=0)
    # vol2._volume_property.set_color(ctf)
    # # vol2._ctf = ctf
    # # vol2.update_ctf = True
    # vol2._otf = otf
    # vol2._volume_property.set_scalar_opacity(otf)

    mlab.view(azimuth=(-150+(i*5)), elevation=60, distance=500, focalpoint=(49,49,49))
    # mlab.move(right=30)
    mlab.savefig(filename='./anim/fig%03d.png'%i)
    # mlab.show()
    mlab.close()


# mlab.show()
#ffmpeg -r 10 -i fig%03d.png -c:v libx264 -vf fps=60 -pix_fmt yuv420p out.mp4
#ffmpeg -f concat -i list.txt -c copy output.mp4
