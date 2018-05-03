import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib import animation
from matplotlib.animation import FuncAnimation

timestep = 500
n= 9 #number of steps
w, h = 3,3 #width, height
y0, x0 = np.mgrid[:h, :w] #intialising mean positions
z = np.ones(x0.shape) #setting color to white at each face

x,y = np.ndarray(shape=(w,h)), np.ndarray(shape=(w,h))

xdelta=np.loadtxt("xdata.dat")
ydelta=np.loadtxt("ydata.dat")

fig=plt.figure()

def update(frame):
	x=x0+xdelta[3*frame:3*frame+3]
	y=y0+ydelta[3*frame:3*frame+3]
	plt.pcolormesh(x, y, z, edgecolor=(0,0,0,1), cmap='Greys', hold=False)
	plt.title("Time="+str(frame))
	plt.xlim(-0.5,w-0.5)
	plt.ylim(-0.5,h-0.5)

anim=FuncAnimation(fig,update, frames=n,interval=timestep, repeat=False)

plt.show()
