import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib.animation import FuncAnimation

timestep = 1000 		#time period between iterations
n= 50 				#number of steps
w, h = 20,20 			#number of nodes in x direction, number of nodes in y direction
y0, x0 = np.mgrid[:h, :w]	#intialising mean positions
z = np.ones(x0.shape) 		#setting color to white at each face

x,y = np.ndarray(shape=(w,h)), np.ndarray(shape=(w,h))

xdelta=np.loadtxt("u_x.txt")
ydelta=np.loadtxt("u_y.txt")

fig=plt.figure()
plt.axes(aspect='equal')

def update(frame):
	x=x0+xdelta[w*frame:w*frame+w]
	y=y0+ydelta[h*frame:h*frame+h]
	plt.pcolormesh(x, y, z, edgecolor=(0,0,0,1), cmap='Greys', hold=False)
	plt.title("Time="+str(frame)+"s")
	plt.axis([-0.5,w-0.5,-0.5,h-0.5])
	
anim=FuncAnimation(fig,update, frames=n,interval=timestep, repeat=False)

plt.show()
