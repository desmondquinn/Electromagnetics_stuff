""""1D FDTD simulation
with Gaussian pulse (soft) source 
perfectly absorbing boundary condition
Ey/Hx mode"""

import numpy as np
from math import exp
from matplotlib import pyplot as plt
from matplotlib import animation


#no of cells in grid
Nz = 200
#initialising Ex and Hy fields
Hx = np.zeros(Nz)
Ey = np.zeros(Nz)
# Input pulse parameters
kc=5        #point at lower end of the grid where source is injected
t0 = 40
spread = 12

#no of time steps
steps = 850

#Properties of Material 1
mu=np.ones((Nz))
er=np.ones((Nz))

#material 2
start = 100
end = 130
er[start:end] = 3.45
#to draw material boundaries
l = np.linspace(-1,2,1000)
x_start = start*np.ones((1000))
x_end = end*np.ones((1000))

#update coefficients
mEy = 0.5/er
mHx = 0.5/mu

boundary_lo = [0,0]
boundary_hi = [0,0]

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, Nz), ylim=(-1, 2))
point, = ax.plot([], [], lw=2)

time_step = 0

# initialization function: plot the background of each frame
def init():
    global l, x_start, x_end
    plt.plot(x_start,l,'r')
    plt.plot(x_end,l,'r')
    point.set_data([], [])
    return point,


def animate(i):
    global kc, t0, spread, steps, Hx, Ey, mHx, mEy, boundary_lo, boundary_hi,time_step,l, x_start, x_end
    #main FDTD loop
    #for i in np.arange(i):
    time_step=time_step+1
    # Calculate the Hx field
    for k in range(0, Nz-1):
        Hx[k] = Hx[k] + mHx[k] * (Ey[k+1] - Ey[k])
    Hx[Nz-1] = Hx[Nz-1] + mHx[k] * (0-Ey[Nz-1])
    # Gaussian pulse
    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)
    Hx[kc] = pulse + Hx[kc]          #over riding the E value at kc
    # Absorbing Boundary Conditions
    Hx[0] = boundary_lo.pop(0)
    boundary_lo.append(Hx[1])
    Hx[Nz - 1] = boundary_hi.pop(0)
    boundary_hi.append(Hx[Nz - 2])
    # Calculate the Ey field
    Ey[0] = Ey[0] + mEy[k] * (Hx[0]-0)
    for k in range(1,Nz):
        Ey[k] = Ey[k] + mEy[k] * (Hx[k] - Hx[k - 1])
    
    point.set_data(np.arange(Nz),Ey)
    return point,


anim = animation.FuncAnimation(fig, animate, init_func = init, frames = steps, interval = 10, blit = True)
anim.save('1d_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()
