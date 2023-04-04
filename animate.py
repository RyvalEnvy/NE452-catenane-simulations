
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import argparse

parser = argparse.ArgumentParser(description='Animate some rings!!!')
parser.add_argument('-d1', '--data1', help='first dataset, np.array(time, axes, particles)')
parser.add_argument('-d2', '--data2', help='second dataset, np.array(time, axes, particles)')
parser.add_argument('-mins', '--minimums', help='plotting display minimums, [x,y,x]')
parser.add_argument('-maxs', '--maximums', help='plotting display maximums, [x,y,x]')
parser.add_argument('-c', '--count', help='number of time steps, int')
parser.add_argument('-di', '--dispinit', help='display initial particles, bool, default=True')
parser.add_argument('-ds', '--dispsurf', help='display surface history, bool, default=True')
parser.add_argument('-dc', '--dispcurr', help='display current particles, bool, default=True')

def animate_data(data1, data2, mins, maxs, n_t, initials = True, hists = True, curr = True):

   # this is the function used in the FuncAnimation function
   def animate_func(num):

       # clear current figure
       ax.clear()     

       if hists:
           # update the data for the line, use num+1 bc Python indexing
           ax.plot_surface(data1[:num+1,0,:], data1[:num+1,1,:],data1[:num+1,2,:], color='blue')
           ax.plot_surface(data2[:num+1,0,:], data2[:num+1,1,:],data2[:num+1,2,:], color='orange')

       if curr:
           # add a point for the current particle point
           ax.scatter(data1[num,0,:], data1[num,1,:], data1[num,2,:], c='blue', marker='o')
           ax.scatter(data2[num,0,:], data2[num,1,:], data2[num,2,:], c='orange', marker='o')

       if initials:
           # add a point for the particle start point
           ax.plot3D(data1[0,0,:], data1[0,1,:], data1[0,2,:], c='black', marker='o')
           ax.plot3D(data2[0,0,:], data2[0,1,:], data2[0,2,:], c='red', marker='o') 
       
       # Setting Axes Limits
       ax.set_xlim3d([mins[0], maxs[0]])
       ax.set_ylim3d([mins[1], maxs[1]])
       ax.set_zlim3d([mins[2], maxs[2]])

       # Adding Figure Labels
       ax.set_title('Trajectory \nStep = ' + str(num))
       ax.set_xlabel('x')
       ax.set_ylabel('y')
       ax.set_zlabel('z')
       
   # plotting the Animation
   fig = plt.figure()
   ax = plt.axes(projection='3d')
   line_ani = animation.FuncAnimation(fig, animate_func, interval=100, frames=n_t-1)
   
   print ("done")
   plt.show()
   
   # Saving the Animation
   f = r"animate_func_test.gif"
   writergif = animation.PillowWriter(fps=n_t/6)
   line_ani.save(f, writer=writergif)
   
   return line_ani

if __name__ == "__main__":
   args=parser.parse_args()
   print(args)
#     animate_data(args)
