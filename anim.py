# import libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import argparse

# need this to display plots in jupyter notebook
# %matplotlib notebook

# setup the argument parser
parser = argparse.ArgumentParser(description='Animate some rings!!!')
parser.add_argument('-d1', '--data1', help='first dataset, np.array(time, axes, particles)')
parser.add_argument('-d2', '--data2', help='second dataset, np.array(time, axes, particles)')
parser.add_argument('-mins', '--minimums', help='plotting display minimums, [x,y,x]')
parser.add_argument('-maxs', '--maximums', help='plotting display maximums, [x,y,x]')
parser.add_argument('-c', '--count', help='number of time steps, int')
parser.add_argument('-n', '--name', help='save file name, string')
parser.add_argument('-di', '--dispinit', help='display initial particles, bool, default=True')
parser.add_argument('-ds', '--dispsurf', help='display surface history, bool, default=True')
parser.add_argument('-dc', '--dispcurr', help='display current particles, bool, default=True')
parser.print_help()

def reload(T = 1, gamma = 1, length = 1, dt = 0.1, skipSteps = 1, 
           f_pdb = "pdb", f_ff = "forcefield", srSize = 10):
    fname = "{}fs_{}ps_{}K_{}ss.xyz".format(dt, length, T, skipSteps)
    totSteps = int(length*1000/(skipSteps * dt))
    print("totSteps", totSteps)

    a = np.loadtxt(fname)
    nParticles = a.shape[-1]//3
    a = a.reshape(totSteps, 3, nParticles)
    
    return a

def animate_data(data1, data2, mins, maxs, n_t, name, initials = True, hists = True, curr = True):

   # this is the function used in the FuncAnimation function
   def animate_func(num):

       # clear current figure
       ax.clear()     

       if hists:
           # update the data for the line, use num+1 bc Python indexing
           ax.plot_surface(data1[:num+1,0,:], data1[:num+1,1,:],
                           data1[:num+1,2,:], color='blue')
           ax.plot_surface(data2[:num+1,0,:], data2[:num+1,1,:],
                           data2[:num+1,2,:], color='orange')

       if curr:
           # add a point for the current particle point
           ax.scatter(data1[num,0,:], data1[num,1,:], data1[num,2,:], 
                      c='blue', marker='o')
           ax.scatter(data2[num,0,:], data2[num,1,:], data2[num,2,:], 
                      c='orange', marker='o')

       if initials:
           # add a point for the particle start point
           ax.plot3D(data1[0,0,:], data1[0,1,:], data1[0,2,:], c='black', 
                     marker='o')
           ax.plot3D(data2[0,0,:], data2[0,1,:], data2[0,2,:], c='red', 
                     marker='o') 


#         # update the data for the line, use num+1 bc Python indexing
#         ax.plot_surface(data1[:num+1,0,:], data1[:num+1,1,:],data1[:num+1,2,:], color='blue')
#         ax.plot_surface(data2[:num+1,0,:], data2[:num+1,1,:],data2[:num+1,2,:], color='orange')


#         # add a point for the current particle point
#         ax.scatter(data1[num,0,:], data1[num,1,:], data1[num,2,:], c='blue', marker='o')
#         ax.scatter(data2[num,0,:], data2[num,1,:], data2[num,2,:], c='orange', marker='o')


#         # add a point for the particle start point
#         ax.plot3D(data1[0,0,:], data1[0,1,:], data1[0,2,:], c='black', marker='o')
#         ax.plot3D(data2[0,0,:], data2[0,1,:], data2[0,2,:], c='red', marker='o') 
       
       # Setting Axes Limits
       ax.set_xlim3d([mins[0], maxs[0]])
       ax.set_ylim3d([mins[1], maxs[1]])
       ax.set_zlim3d([mins[2], maxs[2]])

       # Adding Figure Labels
       ax.set_title('Trajectory \nStep = ' + str(num))
       ax.set_xlabel('x')
       ax.set_ylabel('y')
       ax.set_zlabel('z')
       
   # plotting the animation
   fig = plt.figure()
   ax = plt.axes(projection='3d')
   line_ani = animation.FuncAnimation(fig, animate_func, interval=100, frames=n_t-1)
   
   print ("done")
   plt.show()
   
   # Saving the Animation
   writergif = animation.PillowWriter(fps=n_t/6)
   line_ani.save(name, writer=writergif)
   
   return line_ani

if __name__ == "__main__":
    # load data
    a = reload()

    smSize = 10
    simStart = 0
    simEnd = 500
    simLen = simStart - simEnd

    # parse data into separate loops
    smTraj = a[simStart:simEnd, :, 0:smSize]
    lgTraj = a[simStart:simEnd, :, smSize:]

    print(smTraj.shape)
    print(lgTraj.shape)

    # find the max and min of each dim
    mins = [] # x, y, and z lims
    maxs = []
    for i in range(3):
        minimum = np.min(a[simStart:simEnd, i, :])
        maximum = np.max(a[simStart:simEnd, i, :])
        
        mins.append(minimum)
        maxs.append(maximum)

    print("mins:", mins)
    print("maxs:", maxs)

    mins = [-5, -5, -5]
    mins = [5, 5, 5]

    name = "trial1.gif"

    animate_data(data1 = smTraj, data2 = lgTraj, mins = mins, maxs = maxs, 
                 n_t = simLen, name = name, hists = False)