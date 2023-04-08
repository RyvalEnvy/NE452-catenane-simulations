import numpy as np
import math
from matplotlib import pyplot as plt

def reload(fname,T = 50, gamma = 1/50, length = 1, dt = 0.1, skipSteps = 1, 
           f_pdb = "pdb", f_ff = "forcefield", srSize = 20):
    #fname = "{}fs_{}ps_{}K_{}ss.xyz".format(dt, length, T, skipSteps)
    totSteps = int(length*1000/(skipSteps * dt))
    print("totSteps", totSteps)

    a = np.loadtxt(fname)
    nParticles = a.shape[-1]//3
    a = a.reshape(totSteps, 3, nParticles)
    print(a.shape)

    return a

def euclidean_distance(r1,r2):
    d = np.sqrt(
        np.sum(
            np.square(r1-r2)
        )
    )
    return d

def find_closest_monomer(p_small,p_large):
    n_small = np.shape(p_small)[0]
    n_large = np.shape(p_large)[0]

    distance_matrix = np.zeros((n_small,n_large))

    # Calculate the distance from each unit of the smaller ring to each unit on 
    # the longer ring
    for i in range(n_small):
        for j in range(n_large):
            r1 = p_small[i,:]
            r2 = p_large[j,:]
            distance_matrix[i,j]= euclidean_distance(r1,r2)

    #Return the distance 
    closest_monomer_to_each_bead = np.argmin(distance_matrix,axis=1)
    return closest_monomer_to_each_bead

def calculate_circle_center_of_mass(ring):
    return np.mean(ring,axis=0)

def calculate_vector_magnitude(v):
    return np.sqrt(np.sum(np.square(v)))

def calculate_angular_displacement():
   
    position_time_tensor =reload("animate_bad.xyz")

    small_ring_carbon_indices=(0,20)
    large_ring_carbon_indices=(20,220)
    
    #Extract the small ring positions
    p_small = position_time_tensor[:,:,
                    small_ring_carbon_indices[0]:small_ring_carbon_indices[1]]
        
    #Extract the large ring positions
    p_large = position_time_tensor[:,:,
                    large_ring_carbon_indices[0]:large_ring_carbon_indices[1]]
    
    # Find the closest monomer to each bead for the initial position, 
    # transposing just to give the matrix the right shape 
    initial_closest_monomers_to_each_bead = find_closest_monomer(
        p_small[0,:,:].T,
        p_large[0,:,:].T
    )
    
    initial_index = int(np.mean(initial_closest_monomers_to_each_bead))
    
    nt = np.shape(position_time_tensor)[0]
    angular_displacement = np.zeros(nt)    

    # find the initial radius of the beads
    theta = 2*np.pi/p_large.shape[-1]
    
    for i in range(nt):
        # find the angle within the plane of r_s_CoM and r_l_CoM for each node
        small_center = calculate_circle_center_of_mass(p_small[i,:,:].T)

        # find the current monomer index
        sm_mono_index = find_closest_monomer(np.array([small_center]), 
                                             p_large[i, :, :].T)
        
        # difference in index between init and current 
        delta = sm_mono_index - initial_index

        # angular distance according to the initial monomer geometry
        ang_disp = delta * theta

        angular_displacement[i] = ang_disp
        
    return angular_displacement
    
if __name__ == "__main__":
    angs = calculate_angular_displacement()
    iters = np.arange(angs.shape[0])
    print("shape of angles:", angs.shape)
    print("shape of iters:", iters.shape)
    plt.plot(iters, angs)
    plt.show()
    plt.close()