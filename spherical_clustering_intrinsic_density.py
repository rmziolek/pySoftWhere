import numpy as np
import matplotlib 
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import seaborn as sns

import MDAnalysis as mda
from MDAnalysis.analysis import pca
from MDAnalysis.analysis.align import rotation_matrix

import matplotlib.colors

import time
from numpy.linalg import norm
from scipy.interpolate import RectBivariateSpline
import scipy.stats as stats
import networkx as nx
from collections import defaultdict
import MDAnalysis.analysis.rms
import MDAnalysis
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
import MDAnalysis.analysis.distances
import MDAnalysis.analysis.rdf


###
topology_ = '../production.tpr' 

trajectory_ = '../production.trr'

start_ = 8000
skip_ =5



u=mda.Universe(topology_,trajectory_)
print(u.trajectory)
####select micelle
micelle=u.select_atoms('not name HW? OW')
micelle.unwrap()
####select hydrophobic core and hydrophilic shell
micelle_hydrophilic=micelle.select_atoms('resname TX4 and not name C1D C1C C1B C1A C1E C1F C2 C3A C3B C4 C5 C6C C6B C6A ')
micelle_hydrophobic=micelle.select_atoms('resname TX4 and name C1D C1C C1B C1A C1E C1F C2 C3A C3B C4 C5 C6C C6B C6A ')
water=u.select_atoms('name OW')

measure_intrinsic_density=water # choose micelle_hydrophobic, micelle_hydrophilic, water


no_bins=31 #this actually gives (N-1)**2 bins (this is the number of edges in 1D)
cluster_atoms=u.select_atoms('name C1A')
cutoff_cluster=10.0
###







def cartesian_to_spherical(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
    ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew


no_polymers=len(cluster_atoms)

intrinsic_profile=np.zeros(1499)
intrinsic_rand=np.zeros(1499)



####


output_file_water3="tx114_cluster_size2.txt"

output_file_water4="tx114_core_rad2.txt"

####

for output_file in [output_file_water3,output_file_water4]:
    try:
        open(output_file,'x')
    except FileExistsError:
        open(output_file,'w').close()


for ts in u.trajectory[start_::skip_]:
    
    ######clustering check
    
    cluster_sizes=[]
    
    #distance_array_com=np.zeros([len(ht304_residues),3],dtype=np.float32)
    box_=u.dimensions

#    for res in ht304_residues:
#        distance_array_com[res-1]=ht304.select_atoms('name N1 and resid '+str(res)).center_of_mass() ##can I take this out of the for loop? does it update anyway?
    dist_array=MDAnalysis.analysis.distances.contact_matrix(cluster_atoms.atoms.positions,cutoff=cutoff_cluster,box=box_)  #does this work properly with pbc?

    pairs=[]
    source_id=0
    for source in dist_array:

        for target in range(len(dist_array)):
            if source[target]==True and source_id!=target:
                pairs.append([source_id,target])

        source_id+=1

    G=nx.Graph()
    for pair in pairs:
        G.add_edge(int(pair[0]),int(pair[1]))


    clusters=[h.nodes for h in nx.connected_component_subgraphs(G)]

    cluster_sizes=[len(h) for h in nx.connected_component_subgraphs(G)]
    
    no_ones=no_polymers-sum(cluster_sizes)

    [cluster_sizes.append(1) for i in range(no_ones)]

    
    for cluster in clusters:
        print(str(' '.join(list([str(i) for i in cluster])))) #resids start from zero
        
        if len([(i) for i in cluster])>10:
    ######

            water_dr_w=[]
            water_r_w = []
            monte_carlo_dr_w=[]
            av_core_rad_w=[]
            


            box_  = u.dimensions
            box   = u.dimensions

            micelle=u.select_atoms('not name HW? OW and prop mass > 2 and resid '+str(' '.join(list([str(i) for i in cluster]))))
            micelle_hydrophilic=micelle.select_atoms('resname TX4 and prop mass > 2 and not name C1D C1C C1B C1A C1E C1F C2 C3A C3B C4 C5 C6C C6B C6A and resid '+str(' '.join(list([str(i) for i in cluster]))))
            micelle_hydrophobic=micelle.select_atoms('resname TX4 and prop mass > 2 and name C1D C1C C1B C1A C1E C1F C2 C3A C3B C4 C5 C6C C6B C6A and resid '+str(' '.join(list([str(i) for i in cluster]))))
            
            micelle.wrap()
            
            com = sum([(1/(sum(micelle.masses)))*micelle.masses[i]*micelle.positions[i] for i in range(len(micelle))])
            print(com)




            test_array=cartesian_to_spherical(micelle_hydrophobic.positions-com)
            r_=test_array[:,3]
            theta_=test_array[:,4]
            phi_=test_array[:,5]




            ###change this for the different densities!
            water_positions=water.positions-com ##water_o_positions,po_new,eo_new,n_new

            for dimension in range(3):
                for atom in range(len(water_positions)):
                    if water_positions[atom][dimension]>box_[dimension]/2:
                        water_positions[atom][dimension]=water_positions[atom][dimension]-box_[dimension]
                    elif water_positions[atom][dimension]<-box_[dimension]/2:
                        water_positions[atom][dimension]=water_positions[atom][dimension]+box_[dimension]

            test_array3 = cartesian_to_spherical(water_positions)
            r_3_tmp=test_array3[:,3]
            theta_3_tmp=test_array3[:,4]
            phi_3_tmp=test_array3[:,5]

            r_3_sphere = [i for i in r_3_tmp]# if i<r_max]
            theta_3__sphere=[theta_3_tmp[i] for i in range(len(r_3_tmp))]# if r_3_tmp[i]<r_max]
            phi__sphere=[phi_3_tmp[i] for i in range(len(r_3_tmp))]# if r_3_tmp[i]<r_max]


            interface_ppo=stats.binned_statistic_2d([np.cos(i) for i in theta_],phi_,r_,bins=[np.linspace(-1,1,no_bins),np.linspace(-np.pi,np.pi,no_bins)],expand_binnumbers=True,statistic=np.max)#,statistic=find_half_density) ###np.max() works fine
            W=stats.binned_statistic_2d([np.cos(i) for i in theta_3__sphere],phi__sphere,r_3_sphere,bins=[np.linspace(-1,1,no_bins),np.linspace(-np.pi,np.pi,no_bins)],expand_binnumbers=True,statistic='count')

            interface_ppo_vals=interface_ppo.statistic.copy()


            while np.count_nonzero(~np.isnan(interface_ppo_vals))!=(len(np.linspace(0,box_[0],num=no_bins))-1)*(len(np.linspace(0,box_[0],num=no_bins))-1):
                for i in range(len(interface_ppo_vals)):
                    for j in range(len(interface_ppo_vals)):

                        if np.isnan(interface_ppo_vals[i][j]):
                            n_i=[i-1,i,i+1]
                            n_j=[j-1,j,j+1]
                            interface_ppo_vals[i][j]=np.nanmean([interface_ppo_vals[divmod(ip,len(interface_ppo_vals))[1]][divmod(jp,len(interface_ppo_vals))[1]] for ip,jp in product(n_i,n_j)])




            for i in range(len(r_3_sphere)):
                bin_X=W.binnumber[0][i]-1
                bin_Y=W.binnumber[1][i]-1


                interface_pos=interface_ppo_vals[bin_X][bin_Y]
                water_dr_w.append(-interface_pos+r_3_sphere[i])
                water_r_w.append(r_3_sphere[i])


            ### monte carlo normalisation
            rand_xyz=np.random.uniform(-u.dimensions[0]/2,u.dimensions[0]/2,size=(10*len(r_3_sphere),3)) #water_o_positions  ###was r_3_sphere


            test_array3 = cartesian_to_spherical(rand_xyz)

            rand_r_tmp=test_array3[:,3]
            rand_cos_theta_tmp=test_array3[:,4]
            rand_phi_tmp=test_array3[:,5]

            rand_r = [i for i in rand_r_tmp][:5*len(r_3_sphere)]# if i<r_max][:len(r_3_sphere)]
            rand_cos_theta=rand_cos_theta_tmp[:5*len(r_3_sphere)]
            rand_phi=rand_phi_tmp[:5*len(r_3_sphere)]



            interface=stats.binned_statistic_2d([np.cos(i) for i in theta_],phi_,r_,bins=[np.linspace(-1,1,no_bins),np.linspace(-np.pi,np.pi,no_bins)],expand_binnumbers=True,statistic=np.max)
            interface_vals=interface.statistic.copy()

            W=stats.binned_statistic_2d([np.cos(i) for i in rand_cos_theta],rand_phi,rand_r,bins=[np.linspace(-1,1,no_bins),np.linspace(-np.pi,np.pi,no_bins)],expand_binnumbers=True,statistic='count')

            while np.count_nonzero(~np.isnan(interface_vals))!=(len(np.linspace(0,box_[0],num=no_bins))-1)*(len(np.linspace(0,box_[0],num=no_bins))-1):
                for i in range(len(interface_vals)):
                    for j in range(len(interface_vals)):

                        if np.isnan(interface_vals[i][j]):
                            n_i=[i-1,i,i+1]
                            n_j=[j-1,j,j+1]
                            interface_vals[i][j]=np.nanmean([interface_vals[divmod(ip,len(interface_vals))[1]][divmod(jp,len(interface_vals))[1]] for ip,jp in product(n_i,n_j)])

            for i in range(0,len(rand_r)):
                bin_X=divmod(W.binnumber[0][i]-1,no_bins-1)[1]
                bin_Y=divmod(W.binnumber[1][i]-1,no_bins-1)[1]
                interface_pos=interface_vals[bin_X][bin_Y]

                monte_carlo_dr_w.append(-interface_pos+rand_r[i])

            av_core_rad_w.append(np.nanmean(interface_vals))
            print(ts)
            
            
            ####save the profile at each timestep
            
            print('calculating intrinsic profile')
            
            intrinsic_density,a,b=plt.hist(water_dr_w,bins=np.arange(-50,100,0.1),density=False)
            intrinsic_density_r,a,b=plt.hist(monte_carlo_dr_w,bins=np.arange(-50,100,0.1),density=False)
            plt.close()
            
            for i in range(len(intrinsic_profile)):

                intrinsic_profile[i]+=intrinsic_density[i]
                intrinsic_rand[i]+=intrinsic_density_r[i]
            
            np.savetxt('profile_tx114_water.npy',intrinsic_profile)
            np.savetxt('profile_tx114_water_rand.npy',intrinsic_rand)


            with open(output_file_water3, 'a') as myfile:
                    myfile.write(str(len(cluster))+"\n")


            with open(output_file_water4, 'a') as myfile:
                    myfile.write(str(np.nanmean(interface_vals))+"\n")



