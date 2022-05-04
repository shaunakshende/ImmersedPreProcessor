#!/usr/bin/env python
# coding: utf-8

# In[1]:


import netCDF4 # a module to read the mesh file 
import numpy as np # for array manipulation
from scipy.spatial import ConvexHull # provides a function to calculate volume of polyhedra (using convex hull)


# In[2]:


import os
output_dir = './point_cloud' # the converted mesh will be output here
#Create the directory if doesn't exist
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)


# In[3]:


#Read the mesh file (genesis/exodus)
nc = netCDF4.Dataset('./mesh0_w_void0.g')


# In[4]:


nc.variables.keys()


# In[5]:


#Fetch the FEM mesh info
FEM_node_coordx = nc['coordx']
FEM_node_coordy = nc['coordy']
FEM_node_coordz = nc['coordz']

#Find the number of blocks and elements in the mesh
num_blocks = nc['eb_status'].shape[0]

#Iterate over all the blocks 
for block_id in range(1,num_blocks+1):
    
    #Store the connectivity map for this block
    connect_this_block = nc['connect'+str(block_id)]

    num_elem = connect_this_block.shape[0]

    #Initialize arrays for the point coordinates and volumes
    point_coordinates = np.zeros([num_elem, 3]) # 3d mesh
    point_volume = np.zeros(num_elem)

    #Iterate over each element in the block
    #compute its volume and centroid
    for elem_id, elem_connectivity in enumerate(connect_this_block[:]):
    
        elem_connectivity -= 1 # list starts from 1; python arrays start from 0

        elem_nodes_coordx = FEM_node_coordx[elem_connectivity]
        elem_nodes_coordy = FEM_node_coordy[elem_connectivity]
        elem_nodes_coordz = FEM_node_coordz[elem_connectivity]
        elem_nodes_coordinates = np.array([elem_nodes_coordx, elem_nodes_coordy, elem_nodes_coordz]).T

        #Create a Convex Hull using the coordinates and compute the element volume
        elem_volume = ConvexHull(elem_nodes_coordinates).volume
        point_volume[elem_id] = elem_volume

        #Compute the element centroid to place a meshfree node
        centroid_coordinates = elem_nodes_coordinates.mean(axis=0)
        point_coordinates[elem_id] = centroid_coordinates
        
    ###Write this block to a file
    
    ##Coordinates
    file_name = output_dir + '/input_coor_' + str(block_id) + '.dat'

    #Open the output file
    f=open(file_name,'wb')

    #Write the header
    f.write((str(num_elem)+'\n').encode())

    #Write the point cloud coordinates
    #point_num coord_x coord_y coord_z 
    np.savetxt(f, np.c_[np.arange(1,num_elem+1),point_coordinates], fmt="%d %f %f %f")

    #Close the file
    f.close()

    ##Volume
    file_name = output_dir + '/XVOL_' + str(block_id) + '.dat'

    #Open the output file
    f=open(file_name,'wb')

    #Write the header
    f.write((str(num_elem)+'\n').encode())

    #Write the point cloud volume
    np.savetxt(f, point_volume, fmt="%f")

    #Close the file
    f.close()


# In[ ]:




