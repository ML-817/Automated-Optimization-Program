#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import sys, os
sys.path.append(os.path.dirname(os.path.expanduser('~/Tools/')))
from read_xyz import Read_xyz
from write_gjf import write_gjf

pd.options.display.max_rows = None
pd.options.display.max_colwidth = 100
pd.options.display.float_format = '{:.8f}'.format


# ### Load data from structures

# In[2]:




# ### Write .gjf files
current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]

input_information = '{}/Edit_info_file/gene-v05.dat'.format(current_path)
with open("{}".format(input_information),"r") as f:
         Info_lines = f.readlines()

each_ele_type               = [Info_lines[0].split()[0], Info_lines[1].split()[0], Info_lines[2].split()[0] ]

each_ele_type_num           = [int(Info_lines[0].split()[2]), int(Info_lines[1].split()[2]), int(Info_lines[2].split()[2])]    #团簇中个元素数目如：单一元素：[5]；两种元素：[2,7] ;三种元素：[2,3,4]

each_ele_type_num_without_0 = list(filter(lambda  x: x!=0, each_ele_type_num))

num_atoms                   = sum(each_ele_type_num_without_0)   #团簇中的原子数目

ele_atomic_num              = [int(Info_lines[0].split()[1]), int(Info_lines[1].split()[1]), int(Info_lines[2].split()[1])]

ele_atomic_num_dict         = {each_ele_type[0]:Info_lines[0].split()[1], each_ele_type[1]:Info_lines[1].split()[1], each_ele_type[2]:Info_lines[2].split()[1]}


cluster_name  = Info_lines[24].split()[2]
M             = Info_lines[24].split()[3]
functional    = Info_lines[28].split('#')[1].split('/')[0]
basis_set     = Info_lines[28].split('/')[1].split()[0]

Initial_seeds_gjf_path = '{}/Initial_seeds_gjf/'.format(current_path)
num_logs      = int(len(os.listdir(Initial_seeds_gjf_path)))



num_before_coord = 8
Atomic_Type   = 0

# In[5]:
input_path='{}/Initial_seeds_gjf'.format(current_path)
output_path='{}/'.format(current_path)


for i in range(1, num_logs + 1):
    input_file = os.path.join(input_path, '{}-{}-{}-{}-{}.gjf'.format(cluster_name, M, functional, basis_set, i))
    output_file = os.path.join(output_path, '{}-{}-{}.log'.format(cluster_name, M, i))
    with open(output_file, 'w') as fwrite:
        fwrite.write(' ---------------------------------------------------------------------\n')
        fwrite.write(' #hf opt(maxcycle=999)\n')
        fwrite.write(' ---------------------------------------------------------------------\n')
        fwrite.write('\n')
        fwrite.write('                         Standard orientation:\n')
        fwrite.write(' ---------------------------------------------------------------------\n')
        fwrite.write(' Center     Atomic      Atomic             Coordinates (Angstroms)\n')
        fwrite.write(' Number     Number       Type             X           Y           Z\n')
        fwrite.write(' ---------------------------------------------------------------------\n')
        xyz_x=''
        for j in range(1,num_atoms+1):
            with open(input_file, 'r') as fread:
            
                for xyz_line in fread.readlines()[num_before_coord + j - 1:num_before_coord + j]:
                    element = xyz_line.split()[0]
                    xyz_x=xyz_line.split()[1]
                    xyz_y=xyz_line.split()[2]
                    xyz_z=xyz_line.split()[3]
        
            fwrite.write('     {}       {}           {}     {}  {}  {}\n'.format(j,ele_atomic_num_dict[element],Atomic_Type,xyz_x,xyz_y,xyz_z))
        fwrite.write(' ---------------------------------------------------------------------\n')
        fwrite.write(' SCF Done:  E({}) =     0.00000000            A.U. after    5 cycles'.format(functional))






