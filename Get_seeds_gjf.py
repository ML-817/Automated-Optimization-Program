#!/usr/bin/env python
# coding: utf-8

import os
import time

current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]

input_information_get_cluster = '{}/Edit_info_file/gene-v05.dat'.format(current_path)
with open("{}".format(input_information_get_cluster),"r") as f:
     inf_line = f.readlines()

cluster      = inf_line[24].split()[2]

cluster_name = inf_line[24].split()[2]
if cluster_name[-1] == '-' or cluster_name[-1] == '+':
    cluster_name = cluster_name[0:-1]
else:
    cluster_name = cluster_name

element_1 = inf_line[0].split()[0]
element_2 = inf_line[1].split()[0]
element_3 = inf_line[2].split()[0]
num_element_1 = inf_line[0].split()[2]
num_element_2 = inf_line[0].split()[2]
num_element_3 = inf_line[0].split()[2]

element_list = [element_1, element_2, element_3]
num_element_list = [num_element_1, num_element_2, num_element_3]

element = []
num_each_element=[]

num_element_j = 0
for i in element_list:
    if i != 'XX':
        element.append(i)
        num_each_element.append(num_element_list[num_element_j])
    num_element_j += 1

Min_Multi                = int(inf_line[29].split()[1])
Max_Multi                = int(inf_line[29].split()[1])
#num_initial_structs_list = inf_line[6].split()[2].split(',')
#num_structs_Min_Multi    = int(num_initial_structs_list[0])
#num_structs_other_Multi  = int(num_initial_structs_list[1])
d_length                 = inf_line[17].split()[1]
#server_num               = inf_line[8].split()[2].split(',')
user_name                = inf_line[24].split()[1]
#Guassian_tpye            = inf_line[24].split()[0]
xc_function              = inf_line[28].split('#')[1].split('/')[0]
basis_set                = inf_line[28].split('/')[1].split()[0]
DFT_steps                = int(inf_line[28].split()[1].split('=')[1].split(')')[0])
charge                   = int(inf_line[29].split()[0])
key_word_list            = inf_line[15].split()[3:]
key_word                 = ' '.join(str(i) for i in key_word_list)
memory                   = inf_line[26].split('=')[1].split('M')[0]
nprocshared              = inf_line[27].split('=')[1].split()[0]

############################################## Define mkdir function ###################################################
# 创建安放初始结构文件夹
def mkdir(path):
    # 去除首位空格
    path = path.strip()
    # 去除尾部 \ 符号
    path = path.rstrip("\\")
    # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists = os.path.exists(path)
    # 判断结果
    if not isExists:
       # 如果不存在则创建目录
       # 创建目录操作函数
       os.makedirs(path)
       print(path + ' 创建成功')
       return True
    else:
       # 如果目录存在则不创建，并提示目录已存在
       print(path + ' 目录已存在')
       return False
#################################################################################################

mkpath = "{}/Initial_seeds_gjf".format(current_path)
# 调用函数
mkdir(mkpath)

###################################### 用神经网络产生的初始结构编辑高斯输入文件 ########################################
#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import pandas as pd
import sys, os
sys.path.append(os.path.dirname(os.path.expanduser('~/Tools/')))
from read_xyz import Read_xyz
from write_xyz import write_xyz
from write_gjf_auto import write_gjf
#if len(element_row) == 2:
 #  from write_gjf_MO_auto import write_gjf

pd.options.display.max_rows = None
pd.options.display.max_colwidth = 100
pd.options.display.float_format = '{:.8f}'.format


# ### Load data from structures

# In[2]:
def is_plane(cords):
    return (cords == 0).all(axis=0).any()


for i in range( 0, int(0.5 * (Max_Multi-Min_Multi) + 1) ):
    st = Read_xyz('{}/Initial_seeds_xyz/fil_structs.xyz.0'.format(current_path))
    en, xyz = st.get_xyz()
    enarr, xyzarr = st.df2array()

    mask = np.empty(0, dtype=np.bool_)
    for cor in xyzarr:
        mask = np.append(mask, is_plane(cor))

    no_plane_xyzarr = xyzarr[~mask, :, :]
    no_plane_enarr = enarr[~mask, :]

    write_xyz(fout='{}/Initial_seeds_xyz/fil_structs.xyz'.format(current_path),
              xyz=no_plane_xyzarr,
              en=no_plane_enarr,
              num_atoms=st.num_atoms)

    st = Read_xyz('{}/Initial_seeds_xyz/fil_structs.xyz.1'.format(current_path))
    en, xyz = st.get_xyz()
    enarr, xyzarr = st.df2array()
   

    # ### Write .gjf files

    # In[5]:

    
    write_gjf(output_path      ='{}/Initial_seeds_gjf'.format(current_path),
              format_file      = input_information_get_cluster,
              cluster          = cluster,
              basis_set_info_line  = 32 , # in this case, let's copy line 18-end
              num_atoms        = st.num_atoms,
              num_frames       = st.num_frames,
              xc_function      = xc_function,
              basis_set        = basis_set,
              DFT_steps        = DFT_steps,
              charge           = charge,
              key_word         = key_word,
              memory           = memory,
              nprocshared      = nprocshared,   
              coord            = xyz,
              el               = element, 
              num_each_element = num_each_element,
              M                = (Min_Multi + 2 * i)
              )
########################################################################################################################



