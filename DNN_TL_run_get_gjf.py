#!/usr/bin/env python
# coding: utf-8

from collections import defaultdict
import os,re
import time

current_path = os.getcwd()
if current_path[-1] == '/':
    current_path = current_path[0:-1]
user_name = current_path.split('/')[2]

os.chdir('{}/Edit_info_file/'.format( os.getcwd() ) )

input_information_get_cluster = 'DNN_TL_info_for_gjf.txt'
with open("{}".format(input_information_get_cluster),"r") as f:
     inf_line_get_cluster = f.readlines()

cluster_0                  = inf_line_get_cluster[1].split()[2]
os.chdir(current_path )
with open ('Run_YQ_PGOPT_AUTO_Record.txt', 'a+') as f_re:
    f_re.write('Run time:{} python DNN_TL_run_get_gjf.py for {} \n'.format( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), cluster_0  ) )
    f_re.write('\n')

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
mkpath = "{}/DNN_TL_{}/structure_generation".format(current_path, cluster_0)
# 调用函数
mkdir(mkpath)

mkpath = "{}/DNN_TL_{}/nn_fitting/Input_information_for_{}_gjf".format(current_path, cluster_0, cluster_0)
# 调用函数
mkdir(mkpath)
############################################ Define copy function ######################################################
# srcfile 需要复制、移动的文件   
# dstpath 目的地址
import glob
import shutil

def mycopyfile(srcfile,dstpath):                       # 复制函数
    if not os.path.isfile(srcfile):
        print ("%s not exist!"%(srcfile))
    else:
        fpath,fname=os.path.split(srcfile)             # 分离文件名和路径
        if not os.path.exists(dstpath):
            os.makedirs(dstpath)                       # 创建路径
        shutil.copy(srcfile, dstpath + fname)          # 复制文件
        print ("copy %s -> %s"%(srcfile, dstpath + fname))
#######################################################################################################################

src_dir_info = '{}/Edit_info_file/'.format(current_path)

dst_dir_info = '{}/DNN_TL_{}/'.format(current_path, cluster_0)

src_file_list =  glob.glob(src_dir_info + '{}'.format(input_information_get_cluster) )

for srcfile in src_file_list:
    mycopyfile(srcfile, dst_dir_info)


input_information = '{}/DNN_TL_{}/DNN_TL_info_for_gjf.txt'.format( current_path, cluster_0  )
with open("{}".format(input_information),"r") as f:
     inf_line = f.readlines()


########################################################################################################

cluster                  = inf_line[1].split()[2]
if cluster[-1] == '-' or cluster[-1] == '+':
    cluster_name = cluster[0:-1]
else:
    cluster_name = cluster

atom_counts = defaultdict(int)
elements = re.findall(r'([A-Z][a-z]*)(\d*)', cluster)
print('elements:',elements)
for element, count in elements:
    if count == '':
        count = 1
    else:
        count = int(count)
    atom_counts[element] += count
atom_counts = dict(atom_counts)
print('atom_counts',atom_counts)
element = atom_counts.keys()
num_each_element = atom_counts.values()
num_atom_cluster = sum(num_each_element)

Min_Multi                = int(inf_line[2].split()[2])
Max_Multi                = int(inf_line[3].split()[2])
num_initial_structs_list = inf_line[4].split()[2].split(',')
num_structs_Min_Multi    = int(num_initial_structs_list[0])
num_structs_other_Multi  = int(num_initial_structs_list[1])
DFT_steps                = int(inf_line[5].split()[-1])
d_length                 = inf_line[6].split()[-1]
server_num               = inf_line[7].split()[2].split(',')
Guassian_tpye            = inf_line[8].split()[-1]
xc_function              = inf_line[9].split()[2]
basis_set                = inf_line[10].split()[2]
charge                   = int(inf_line[11].split()[-1])
key_word_list            = inf_line[12].split()[2:]
key_word                 = ' '.join(str(i) for i in key_word_list)
memory                   = inf_line[13].split()[-1]
nprocshared              = inf_line[14].split()[2]

############################################### Copy input information ################################################

src_dir = '{}/DNN_TL_{}/'.format(current_path, cluster)

dst_dir = '{}/DNN_TL_{}/nn_fitting/Input_information_for_{}_gjf/'.format(current_path, cluster, cluster)

src_file_list  = glob.glob(src_dir + '{}'.format(input_information_get_cluster))

for srcfile in src_file_list:
    mycopyfile(srcfile, dst_dir)

os.rename('{}/DNN_TL_{}/nn_fitting/Input_information_for_{}_gjf/{}'.format(current_path, cluster, cluster, input_information_get_cluster), '{}/DNN_TL_{}/nn_fitting/Input_information_for_{}_gjf/{}_{}'.format(current_path, cluster, cluster,input_information_get_cluster, time.strftime("%Y-%m-%d-%H_%M", time.localtime())))

########################################################################################################################


os.chdir('{}/DNN_TL_{}/structure_generation'.format(current_path, cluster))
for i in range( 0, int(0.5 * (Max_Multi-Min_Multi) + 1) ):
################################################## 产生 creat.json #####################################################
# 定义要创建的目录
    mkpath = "{}/DNN_TL_{}/structure_generation/OUT-{}-gas/M{}".format(current_path, cluster, cluster, (Min_Multi + 2 * i))
# 调用函数
    mkdir(mkpath)

    Creat_file = '{}/DNN_TL_{}/structure_generation/{}-M{}-create.json'.format(current_path, cluster, cluster, (Min_Multi + 2 * i))

    with open (Creat_file, 'w') as f:
        f.truncate()
    with open (Creat_file, 'a+') as f:
        f.write('{\n')
        f.write('"tasks": ["create"],\n')
        f.write('"output_dir": "{}/DNN_TL_{}/structure_generation/OUT-{}-gas/M{}",\n'.format(current_path, cluster, cluster, (Min_Multi + 2 * i)))
        f.write('"random_seed": 0,\n')
        f.write('"creation": {\n')
        print('cluster_name:{}'.format(cluster_name))
        f.write('"name": "{}",\n'.format(cluster_name))
        if (Min_Multi + 2 * i) == Min_Multi :
           f.write('"number": {},\n'.format(num_structs_Min_Multi) )
        else:
           f.write('"number": {},\n'.format(num_structs_other_Multi) )
        f.write('"method": "blda",\n')
        f.write('"order": 2,\n')
        f.write('"2d": {}\n'.format(d_length))
        f.write('\t},\n')
        f.write('"filtering-create": {"max_diff": 0.25,"max_diff_report": 1.00}\n')
        f.write('}')
########################################################################################################################


################################################ 产生初始结构 ##########################################################
    os.system('acnnmain {}-M{}-create.json'.format(cluster, (Min_Multi + 2 * i)))

############################################ 创建生成的 .gjf 放置路径 ##################################################
    mkpath = "{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}".format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i))
# 调用函数
    mkdir(mkpath)
os.chdir('..')
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
    st = Read_xyz('{}/DNN_TL_{}/structure_generation/OUT-{}-gas/M{}/fil_structs.xyz.0'.format(current_path, cluster, cluster, (Min_Multi + 2 * i)))
    en, xyz = st.get_xyz()
    enarr, xyzarr = st.df2array()

    mask = np.empty(0, dtype=np.bool_)
    for cor in xyzarr:
        mask = np.append(mask, is_plane(cor))

    no_plane_xyzarr = xyzarr[~mask, :, :]
    no_plane_enarr = enarr[~mask, :]

    write_xyz(fout='{}/DNN_TL_{}/structure_generation/OUT-{}-gas/M{}/fil_structs.xyz'.format(current_path, cluster, cluster, (Min_Multi + 2 * i)),
              xyz=no_plane_xyzarr,
              en=no_plane_enarr,
              num_atoms=st.num_atoms)

    st = Read_xyz('{}/DNN_TL_{}/structure_generation/OUT-{}-gas/M{}/fil_structs.xyz.1'.format(current_path, cluster, cluster, (Min_Multi + 2 * i)))
    en, xyz = st.get_xyz()
    enarr, xyzarr = st.df2array()
   

    # ### Write .gjf files

    # In[5]:

    
    write_gjf(output_path      ='{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}'.format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i)),
              format_file      = input_information,
              cluster          = cluster,
              basis_set_info_line  = 20 , # in this case, let's copy line 18-end
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



#################################### 创建将.gjf文件传输到其他服务器上的文件夹 ##########################################
    mkpath = "{}/DNN_TL_{}/Gaussian/{}/{}-steps-DFT-outputs/M{}".format(current_path, cluster, cluster, DFT_steps, (Min_Multi + 2 * i))
    mkdir(mkpath)

######################################### 将 .gjf 文件复制到指定文件夹 #################################################
#!/usr/bin/env python
# coding: utf-8

# srcfile 需要复制、移动的文件   
# dstpath 目的地址
import glob
import shutil
import os, sys, stat

def mycopyfile(srcfile,dstpath):                       # 复制函数
    if not os.path.isfile(srcfile):
        print ("%s not exist!"%(srcfile))
    else:
        fpath,fname=os.path.split(srcfile)             # 分离文件名和路径
        if not os.path.exists(dstpath):
            os.makedirs(dstpath)                       # 创建路径
        shutil.copy(srcfile, dstpath + fname)          # 复制文件
        print ("copy %s -> %s"%(srcfile, dstpath + fname))
 
for i in range( 0, int(0.5 * (Max_Multi-Min_Multi) + 1) ): 
   
    src_dir = '{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}/'.format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i))
     
    l  = 0
    if i == 0:
       initial_gjf_block = int(num_structs_Min_Multi) // int(len(server_num))
       run_block         = int(num_structs_Min_Multi) // int((2*len(server_num)))
    else:
       initial_gjf_block = int(num_structs_other_Multi) // int(len(server_num))
       run_block         = int(num_structs_other_Multi) // int((2*len(server_num)))
    for j in range(0,len(server_num)):
       src_file_list = []

       mkpath = "{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}/{}-struct-to-{}-steps-DFT_M{}_server{}".format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i), cluster, DFT_steps, (Min_Multi + 2 * i), server_num[j])
# 调用函数
       mkdir(mkpath)       

       dst_dir = '{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}/{}-struct-to-{}-steps-DFT_M{}_server{}/' .format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i), cluster, DFT_steps, (Min_Multi + 2 * i), server_num[j])   
       
       for k in range( ((j*initial_gjf_block)+1) , (((j+1)*initial_gjf_block)+1) ):
           src_file_list =  src_file_list + glob.glob(src_dir + '{}-M{}-{}-{}-{}.gjf'.format(cluster, (Min_Multi + 2 * i), xc_function, basis_set, k))
       for srcfile in src_file_list:
           mycopyfile(srcfile, dst_dir)                # 复制文件

       Gussian_run_file = '{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}/{}-struct-to-{}-steps-DFT_M{}_server{}/'.format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i), cluster, DFT_steps, (Min_Multi + 2 * i), server_num[j])
       for m in range(1,3):
          f = open(Gussian_run_file + 'run-{}_{}'.format((j+1),m), 'w')
          with open(Gussian_run_file + 'run-{}_{}'.format((j+1),m), 'w') as f:
              f.truncate()
          with open (Gussian_run_file + 'run-{}_{}'.format((j+1),m), 'a+') as f:
              for nn in range( ((l*run_block)+1) , (((l+1)*run_block)+1) ):
                 f.write('g{} {}-M{}-{}-{}-{}.gjf\n'.format(Guassian_tpye, cluster, (Min_Multi + 2 * i), xc_function, basis_set, nn ))
          l = l + 1

       os.chdir('{}/DNN_TL_{}/Gaussian/{}/{}-struct-to-{}-steps-DFT/M{}'.format(current_path, cluster, cluster, cluster, DFT_steps, (Min_Multi + 2 * i)))
       os.system('scp -r {}-struct-to-{}-steps-DFT_M{}_server{} {}@server{}:/home/{}'.format(cluster, DFT_steps, (Min_Multi + 2 * i), server_num[j], user_name, server_num[j], user_name ))
    

#with open ('Run_YQ_PGOPT_AUTO_Record.txt', 'a+') as f_re:
#    f_re.write('*Finished time:{}* \n'.format( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ) )

