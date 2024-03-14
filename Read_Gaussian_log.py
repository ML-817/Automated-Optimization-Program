#!/usr/bin/env python3
# coding: utf-8

import os,sys
import time
os.system('bash -c "source /home/soft/anaconda3/etc/profile.d/conda.sh && conda activate base" ')

current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]
        
input_information_get_cluster = '{}/Edit_info_file/gene-v05.dat'.format(current_path)
with open("{}".format(input_information_get_cluster),"r") as f:
     inf_line_get_cluster = f.readlines()

element_1 = inf_line_get_cluster[0].split()[0]
element_2 = inf_line_get_cluster[1].split()[0]
element_3 = inf_line_get_cluster[2].split()[0]

element_1_num = int(inf_line_get_cluster[0].split()[2])
element_2_num = int(inf_line_get_cluster[1].split()[2])
element_3_num = int(inf_line_get_cluster[2].split()[2])
num_atom_cluster = element_1_num + element_2_num + element_3_num

DFT_steps = int(inf_line_get_cluster[28].split()[1].split(')')[0].split('=')[1])

#判断元素种类数目
if element_1 != 'XX' and element_2 != 'XX' and element_3 != 'XX':
   element_num = 3
   element     = [element_1,element_2,element_3]
   each_element_num = [element_1_num, element_2_num, element_3_num]

elif (element_1 == 'XX' and element_2 != 'XX' and element_3 != 'XX'):
   element_num = 2
   element     = [element_2,element_3]
   each_element_num = [element_2_num, element_3_num]

elif (element_1 != 'XX' and element_2 == 'XX' and element_3 != 'XX'): 
   element_num = 2
   element     = [element_1,element_3]  
   each_element_num = [element_1_num, element_3_num]

elif (element_1 != 'XX' and element_2 != 'XX' and element_3 == 'XX'):
   element_num = 2
   element     = [element_1,element_2]
   each_element_num = [element_1_num, element_2_num]

elif (element_1 == 'XX' and element_2 == 'XX' and element_3 != 'XX'):
   element_num = 1
   element     = [element_3]
   each_element_num = [element_3_num]

elif (element_1 == 'XX' and element_2 != 'XX' and element_3 == 'XX'): 
   element_num = 1
   element     = [element_2]
   each_element_num = [element_2_num]

elif (element_1 != 'XX' and element_2 == 'XX' and element_3 == 'XX'):
   element_num = 1
   element     = [element_1]
   each_element_num = [element_1_num]

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

data_path = '{}/data'.format(current_path)
if not os.path.exists(data_path):
    os.mkdir( data_path  )

log_path = '{}/All_logs'.format(current_path)
if not os.path.exists(log_path):
    os.mkdir( log_path  )


os.system('cp -f *.log  All_logs')

#src_dir_info = os.getcwd()

#dst_dir_info = log_path

#src_file_list = glob.glob(src_dir_info + '*.log')

#for srcfile in src_file_list:
#    mycopyfile(srcfile, dst_dir_info)


if inf_line_get_cluster[18].find('Iop(5/13=1)') or inf_line_get_cluster[18].find('iop(5/13=1)') or inf_line_get_cluster[18].find('IOP(5/13=1)'):
   iop = True
else:
   iop = False
 
    

#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import pandas as pd
import shutil
#from sklearn import preprocessing
import sys, os, re, glob
sys.path.append(os.path.dirname(os.path.expanduser('~/Tools/')))
from avoid_overwritting import new_file_name

if element_num == 1:
   from read_log_Gaussian import Read_log_Gaussian

if element_num == 2:
   from read_log_Gaussian_MO import Read_log_Gaussian

if element_num == 3:
   from read_log_Gaussian_MO import Read_log_Gaussian

from remove_abnormal_energy import Remove_abnormal_energy

pd.options.display.max_rows = None
pd.options.display.max_colwidth = 100


# ### Write data file from log files

# In[2]:


###### Where is the Gaussian output directory.
dlogs = os.path.join(
        os.path.expandvars('$ACNNHOME'),log_path)
dlogs = os.path.abspath(dlogs)
dlogs = dlogs + '/*.log'
logs = glob.iglob(dlogs)

###### Give a file name to write the structures.
fout = os.path.join(
       os.path.expandvars('$ACNNHOME'),'{}/Logs-en.xyz'.format(data_path))
fout = new_file_name(fout)

for i, log in enumerate(logs):
    print('\n')
    print('The {}th log:'.format(i))
    print('The log path:{}'.format(log))
    popen_input_energy = 'grep "SCF Done" {} | wc -l'.format(log)
    ndone = int(os.popen(popen_input_energy).read())
    if ndone < 3:
        bad_log = '{}'.format(log)
        os.remove(bad_log)
        continue
    t = Read_log_Gaussian(log, element, num_atoms=num_atom_cluster, steps=DFT_steps, IOp= iop )
    if t.type in t.badtypelist:
       continue

    r = Remove_abnormal_energy(t.en, t.xyz)
    en = r.en
    xyz = r.xyz

    if element_num == 1: 
       with open(fout, 'a+') as f:
           for j, (cords, e) in enumerate(zip(xyz, en)):
               f.write('{}\n'.format(t.num_atoms))
               f.write('STR:{}:{} = {:15.8f}\n'.format(i, j, e))
               for cord in cords:
                   f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                           format(t.element, cord[0], cord[1], cord[2]))

    if element_num == 2:
       num_first_atoms = each_element_num[0]

       with open(fout, 'a+') as f:
           for j, (cords, e) in enumerate(zip(xyz, en)):
               f.write('{}\n'.format(t.num_atoms))
               f.write('STR:{}:{} = {:15.8f}\n'.format(i, j, e))
               m = 1

               for cord in cords:
                   if m <= num_first_atoms:
                       f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                               format(t.element[0], cord[0], cord[1], cord[2]))
                       m = m + 1
                   elif m > num_first_atoms:
                       f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                               format(t.element[1], cord[0], cord[1], cord[2]))

    if element_num == 3:
       num_first_atoms = each_element_num[0]
       num_second_atoms = (each_element_num[0] + each_element_num[1])

       with open(fout, 'a+') as f:
           for j, (cords, e) in enumerate(zip(xyz, en)):
               f.write('{}\n'.format(t.num_atoms))
               f.write('STR:{}:{} = {:15.8f}\n'.format(i, j, e))
               m = 1

               for cord in cords:
                  if m <= num_first_atoms:
                      f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'. 
                              format(t.element[0], cord[0], cord[1], cord[2]))
                      m = m + 1   
                  elif num_second_atoms >=  m > num_first_atoms:        
                      f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                              format(t.element[1], cord[0], cord[1], cord[2]))
                      m = m + 1
                  elif m > num_second_atoms:
                      f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                              format(t.element[2], cord[0], cord[1], cord[2]))

