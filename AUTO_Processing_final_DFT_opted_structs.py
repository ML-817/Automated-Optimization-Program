#!/usr/bin/env python
# coding: utf-8

import os
import time

current_path = os.getcwd()
if current_path[-1] == '/':
    current_path = current_path[0:-1]
    
os.chdir('{}/Edit_info_file/'.format( current_path ) )

input_information_get_cluster = 'DNN_TL_input_for_gjf.txt'
with open("{}".format(input_information_get_cluster),"r") as f:
     inf_line_get_cluster = f.readlines()

cluster_name                  = inf_line_get_cluster[1].split()[2]
os.chdir(current_path )

Max_diff_for_filter = 0.25

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
mkpath = "{}/Processing-{}".format(os.getcwd(), time.strftime("%Y-%m-%d-%H_%M", time.localtime()))
# 调用函数
mkdir(mkpath)

os.chdir('{}'.format( current_path ))

os.system('bash -c "source /home/soft/anaconda3/etc/profile.d/conda.sh && conda activate base && python Read_log_last_struct_xyz_M_en_for_filter.py" ')

Creat_file = '{}/Filter-for-DFT.json'.format(os.getcwd())

with open (Creat_file, 'w') as f:
    f.truncate()
with open (Creat_file, 'a+') as f:
    f.write('{\n')
    f.write('"tasks": ["filter"],\n')
    f.write('"output_dir": "{}",\n'.format(os.getcwd(),mkpath))
    f.write('"random_seed": 0,\n')
    f.write('"filtering": {\n')
    f.write('"sort": true,\n')
    f.write('"input_file": "{}/En_structs_info.txt",'.format(os.getcwd()))
    f.write('"align": true,\n')
    f.write('"pre_sort": true,\n')
    f.write('"max_diff_report": 1.0,"max_diff": {}\n'.format(Max_diff_for_filter))
    f.write('\t}\n')
    f.write('}')

os.system('bash -c "source /home/soft/anaconda3/etc/profile.d/conda.sh && conda activate python2 && acnnmain Filter-for-DFT.json" ')


os.system('mv En_structs_info.txt {}'.format(mkpath))

os.system('bash -c "source /home/soft/anaconda3/etc/profile.d/conda.sh && conda activate base && python Processing_finel_DFT_opted_structs.py" ')


Creat_file_2 = '{}/Draw_cluster.json'.format(os.getcwd())

with open (Creat_file_2, 'w') as f:
    f.truncate()
with open (Creat_file_2, 'a+') as f:
    f.write('{\n')
    f.write('"tasks": ["draw"],\n')
    f.write('"output_dir": "{}",\n'.format(mkpath))
    f.write('"random_seed": 0,\n')
    f.write('"report": {\n')
    f.write('"input_file": "{}/fil_result.xyz.0",\n'.format(os.getcwd()))
    f.write('"number": 9999,\n')
    f.write('"ratio": 1.6\n')
    f.write('\t}\n')
    f.write('}')

os.system('bash -c "source /home/soft/anaconda3/etc/profile.d/conda.sh && conda activate python2 && acnnmain Draw_cluster.json" ')

os.system('mv *.0 {}'.format(mkpath))

os.chdir('{}'.format(mkpath ))

os.rename('report.pdf', 'report_{}.pdf'.format(cluster_name))

os.system('scp -r report_{}.pdf server1:~'.format(cluster_name))

