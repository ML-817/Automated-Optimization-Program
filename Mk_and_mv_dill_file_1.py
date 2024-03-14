#!/usr/bin python
# coding: utf-8

import os, shutil, glob

current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]
        
network_1_path = '{}/Fitted-Network/Pretrain-DNN/'.format(current_path)
if not os.path.exists(network_1_path):
   os.mkdir(network_1_path)

src_dir_info  =  '{}/Fitted-Network/'.format(current_path)
src_file_list =  glob.glob(src_dir_info + '*.0')
dst_dir_info  =  network_1_path
for srcfile in src_file_list:
    shutil.move(srcfile,dst_dir_info)

