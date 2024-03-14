#!/usr/bin/env python
# coding: utf-8

import json, os


current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]

with open('{}/Edit_info_file/gene-v05.dat'.format(current_path),'r') as f:
    Information = f.readlines()

Pretrain_cluster = Information[21].split()[3]
Multi = Information[21].split()[4]
user_name = Information[24].split()[1]

load_network_path = "/home/{}/GA_DNN_TL_OPT_STRUCTS/{}/{}/Fitted-Network/fit_network.dill.0".format(user_name, Pretrain_cluster, Multi)

output_dir_path = "{}/Fitted-Network".format(current_path)

input_file_path = "{}/data/Logs-en.xyz.0".format(current_path)

# 读取JSON文件
with open("{}/Program_sub_script/Get_Pretrain_DNN_for_fit.json".format(current_path), "r") as file:
     data = json.load(file)

# 修改JSON数据
data["fitting_net"]["load_network"] = load_network_path
data["output_dir"]                  = output_dir_path
data["fitting"]["input_file"]       = input_file_path
#data["hobbies"].append("swimming")

# 将修改后的数据写回文件中
with open("{}/Program_sub_script/Get_Pretrain_DNN_for_fit.json".format(current_path), "w") as file:
     json.dump(data, file, indent=4)


