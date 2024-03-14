#!/usr/bin/env python
# coding: utf-8

import json, os


current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]

with open('{}/Edit_info_file/gene-v05.dat'.format(current_path), 'r') as f:
    info_lines = f.readlines()

cluster_name = info_lines[24].split()[2]

if cluster_name[-1] == '-' or cluster_name[-1] == '+':
    cluster_name_without_charge_str = cluster_name[0:-1]
else:
    cluster_name_without_charge_str = cluster_name

output_dir_path   = "{}/Initial_seeds_xyz".format(current_path)


# 读取JSON文件
with open("{}/Program_sub_script/BLDA_method_for_generating_initial_seed.json".format(current_path), "r") as file:
     data = json.load(file)

# 修改JSON数据
data["output_dir"] = output_dir_path
data["creation"]["name"] = cluster_name_without_charge_str
#data["hobbies"].append("swimming")

# 将修改后的数据写回文件中
with open("{}/Program_sub_script/BLDA_method_for_generating_initial_seed.json".format(current_path), "w") as file:
     json.dump(data, file, indent=4)


