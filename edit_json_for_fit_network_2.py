#!/usr/bin/env python
# coding: utf-8

import json, os


current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]
        
output_dir_path   = "{}/Fitted-Network".format(current_path)
load_network_path = "{}/Fitted-Network/Network-1/fit_network.dill.0".format(current_path)
input_file_path   = "{}/data/Logs-en.xyz.1".format(current_path)


# 读取JSON文件
with open("{}/Program_sub_script/fit_network_2.json".format(current_path), "r") as file:
     data = json.load(file)

# 修改JSON数据
data["output_dir"] = output_dir_path
data["fitting_net"]["load_network"] = load_network_path
data["fitting"]["input_file"] = input_file_path
#data["hobbies"].append("swimming")

# 将修改后的数据写回文件中
with open("{}/Program_sub_script/fit_network_2.json".format(current_path), "w") as file:
     json.dump(data, file, indent=4)


