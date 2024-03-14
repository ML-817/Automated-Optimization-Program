#!/usr/bin/env python
# coding: utf-8

import json, os


current_path = os.getcwd().split('Program_sub_script')[0]

if current_path[-1] == '/':
    current_path = current_path[0:-1]

output_dir_path   = "{}/".format(current_path)
load_network_path = "{}/Fitted-Network/fit_network.dill.0".format(current_path)

# 读取JSON文件
with open("{}/Program_sub_script/DNN-opt.json".format(current_path), "r") as file:
     data = json.load(file)

input_file_name = data["optimization"]["input_file"].split('/')[-1]

input_file_path   = "{}/{}".format(current_path, input_file_name)

# 修改JSON数据
data["fitting_net"]["load_network"] = load_network_path
data["optimization"]["load_network"] = load_network_path
data["optimization"]["input_file"] = input_file_path
data["output_dir"] = output_dir_path
#data["hobbies"].append("swimming")

# 将修改后的数据写回文件中
with open("{}/Program_sub_script/DNN-opt.json".format(current_path), "w") as file:
     json.dump(data, file, indent=4)


