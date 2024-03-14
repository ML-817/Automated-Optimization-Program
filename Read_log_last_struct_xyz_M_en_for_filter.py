#!/usr/bin/env python
# coding: utf-8


import os
import numpy as np
import sys
import linecache
from collections import defaultdict



def find_n_sub_str(s, sub, times, start=0):
    index = s.find(sub, start)
    if index != -1 and times > 0:
        return find_n_sub_str(s, sub, times - 1, index + 1)
    return index


def get_all_file(base=None):
    result_paths = []
    if base is None:
        base = os.getcwd()
    for file in os.listdir(base):
        path = os.path.join(base, file)
        if os.path.isdir(path):
            result_paths.extend(get_all_file(path))
        else:
            if ".log" in file:
                result_paths.append(path)
    result_paths.sort(key=lambda x: os.path.split(x)[1])
    return result_paths


def file_processor(path):
    block='Input orientation:'
    block_2='NAtoms'
    block_3=' Mulliken charges and spin densities:'
    block_4='Sum of Mulliken charges'
    block_5=' Mulliken atomic charges:'
    block_6=' Sum of Mulliken atomic charges'
    block_7='Multiplicity'
    block_8=' Sum of electronic and zero-point Energies='
    block_9='Mulliken charges:'
    sepst='--'
    file = open(path, "r")
    txt=file.read().splitlines()
    file.close()
    cluster_Coordinates = []
    energy = []
    
    line_num=0
    atom_type_head = 0
    atom_type_1_num = 0
    atom_type_2_num = 0
    atom_type_1 = '*'
    atom_type_2 = '**'

    for line in txt :
        line_num+=1
        if line.find(block_2)>-1:
            nat = int(line.split()[1])   #原子数目

        if line.find(block_7)>-1:
            multi = line.split('=')[2]   #多重度
            charge = line.split()[2]      #电荷数

        if line.find(block_3)>-1 or line.find(block_5)>-1 or line.find(block_9)>-1:
#            print("line_num:",line_num)
            atom_type_head = line_num + 2
            atom_type_tail = atom_type_head + nat
            atom_type_1 = linecache.getline(path, atom_type_head).split()[1]  #第一种元素符号
#            print('atom_type_1:',atom_type_1)
            for line_line in range(atom_type_head,atom_type_tail):
                atom_type = linecache.getline(path, line_line).split()[1]
                if atom_type == atom_type_1:
                    atom_type_1_num += 1    #第一种元素原子数目
                else:
                    atom_type_2_num += 1    #第二种元素原子数目

        if line.find(block_4) > -1 or line.find(block_6)>-1:
#            print('nat:',nat)
#            print('atom_type_1_num:',atom_type_1_num)
            if nat == atom_type_1_num:
                atom_type_2 = None
            else:
                atom_type_2_line = atom_type_head + atom_type_1_num
#                print(atom_type_2_line)
                atom_type_2 = linecache.getline(path, atom_type_2_line).split()[1]  #第二种元素符号
            break

    print(path)
    i = 0
    for line in txt :
        i = i + 1
        if line.find(block)>-1:
            cluster_Coordinates.clear()
            Coordinate_line_head = i + 5
            Coordinate_line_tail = Coordinate_line_head + nat
            for line in range(Coordinate_line_head,Coordinate_line_tail):
                data1 = linecache.getline(path, line)[37:46]    #坐标x
                data2 = linecache.getline(path, line)[49:58]    #坐标y
                data3 = linecache.getline(path, line)[61:70]    #坐标z
                cluster_Coordinates = cluster_Coordinates + [float(data1)] + [float(data2)] + [float(data3)]

    xyz = np.array(cluster_Coordinates).reshape(-1, nat, 3)

    #提取零点矫正能
    for line in txt :
#        print(line)
        if line.find(block_8) > -1:
            energy.clear()
            energy = energy + [float(line.split('=')[1])]

    cluster_Coordinates.clear()

    #得到团簇名称
    if nat == atom_type_1_num and atom_type_1_num != 1:
        cluster_name = str(atom_type_1) + str(atom_type_1_num)
    elif atom_type_1_num == 1 and atom_type_2_num == 1:
        cluster_name = str(atom_type_1) + str(atom_type_2)
    elif atom_type_1_num == 1 and atom_type_2_num != None and atom_type_2_num !=1:
        cluster_name = str(atom_type_1) + str(atom_type_2) + str(atom_type_2_num)
    elif atom_type_1_num != 1 and atom_type_2_num ==1:
        cluster_name = str(atom_type_1) + str(atom_type_1_num) + str(atom_type_2)
    elif atom_type_1_num !=1 and atom_type_2_num != None and atom_type_2_num !=1:
        cluster_name = str(atom_type_1)+str(atom_type_1_num)+str(atom_type_2)+str(atom_type_2_num)


#    print(cluster_name)
#    print(energy)
    return 	xyz,nat,atom_type_1_num,atom_type_2_num,atom_type_1,atom_type_2,multi,energy,charge,cluster_name


def generate_result():
    result_path = get_all_file()
#    csv_output = open("fil_structs.xyz.0", "a+", newline="", encoding="utf-8")

    dic_energy = {}       #定义存储团簇所有能量的字典
    log_path = {}       #定义存储路径的字典
    dic_energy_min = {}   #定义存储团簇最低能量的字典
    cluster_name = []     #定义存储团簇名称的列表

    #得到团簇名称列表和不同团簇的空字典
    for path in result_path:
        result = file_processor(path)
        #记录团簇名称
        cluster_name = cluster_name + [result[9]]
        #定义不同键的空字典
        dic_energy['{}'.format(result[9])] = []
        log_path['{}'.format(result[9])] = []

    for path in result_path:
        result = file_processor(path)
        if len(result[7])  == 0:
            print("The {} have no ' Sum of electronic and zero-point Energies='".format(path))
            sys.exit(0)
        # 相同团簇不同结构的能量归类
        dic_energy['{}'.format(result[9])].append(float(result[7][0]))
        log_path['{}'.format(result[9])].append(path)

    #去掉cluster_name中的重复元素
    all_cluster_name = list(set(cluster_name))
#    print(all_cluster_name)

    #构建存储团簇的最低能量的字典
    dic_energy_min = {}  # 定义存储团簇最低能量的字典
    for cluster_name_v in all_cluster_name:
        energy_min = min(dic_energy[cluster_name_v])
        dic_energy_min['{}'.format(cluster_name_v)] = []
        dic_energy_min['{}'.format(cluster_name_v)].append(energy_min)

    # 对能量排序
    dic_energy_sorted = {}  # 定义存储同一团簇从小到大排序好的能量字典
    for cluster_name_v in all_cluster_name:
        dic_energy_np = np.array(dic_energy[cluster_name_v])
        dic_energy_index = np.argsort(dic_energy_np)
        dic_energy_sorted['{}'.format(cluster_name_v)] = []
        for index_i in dic_energy_index:
            dic_energy_sorted['{}'.format(cluster_name_v)].append(log_path[cluster_name_v][index_i])

#    print(dic_energy_sorted)

    #得到同一团簇按能量从小到大排列的log文件读取新顺序
    new_read_log_list = []
    for cluster_name_v in all_cluster_name:
        new_read_log_list.append(dic_energy_sorted['{}'.format(cluster_name_v)])
#    print('new_read_log_list:',new_read_log_list[0])

    with open("En_structs_info.txt", 'w') as f:
        f.truncate()
    for new_read_log_list_index in range(0,len(all_cluster_name)):        #不同团簇有其对应的log列表，因此根据团簇种类读取其对应的列表
        for path in new_read_log_list[new_read_log_list_index]:
        #        print(path)
            result = file_processor(path)
            Relative_energy = (float(result[7][0]) - float(dic_energy_min['{}'.format(result[9])][0])) * 27.2114
            #输出文件
            with open("En_structs_info.txt", 'a+') as f:
                #f.write('{} (charge = {})\n'.format(result[9],result[8]))
                f.write('{} \n'.format(result[1]))
                if  path.find('\\')>-1:
                    f.write('OPT:{} = {:.3f} eV\n'.format(path.split('\\')[-1] ,Relative_energy))
                if path.find('/') > -1:
                    f.write('OPT:{} = {:.3f} eV\n'.format(path.split('/')[-1], Relative_energy))
                # f.write('M={}\n'.format(result[6]))
                # if path.find('\\')>-1:
                #     f.write('STR:{}\n'.format(path.split('\\')[-1]))
                # if path.find('/')>-1:
                #     f.write('STR:{}\n'.format(path.split('/')[-1]))

                for j,cords in enumerate(result[0]):
                    m = 1
                    cords = cords.tolist()

                    for cord in cords:
                        if m <= result[2]:
                            f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                                       format('{}'.format(result[4]), cord[0], cord[1], cord[2]))
                            m = m + 1
                        elif m > result[2]:
                            f.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.
                                    format('{}'.format(result[5]), cord[0], cord[1], cord[2]))

if __name__ == "__main__":
    generate_result()
