#!/usr/bin/env python
# coding: utf-8
import numpy as np
import os,linecache,csv,math,re

current_path = os.getcwd().split('Program_sub_script')[0]
if current_path[-1] == '/':
    current_path = current_path[0:-1]
        
with open('{}/Edit_info_file/gene-v05.dat'.format(current_path), 'r') as input_info:
    Info_lines = input_info.readlines()

cluster_name                = Info_lines[24].split()[2]    #团簇名称带正负号
each_ele_type               = [Info_lines[0].split()[0], Info_lines[1].split()[0], Info_lines[2].split()[0] ]
each_ele_type_num           = [int(Info_lines[0].split()[2]), int(Info_lines[1].split()[2]), int(Info_lines[2].split()[2])]    #团簇中个元素数目如：单一元素：[5]；两种元素：[2,7] ;三种元素：[2,3,4]
each_ele_type_num_without_0 = list(filter(lambda  x: x!=0, each_ele_type_num))
num_atoms                   = sum(each_ele_type_num_without_0)   #团簇中的原子数目
charge                      = Info_lines[29].split()[0]
multi                       = Info_lines[29].split()[1]    #多重度
num_pick_structs            = int(Info_lines[22].split()[1])    #设定选取的结构数目
out_file_name               = '{}{}{}{}{}{}__ALL.out'.format(each_ele_type[0].ljust(2,'0'), str(each_ele_type_num[0]).zfill(2), each_ele_type[1].ljust(2,'0'), str(each_ele_type_num[1]).zfill(2), each_ele_type[2][0], str(each_ele_type_num[2]).zfill(2))    #out文件名
ele_atomic_num              = [int(Info_lines[0].split()[1]), int(Info_lines[1].split()[1]), int(Info_lines[2].split()[1])]     
ele_atomic_num_dict         = {Info_lines[0].split()[1]:each_ele_type[0], Info_lines[1].split()[1]:each_ele_type[1], Info_lines[2].split()[1]:each_ele_type[2]}


#分配计算任务
if num_pick_structs % 3 == 1:
    Final_GAUSSIAN_run_1_max_serial_num = num_pick_structs/3 + 1
    Final_GAUSSIAN_run_2_max_serial_num = Final_GAUSSIAN_run_1_max_serial_num + num_pick_structs/3
    Final_GAUSSIAN_run_3_max_serial_num = Final_GAUSSIAN_run_2_max_serial_num + num_pick_structs/3
elif num_pick_structs % 3 == 2:
    Final_GAUSSIAN_run_1_max_serial_num = num_pick_structs/3 + 1
    Final_GAUSSIAN_run_2_max_serial_num = Final_GAUSSIAN_run_1_max_serial_num + num_pick_structs/3 + 1
    Final_GAUSSIAN_run_3_max_serial_num = Final_GAUSSIAN_run_2_max_serial_num + num_pick_structs/3


Final_GAUSSIAN_run_1 = '{}/Final_GAUSSIAN_run_1'.format(current_path)
Final_GAUSSIAN_run_2 = '{}/Final_GAUSSIAN_run_2'.format(current_path)
Final_GAUSSIAN_run_3 = '{}/Final_GAUSSIAN_run_3'.format(current_path)

with open(Final_GAUSSIAN_run_1, 'w') as f_GAUSSIAN_run_1:
    f_GAUSSIAN_run_1.truncate()

with open(Final_GAUSSIAN_run_2, 'w') as f_GAUSSIAN_run_2:
    f_GAUSSIAN_run_2.truncate()

with open(Final_GAUSSIAN_run_3, 'w') as f_GAUSSIAN_run_3:
    f_GAUSSIAN_run_3.truncate()

#读取输入文件Read_GA_out_xyz_for_gjf.txt
with open('{}/Edit_info_file/Read_GA_out_xyz_for_gjf.txt'.format(current_path), 'r') as input_info_2:
    Info_lines_2 = input_info_2.readlines()
#GAUSSIAN输入文件信息
memory         = Info_lines[26].split('=')[1].split('M')[0]
nprocshared    = Info_lines[27].split('=')[1].split()[0]
Gaussian_type  = Info_lines_2[2].split()[2]
functional     = Info_lines_2[3].split()[2]
basis_set      = Info_lines_2[4].split()[2]
opt_maxcycle   = Info_lines_2[5].split()[2]
key_words_list = Info_lines_2[6].split()[2:]
print(key_words_list)
key_words      = ' '.join(str(i) for i in key_words_list)
##############################################################################

#获取out文件中总的结构数
with open(out_file_name, 'r') as out_file:
    out_file_lines = out_file.readlines()
Step_number = []
for line in out_file_lines:
    if line.find('Step number') > -1:
        Step_number.append(int(line.split()[2]))
max_num_structs_in_out = max(Step_number)
#########################################################

def file_processor(path):

    #确定需要选取的结构所在行数
    
    pick_struct_begain_line = 4 + (12 + num_atoms) * (max_num_structs_in_out - num_pick_structs)

    with open(path, 'r')  as f:
        txt = f.readlines()

        cluster_Coordinates = []
        lines_input_1 = ''
        with open('{}/Edit_info_file/Read_GA_out_xyz_for_gjf.txt'.format(current_path), 'r') as fread:
            for line in fread.readlines()[10:]:  # in this case, let's copy line 6-结尾
                lines_input_1 += line

        if linecache.getline(path, pick_struct_begain_line + 1).split()[0] == 'Standard' and linecache.getline(path, pick_struct_begain_line + 1).split()[1] == 'orientation:':

            for struct_num_i in range(num_pick_structs):
                xyz_begain_line = pick_struct_begain_line + (12 + num_atoms) * struct_num_i + 6
                xyz_end_line      = xyz_begain_line + num_atoms

                if (num_pick_structs - struct_num_i) <= Final_GAUSSIAN_run_3_max_serial_num and (num_pick_structs - struct_num_i) > Final_GAUSSIAN_run_2_max_serial_num:
                    with open(Final_GAUSSIAN_run_3, 'a+') as f_GAUSSIAN_run_3:
                        f_GAUSSIAN_run_3.write('{} {}-M{}-{}.gjf\n'.format(Gaussian_type, cluster_name, multi, num_pick_structs - struct_num_i))

                if (num_pick_structs - struct_num_i) <= Final_GAUSSIAN_run_2_max_serial_num and (num_pick_structs - struct_num_i) > Final_GAUSSIAN_run_1_max_serial_num:
                    with open(Final_GAUSSIAN_run_2, 'a+') as f_GAUSSIAN_run_2:
                        f_GAUSSIAN_run_2.write('{} {}-M{}-{}.gjf\n'.format(Gaussian_type, cluster_name, multi, num_pick_structs - struct_num_i))

                if (num_pick_structs - struct_num_i) <= Final_GAUSSIAN_run_1_max_serial_num:
                    with open(Final_GAUSSIAN_run_1, 'a+') as f_GAUSSIAN_run_1:
                        f_GAUSSIAN_run_1.write('{} {}-M{}-{}.gjf\n'.format(Gaussian_type, cluster_name, multi, num_pick_structs - struct_num_i))

                if len(each_ele_type_num_without_0) == 1:
                    ele_1 = ele_atomic_num_dict[linecache.getline(path, xyz_begain_line).split()[1]]

                if len(each_ele_type_num_without_0) == 2:
                    ele_1 = ele_atomic_num_dict[linecache.getline(path, xyz_begain_line).split()[1]]
                    ele_2 = ele_atomic_num_dict[linecache.getline(path, xyz_begain_line + each_ele_type_num_without_0[0]).split()[1]]

                if len(each_ele_type_num_without_0) == 3:
                    ele_1 = ele_atomic_num_dict[linecache.getline(path, xyz_begain_line).split()[1]]
                    ele_2 = ele_atomic_num_dict[linecache.getline(path, xyz_begain_line + each_ele_type_num_without_0[0]).split()[1]]
                    ele_3 = ele_atomic_num_dict[linecache.getline(path, xyz_begain_line + each_ele_type_num_without_0[0] + each_ele_type_num_without_0[1]).split()[1]]

                output_file = os.path.join(os.getcwd(), '{}-M{}-{}.gjf'.format(cluster_name, multi, num_pick_structs - struct_num_i))
                with open(output_file, 'w') as f:
                    f.truncate()
                with open(output_file, 'a+') as fwrite:
                    fwrite.write('%chk={}-M{}-{}.chk\n'.format(cluster_name, multi, num_pick_structs - struct_num_i))
                    fwrite.write('%mem={}MB\n'.format(memory))
                    fwrite.write('%nprocshared={}\n'.format(nprocshared))
                    fwrite.write('# {}/{}  opt(MaxCycle={})  {}\n'.format(functional, basis_set, opt_maxcycle, key_words))
                    fwrite.write('\n')
                    fwrite.write('Title Card Required\n')
                    fwrite.write('\n')
                    fwrite.write('{}  {}\n'.format(charge, multi))
                    cluster_Coordinates = []
                    for line_line in range(xyz_begain_line,xyz_end_line):
                        data1 = linecache.getline(path, line_line).split()[3]  # 坐标x
                        data2 = linecache.getline(path, line_line).split()[4]  # 坐标y
                        data3 = linecache.getline(path, line_line).split()[5]  # 坐标z
                        cluster_Coordinates = cluster_Coordinates + [float(data1)] + [float(data2)] + [float(data3)]
                    xyz = np.array(cluster_Coordinates).reshape(-1, num_atoms, 3)
                    for j, cords in enumerate(xyz):
                        m = 1
                        cords = cords.tolist()
                        for cord in cords:
                            if len(each_ele_type_num_without_0) == 1:
                                fwrite.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.format('{}'.format(ele_1), cord[0], cord[1], cord[2]))

                            if len(each_ele_type_num_without_0) == 2:
                                if m <= each_ele_type_num_without_0[0]:
                                     fwrite.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.format('{}'.format(ele_1), cord[0], cord[1], cord[2]))
                                     m = m + 1
                                else:
                                    fwrite.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.format('{}'.format(ele_2), cord[0], cord[1], cord[2]))

                            if len(each_ele_type_num_without_0) == 3:
                                if m <= each_ele_type_num_without_0[0]:
                                     fwrite.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.format('{}'.format(ele_1), cord[0], cord[1], cord[2]))
                                     m = m + 1
                                elif m > each_ele_type_num_without_0[0] and m <= (each_ele_type_num_without_0[0] + each_ele_type_num_without_0[1]):
                                    fwrite.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.format('{}'.format(ele_2), cord[0], cord[1], cord[2]))
                                    m = m + 1
                                elif m >   (each_ele_type_num_without_0[0] + each_ele_type_num_without_0[1]):
                                    fwrite.write('{:<5s}{:15.8f}{:15.8f}{:15.8f}\n'.format('{}'.format(ele_3), cord[0], cord[1], cord[2]))
                    fwrite.write('\n')
                    fwrite.write(lines_input_1)


def generate_result():
     path = '{}/{}'.format(current_path,out_file_name)
     file_processor(path)



if __name__ == "__main__":
    generate_result()
