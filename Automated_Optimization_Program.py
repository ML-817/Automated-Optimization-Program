#!/usr/bin python
# coding: utf-8

import os
import shutil
import re
import pandas as pd
from collections import defaultdict

import statsmodels.datasets


def calculate_min_multiplicity(total_electrons):
    if total_electrons % 2 == 1:
        min_multiplicity = 1
    else:
        min_multiplicity = 2
    
    return min_multiplicity

def extract_elements_and_charge(cluster_name):
    pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    elements = pattern.findall(cluster_name)

    element_counts = {}
    total_charge = 0
    
    for element, count in elements:
        if count:
            count = int(count)
        else:
            count = 1

        if cluster_name[-1] in ['+', '-']:
            charge = int(cluster_name[-1] + '1')
            # total_charge += charge
            # element = element[:-1]
        else:
            charge = 0

        if element in element_counts:
            element_counts[element] += count
        else:
            element_counts[element] = count

    return element_counts, charge

def calculate_total_electrons(element_counts):
    total_electrons = 0
    for element, count in element_counts.items():
        atomic_number = element_atomic_numbers.get(element, 0)
        total_electrons += atomic_number * count
    
    return atomic_number, total_electrons

def get_atomic_numbers(element_counts, element_atomic_numbers):
    elements_info = {}
    for element, count in element_counts.items():
        atomic_number = element_atomic_numbers.get(element)
        if atomic_number is not None:
            elements_info[element] = {'atomic_number': atomic_number, 'count': count}

    return elements_info

def create_folders_and_copy_files(cluster_name, elements_info, multiplicity, charge, pretrain_cluster,
                                  pretrain_cluster_multi, Optimization_type, Num_initial_structs_for_DNN_TL,
                                  Num_initial_structs_for_DNN_TL_GA, DFT_steps):
    atom_counts = defaultdict(int)
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', cluster_name)
    element_list = []
    for element, count in elements:
        if count == '':
            count = 1
        else:
            count = int(count)
        element_list.append(element)
        atom_counts[element] += count
    atom_counts = dict(atom_counts)
    Natoms = sum(atom_counts.values())
    current_path = os.getcwd()

    if Optimization_type == 'GA':
        try:
            path = '{}/GA_{}/M{}'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(path, exist_ok=True)

            Input_info_format_file_path = '{}/Input_info_format_file'.format(os.getcwd())
            print(os.getcwd())
            Main_program_path = '{}/Main_program'.format(os.getcwd())
            Program_sub_script_path = '{}/Program_sub_script'.format(os.getcwd())
            Program_sub_script_for_cluster_path = '{}/GA_{}/M{}/Program_sub_script'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(Program_sub_script_for_cluster_path, exist_ok=True)

            Edit_info_file_path = '{}/GA_{}/M{}/Edit_info_file'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(Edit_info_file_path, exist_ok=True)

            files_to_copy = ["GA_input.dat"]
            for file_to_copy in files_to_copy:
                file_path = os.path.join(Input_info_format_file_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, Edit_info_file_path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            Edit_files_to_copy = ["GA_out_xyz_for_gjf.txt"]
            for file_to_copy in Edit_files_to_copy:
                file_path = os.path.join(Input_info_format_file_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, Edit_info_file_path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            files_to_copy = ["GA.exe"]
            for file_to_copy in files_to_copy:
                file_path = os.path.join(Main_program_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            Program_sub_script_files_to_copy = ["Get_seeds_gjf.py", "Get_seeds_log.py", "BLDA_method_for_generating_initial_seed.json", "edit_json_for_BLDA_method.py", "Read_GA_out_xyz_for_gjf.py"]
            for file_to_copy in Program_sub_script_files_to_copy:
                file_path = os.path.join(Program_sub_script_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, Program_sub_script_for_cluster_path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            line_input_1 = ''
            line_input_2 = ''
            with open('{}/Input_info_format_file/GA_input.dat'.format(os.getcwd()), 'r') as fread:
                lines_txt = fread.readlines()[:]
            for line in lines_txt[3:29]:
                line_input_1 += line
            for line in lines_txt[30:]:
                line_input_2 += line

            dat_file_path = os.path.join(Edit_info_file_path, "GA_input.dat")
            with open(dat_file_path, 'w') as data_file:
                data_file.truncate()
            with open(dat_file_path, 'a+') as data_file:
                if len(element_list) == 3:
                    data_file.write(' {}   {}   {}                     *** Metal 1, M1 | Format(A3) free for 4-20     \n'.format(elements[0][0], elements_info[elements[0][0]]['atomic_number'], elements_info[elements[0][0]]['count']))
                    data_file.write(' {}   {}   {}                     *** Metal 2, M2 | For three elements     \n'.format(elements[1][0], elements_info[elements[1][0]]['atomic_number'], elements_info[elements[1][0]]['count']))
                    data_file.write(' {}   {}   {}                     *** Oxygen, O   | in general        \n'.format(elements[2][0], elements_info[elements[2][0]]['atomic_number'], elements_info[elements[2][0]]['count']))
                elif len(element_list) == 2:
                    data_file.write(' {}   {}   {}                     *** Metal 1, M1 | Format(A3) free for 4-20     \n'.format(elements[0][0], elements_info[elements[0][0]]['atomic_number'], elements_info[elements[0][0]]['count']))
                    data_file.write(' XX   XX   0                     *** Metal 2, M2 | For three elements     \n')
                    data_file.write(' {}   {}   {}                     *** Oxygen, O   | in general        \n'.format(elements[1][0], elements_info[elements[1][0]]['atomic_number'], elements_info[elements[1][0]]['count']))
                elif len(element_list) == 1:
                    data_file.write(' XX   XX   0                     *** Metal 1, M1 | Format(A3) free for 4-20     \n')
                    data_file.write(' XX   XX   0                     *** Metal 2, M2 | For three elements     \n')
                    data_file.write(' {}   {}   {}                     *** Oxygen, O   | in general        \n'.format(elements[0][0], elements_info[elements[0][0]]['atomic_number'], elements_info[elements[0][0]]['count']))
                data_file.write(line_input_1)
                data_file.write('{}  {}    \n'.format(charge, multiplicity))
                data_file.write(line_input_2)
                data_file.write('\n')
                data_file.write('\n')

            os.chdir(path)
            os.system('./GA.exe')
            os.chdir(current_path)
        except Exception as e:
             print(f"Error creating folder, copying files, or modifying gene-v05.dat: {e}")

    elif Optimization_type == 'DNN_TL':
        try:
            path = '{}/DNN_TL_{}/M{}'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(path, exist_ok=True)

            Input_info_format_file_path = '{}/Input_info_format_file'.format(os.getcwd())
            Edit_info_file_path = '{}/DNN_TL_{}/M{}/Edit_info_file'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(Edit_info_file_path, exist_ok=True)

            Main_program_path = '{}/Main_program'.format(os.getcwd())
            Program_sub_script_path = '{}/Program_sub_script'.format(os.getcwd())

            files_to_copy = ["DNN_TL_Auto_Run"]
            for file_to_copy in files_to_copy:
                file_path = os.path.join(Main_program_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            files_to_copy = [ "DNN_TL_run_get_gjf.py", "DNN_TL_run_get_DNN_opt_structs.py", "DNN_TL_get_structs_and_energies_from_Gaussian_log.py", "AUTO_Processing_final_DFT_opted_structs.py", "Read_log_last_struct_xyz_M_en_for_filter.py", "Filter-for-DFT.json", "Draw_cluster.json", "Processing_finel_DFT_opted_structs.py"]
            for file_to_copy in files_to_copy:
                file_path = os.path.join(Program_sub_script_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            Edit_info_files_to_copy = ["DNN_TL_input_for_gjf.txt", "DNN_TL_input_for_opt_structs.txt"]
            for file_to_copy in Edit_info_files_to_copy:
                file_path = os.path.join(Input_info_format_file_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, Edit_info_file_path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            line_input_1 = ''
            with open('{}/DNN_TL_input_for_gjf.txt'.format(Input_info_format_file_path), 'r') as fread:
                lines_txt = fread.readlines()[:]
            for line in lines_txt[5:]:
                line_input_1 += line
                #print('***************************************************************************')

            dat_file_path = os.path.join(Edit_info_file_path, "DNN_TL_input_for_gjf.txt")
            with open(dat_file_path, 'w') as data_file:
                data_file.truncate()
            with open(dat_file_path, 'a+') as data_file:
                data_file.write("# information for DNN_TL_run_get_gjf.py#\n")
                data_file.write('Cluster                =   {}   \n'.format(cluster_name))
                data_file.write('Min_Multi              =   {}    \n'.format(multiplicity))
                data_file.write('Max_Multi              =   {}    \n'.format(multiplicity + 4))
                if len(element_list) == 1:
                   data_file.write('Num_initial_structs    =   {},{}    \n'.format(Num_initial_structs_for_DNN_TL['{}'.format(Natoms)],
                                                                                   0.5 * Num_initial_structs_for_DNN_TL['{}'.format(Natoms)]))
                elif len(element_list) == 2:
                    data_file.write('Num_initial_structs    =   {},{}    \n'.format(1.5 * Num_initial_structs_for_DNN_TL['{}'.format(Natoms)],
                                                                                    0.5 * 1.5* Num_initial_structs_for_DNN_TL['{}'.format(Natoms)]))
                elif len(element_list) == 3:
                    data_file.write('Num_initial_structs    =   {},{}    \n'.format(2 * Num_initial_structs_for_DNN_TL['{}'.format(Natoms)], 0.5 * 2 * Num_initial_structs_for_DNN_TL['{}'.format(Natoms)]))
                data_file.write('DFT_steps              =   {}\n'.format(DFT_steps['{}'.format(Natoms)]))
                data_file.write(line_input_1)
                data_file.write('\n')
                data_file.write('\n')

            line_input_2 = ''
            line_input_3 = ''
            line_input_4 = ''
            line_input_5 = ''
            with open('{}/DNN_TL_input_for_opt_structs.txt'.format(Input_info_format_file_path), 'r') as fread:
                lines_txt = fread.readlines()[:]
            for line in lines_txt[5:15]:
                line_input_2 += line
            for line in lines_txt[17:37]:
                line_input_3 += line
            line_input_4 = lines_txt[39]
            for line in lines_txt[41:]:
                line_input_5 += line

            dat_file_path_2 = os.path.join(Edit_info_file_path, "DNN_TL_input_for_opt_structs.txt")
            with open(dat_file_path_2, 'w') as data_file_2:
                data_file_2.truncate()
            with open(dat_file_path_2, 'a+') as data_file_2:
                data_file_2.write("# information for DNN_TL_run_get_DNN_opt_structs.py #\n")
                data_file_2.write('Cluster                =   {}   \n'.format(cluster_name))
                data_file_2.write('Min_Multi              =   {}    \n'.format(multiplicity))
                data_file_2.write('Max_Multi              =   {}    \n'.format(multiplicity + 4))
                data_file_2.write(line_input_2)
                data_file_2.write('Offer_pretrain_NN_cluster            =  {}\n'.format(pretrain_cluster))
                data_file_2.write('Pretrain_cluster_Min_Multi        =  {}\n'.format(pretrain_cluster_multi))
                data_file_2.write(line_input_3)
                data_file_2.write('Final_NN_Multi             =   {}\n'.format(multiplicity))
                data_file_2.write('Which_Multi_Transfer_to_Final_NN             =   {}\n'.format(multiplicity + 4))
                data_file_2.write(line_input_4)
                data_file_2.write('Ref_NN_Multi             =   {}\n'.format(multiplicity))
                data_file_2.write(line_input_5)
                data_file_2.write('\n')
                data_file_2.write('\n')

            os.chdir(path)
            os.system('./DNN_TL_Auto_Run')
            os.chdir(current_path)
        except Exception as e:
             print(f"Error creating folder, copying files, or modifying gene-v05.dat: {e}")

    elif Optimization_type == 'DNN_TL_GA':
        try:
            path = '{}/DNN_TL_GA_{}/M{}'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(path, exist_ok=True)

            Input_info_format_file_path = '{}/Input_info_format_file'.format(os.getcwd())
            Main_program_path = '{}/Main_program'.format(os.getcwd())

            Edit_info_file_path = '{}/DNN_TL_GA_{}/M{}/Edit_info_file'.format(os.getcwd(), cluster_name, multiplicity)
            os.makedirs(Edit_info_file_path, exist_ok=True)

            files_to_copy = ["DNN_TL_GA.exe"]
            for file_to_copy in files_to_copy:
                file_path = os.path.join(Main_program_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            folder_to_copy = "Program_sub_script"
            folder_path = os.path.join(os.getcwd(), folder_to_copy)
            shutil.copytree(folder_path, "{}/Program_sub_script".format(path))

            Edit_files_to_copy = ["DNN_TL_GA_input.dat", "DNN_TL_GA_out_xyz_for_gjf.txt"]
            for file_to_copy in Edit_files_to_copy:
                file_path = os.path.join(Input_info_format_file_path, file_to_copy)
                if os.path.exists(file_path):
                    shutil.copy(file_path, Edit_info_file_path)
                else:
                    print(f"File '{file_to_copy}' not found in the current directory.")

            line_input_1 = ''
            line_input_3 = ''
            line_input_5 = ''
            line_input_6 = ''
            line_input_7 = ''
            with open('{}/Input_info_format_file/DNN_TL_GA_input.dat'.format(os.getcwd()), 'r') as fread:
                lines_txt = fread.readlines()[:]
            for line in lines_txt[3:20]:
                line_input_1 += line

            line_input_2 = lines_txt[21].split()[0:3]

            for line in lines_txt[22:24]:
                line_input_3 += line

            line_input_4 = lines_txt[24]

            for line in lines_txt[25:28]:
                line_input_5 += line

            line_input_6 = lines_txt[28].split()

            for line in lines_txt[30:]:
                line_input_7 += line

            #3 3 3                Largest pairs of parents    分配三次DFT计算时的父代最大成对数
            number = Num_initial_structs_for_DNN_TL_GA['{}'.format(Natoms)] / 5
            for pairs_1 in range(int(number // 2), 1, -1):
                remainder = number - pairs_1
                if remainder % 2 == 0:
                    pairs_2 = pairs_3 = int(remainder // 2)

            dat_file_path = os.path.join(Edit_info_file_path, "DNN_TL_GA_input.dat")
            with open(dat_file_path, 'w') as data_file:
                data_file.truncate()
            with open(dat_file_path, 'a+') as data_file:
                print(len(element_list))
                if len(element_list) == 3:
                    data_file.write(' {}   {}   {}                     *** Metal 1, M1 | Format(A3) free for 4-20     \n'.format(elements[0][0], elements_info[elements[0][0]]['atomic_number'], elements_info[elements[0][0]]['count']))
                    data_file.write(' {}   {}   {}                     *** Metal 2, M2 | For three elements     \n'.format(elements[1][0], elements_info[elements[1][0]]['atomic_number'], elements_info[elements[1][0]]['count']))
                    data_file.write(' {}   {}   {}                     *** Oxygen, O   | in general        \n'.format(elements[2][0], elements_info[elements[2][0]]['atomic_number'], elements_info[elements[2][0]]['count']))
                elif len(element_list) == 2:
                    data_file.write(' {}   {}   {}                     *** Metal 1, M1 | Format(A3) free for 4-20     \n'.format(elements[0][0], elements_info[elements[0][0]]['atomic_number'], elements_info[elements[0][0]]['count']))
                    data_file.write(' XX   XX   0                     *** Metal 2, M2 | For three elements     \n')
                    data_file.write(' {}   {}   {}                     *** Oxygen, O   | in general        \n'.format(elements[1][0], elements_info[elements[1][0]]['atomic_number'], elements_info[elements[1][0]]['count']))
                elif len(element_list) == 1:
                    data_file.write(' XX   XX   0                     *** Metal 1, M1 | Format(A3) free for 4-20     \n')
                    data_file.write(' XX   XX   0                     *** Metal 2, M2 | For three elements     \n')
                    data_file.write(' {}   {}   {}                     *** Oxygen, O   | in general        \n'.format(elements[0][0], elements_info[elements[0][0]]['atomic_number'], elements_info[elements[0][0]]['count']))
                data_file.write(line_input_1)
                data_file.write('{} {} {}                Largest pairs of parents    \n'.format(pairs_1, pairs_2, pairs_3))
                data_file.write('{} {} {} {} M{}\n'.format(line_input_2[0], line_input_2[1], line_input_2[2], pretrain_cluster, pretrain_cluster_multi))
                data_file.write(line_input_3)
                data_file.write('{}  {}  {}  M{}  {}  {}  {}\n'.format(line_input_4.split()[0], line_input_4.split()[1], cluster_name, multiplicity, line_input_4.split()[4], line_input_4.split()[5], line_input_4.split()[6]))
                data_file.write(line_input_5)
                data_file.write('{} Opt(MaxCycle={})  {}\n'.format(line_input_6[0], DFT_steps[Natoms], line_input_6[2:]))
                data_file.write('{}  {}    \n'.format(charge, multiplicity))
                data_file.write(line_input_7)
                data_file.write('\n')
                data_file.write('\n')
            os.chdir(path)
            os.system('./DNN_TL_GA.exe')
            os.chdir(current_path)
        except Exception as e:
             print(f"Error creating folder, copying files, or modifying gene-v05.dat: {e}")

element_atomic_numbers = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6,
    'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
    'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ga': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
    'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
    'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
    'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86
}

Num_initial_structs_for_DNN_TL = {'5':50, '6':50, '7':50, '8':50, '9':50, '10':60, '11':90, '12':110, '13':80, '14':90}
Num_initial_structs_for_DNN_TL_GA = {'5':30, '6':30, '7':30, '8':30, '9':30, '10':45, '11':60, '12':90, '13':60, '14':60, '15':75, '16':90, '17':90}
DFT_steps = {'5':10, '6':10, '7':10, '8':10, '9':10, '10':10, '11':10, '12':10, '13':20, '14':20, '15':20, '16':20, '17':20}

element_charge = {'+': 1, '-': -1}

Select_type_optimization   = 'GA'  # 可以选择【'GA'、'DNN_TL'、'DNN_TL_GA'】 这3种类型
if __name__ == "__main__":
    cluster_name_list = ["Rh2V6O1-", "Rh2V2O2-", "Rh2V2O3-"]
    pretrain_cluster = ''
    pretrain_cluster_multi = 0
    if Select_type_optimization == 'DNN_TL' or Select_type_optimization == 'DNN_TL_GA':
        Judge_if_first_run = 0
    for cluster_name in cluster_name_list:
        elements_and_charge_result = extract_elements_and_charge(cluster_name)
        total_electrons_result = calculate_total_electrons(elements_and_charge_result[0])
        min_multiplicity_result = calculate_min_multiplicity(total_electrons_result[1])
        multiplicity_list = [min_multiplicity_result, min_multiplicity_result+2, min_multiplicity_result+4]
        atomic_numbers_result = get_atomic_numbers(elements_and_charge_result[0], element_atomic_numbers)

        for multiplicity in multiplicity_list:
            if Select_type_optimization == 'DNN_TL_GA':
                if Judge_if_first_run == 0:
                    with open('./Input_info_format_file/DNN_TL_GA_input.dat', 'r') as f:
                        txt = f.readlines()[21].split()[3:5]
                        pretrain_cluster = txt[0]
                        pretrain_cluster_multi = txt[1][-1]
                elif Judge_if_first_run == 1:
                    pretrain_cluster_multi = multiplicity_list[0]
                elif Judge_if_first_run == 2:
                    pretrain_cluster_multi = multiplicity_list[1]

            elif Select_type_optimization == 'DNN_TL':
                if Judge_if_first_run == 0:
                    with open('./Input_info_format_file/DNN_TL_input_for_opt_structs.txt', 'r') as f:
                        pretrain_cluster = f.readlines()[15].split()[2]
                    with open('./Input_info_format_file/DNN_TL_input_for_opt_structs.txt', 'r') as f:
                        pretrain_cluster_multi = f.readlines()[16].split()[2]

            create_folders_and_copy_files(cluster_name,  atomic_numbers_result,multiplicity,
                                          elements_and_charge_result[1], pretrain_cluster, pretrain_cluster_multi,
                                          Select_type_optimization, Num_initial_structs_for_DNN_TL,
                                          Num_initial_structs_for_DNN_TL_GA, DFT_steps)

            if Select_type_optimization == 'DNN_TL' or Select_type_optimization == 'DNN_TL_GA':
                Judge_if_first_run += 1
                pretrain_cluster = cluster_name
                pretrain_cluster_multi = multiplicity_list[-1]










