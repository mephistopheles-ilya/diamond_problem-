import os
import itertools as it
import subprocess

DIR_OUT = "new_ChangedModels"

os.makedirs(DIR_OUT, exist_ok=True)

input_files = ["Reflect_00001"]
'''
        "Reflect_10001",
        "Reflect_20001", 
        "Reflect_30002",
        "Reflect_40001",
        "Reflect_50001",
        "Reflect_60001",
        "Reflect_70001"]
'''



all_paremtrs = []

for file in input_files:
    for sigma_slope, sigma_azimuth, sigma_shift in it.permutations({0, 1, 2, 3}, 3):
        if (sigma_slope, sigma_azimuth, sigma_shift) != (0, 0, 0):
            for is_change_pavilion, is_change_crown in it.permutations({0, 1}, 2):
                if (is_change_pavilion, is_change_crown) != (0, 0) :
                    parametrs = dict()
                    parametrs["--input_file"] = "InitialModels/" + file + ".txt"
                    new_file = file + "_f"
                    if is_change_crown == True:
                        new_file += "c"
                    if is_change_pavilion == True:
                        new_file += "p"
                    if sigma_slope > 0:
                        new_file += "s"
                        new_file += str(sigma_slope )
                    if sigma_azimuth > 0 == True:
                        new_file += "a"
                        new_file += str(sigma_azimuth)
                    if sigma_shift > 0 == True:
                        new_file += "d"
                        new_file += str(sigma_shift)
                    new_file += "G"
                    parametrs["--type_of_input_file"] = "txt"
                    parametrs["--output_file"] = DIR_OUT + new_file
                    parametrs["--type_of_input_file"] = "txt obj ply"
                    parametrs["--is_change_slope"] = str(bool(sigma_slope))
                    parametrs["--is_change_azimuth"] = str(bool(sigma_azimuth))
                    parametrs["--is_change_shift"] = str(bool(sigma_shift))
                    parametrs["--is_change_pavilion"] = str(is_change_pavilion)
                    parametrs["--is_change_crown"] = str(is_change_crown)
                    parametrs["--sigma_slope"] = str(sigma_slope)
                    parametrs["--sigma_azimuth"] = str(sigma_azimuth)
                    parametrs["--sigma_shift"] = str(sigma_shift / 1000)
                    parametrs["--distribution"] = "normal"
                    all_paremtrs.append(parametrs)

for file in input_files:
    if file != "Reflect_40001":
        for procent in 98, 99, 101, 102:
            for is_change_height_crown, is_change_height_pavilion in it.permutations({0, 1}, 2):
                if (is_change_height_crown, is_change_height_pavilion) != (0, 0) and (is_change_height_crown, is_change_height_pavilion) != (1, 1):
                    parametrs = dict()
                    parametrs["--input_file"] = "InitialModels/" + file + ".txt"
                    new_file = file + "_e"
                    if is_change_height_crown == 1 :
                        new_file += "c"
                        new_file += str(int(100 - procent))
                    if is_change_height_pavilion == 1:
                        new_file += "p"
                        new_file += str(int(100 - procent))
                    parametrs["--type_of_input_file"] = "txt"
                    parametrs["--output_file"] = DIR_OUT + new_file
                    parametrs["--type_of_input_file"] = "txt obj ply"
                    parametrs["--is_change_height_crown"] = str(int(is_change_height_crown))
                    parametrs["--is_change_height_pavilion"] = str(int(is_change_height_pavilion))
                    parametrs["--parametr_of_height_crown"] = str(int(procent))
                    parametrs["--parametr_of_height_pavilion"] = str(int(procent))
                    all_paremtrs.append(parametrs)


for parametrs in all_paremtrs:
    command = ["./transform_model.out"] + [f"{key}={value}" for key, value in parametrs.items()]
    print("Running:", " ".join(command))

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stdout:
        print("stdout:", result.stdout)
    if result.stderr:
        print("stderr:", result.stderr)

    if result.returncode != 0:
        print(f"Error: Program returned non-zero exit code {result.returncode}")
        break
    print("----------------------------------------")





