import os
import itertools as it

os.makedirs("ChangedModels", exist_ok=True)

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
    for sigma_slope, sigma_azimuth, sigma_shift in it.product([1, 2], repeat=3):
        for is_change_slope, is_change_azimuth, is_change_shift in it.product([0, 1], repeat=3):
            if (is_change_slope, is_change_azimuth, is_change_azimuth) != (0, 0, 0) :
                for is_change_pavilion, is_change_crown in it.product([0, 1], repeat=2):
                    if (is_change_pavilion, is_change_crown) != (0, 0) :
                        parametrs = dict()
                        parametrs["--input_file"] = "InitialModels/" + file + ".txt"
                        new_file = file + "_f"
                        if is_change_crown == True:
                            new_file += "c"
                        if is_change_pavilion == True:
                            new_file += "p"
                        if is_change_slope == True:
                            new_file += "s"
                            new_file += str(sigma_slope )
                        if is_change_azimuth == True:
                            new_file += "a"
                            new_file += str(sigma_azimuth)
                        if is_change_shift == True:
                            new_file += "d"
                            new_file += str(sigma_shift)
                        new_file += "G"
                        parametrs["--type_of_input_file"] = "txt"
                        parametrs["--output_file"] = "ChangedModels/" + new_file
                        parametrs["--type_of_input_file"] = "txt obj ply"
                        parametrs["--is_change_slope"] = str(is_change_slope)
                        parametrs["--is_change_azimuth"] = str(is_change_azimuth)
                        parametrs["--is_change_shift"] = str(is_change_shift)
                        parametrs["--is_change_pavilion"] = str(is_change_pavilion)
                        parametrs["--is_change_crown"] = str(is_change_crown)
                        parametrs["--sigma_of_slope"] = str(sigma_slope)
                        parametrs["--sigma_of_azimuth"] = str(sigma_azimuth)
                        parametrs["--sigma_of_shift"] = str(sigma_shift / 1000)
                        all_paremtrs.append(parametrs)

for parametrs in all_paremtrs:
    for key, value in parametrs.items():
        print(key, value)
    print()

print(len(all_paremtrs))





