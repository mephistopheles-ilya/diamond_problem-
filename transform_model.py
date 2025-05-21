import os
import itertools as it
import subprocess
import time

#Elapsed : 686.2949342727661 for  1560
start = time.time()

DIR_OUT = "new_new_ChangedModels"

os.makedirs(DIR_OUT, exist_ok=True)

input_files = ["Reflect_30002"]




all_parametrs_model = []
all_parametrs_contour = []

for file in input_files:
    for sigma_slope, sigma_azimuth, sigma_shift in it.product({0, 1, 2, 3}, repeat=3):
        if (sigma_slope, sigma_azimuth, sigma_shift) != (0, 0, 0):
            for is_change_pavilion, is_change_crown in it.product({0, 1}, repeat=2):
                if (is_change_pavilion, is_change_crown) != (0, 0) :
                    parametrs_model = dict()
                    parametrs_model["--input_file"] = "InitialModels/" + file + ".txt"
                    new_file = file + "_f"
                    if is_change_crown == 1:
                        new_file += "c"
                    if is_change_pavilion == 1:
                        new_file += "p"
                    if sigma_slope > 0:
                        new_file += "s"
                        new_file += str(sigma_slope )
                    if sigma_azimuth > 0:
                        new_file += "a"
                        new_file += str(sigma_azimuth)
                    if sigma_shift > 0:
                        new_file += "d"
                        new_file += str(sigma_shift)
                    new_file += "G"
                    parametrs_model["--type_of_input_file"] = "txt"
                    parametrs_model["--output_file"] = DIR_OUT + "/" + new_file
                    parametrs_model["--type_of_output_file"] = "txt_obj_ply"
                    parametrs_model["--is_change_slope"] = str(int(bool(sigma_slope)))
                    parametrs_model["--is_change_azimuth"] = str(int(bool(sigma_azimuth)))
                    parametrs_model["--is_change_shift"] = str(int(bool(sigma_shift)))
                    parametrs_model["--is_change_pavilion"] = str(is_change_pavilion)
                    parametrs_model["--is_change_crown"] = str(is_change_crown)
                    parametrs_model["--sigma_slope"] = str(sigma_slope/10)
                    parametrs_model["--sigma_azimuth"] = str(sigma_azimuth/10)
                    parametrs_model["--sigma_shift"] = str(sigma_shift / 1000)
                    parametrs_model["--distribution"] = "normal"
                    all_parametrs_model.append(parametrs_model)

                    parametrs_contour = dict()
                    parametrs_contour["--init_file"] = DIR_OUT + "/" + new_file + ".txt"
                    parametrs_contour["--directory_out"] = DIR_OUT + "/" + new_file + "_grid"
                    os.makedirs(parametrs_contour["--directory_out"], exist_ok=True)
                    all_parametrs_contour.append(parametrs_contour)


for file in input_files:
    #if file != "Reflect_40001":
    for procent in 101, 102, 103:
        for is_change_height_crown, is_change_height_pavilion in it.product({0, 1}, repeat=2):
            if (is_change_height_crown, is_change_height_pavilion) != (0, 0) and (is_change_height_crown, is_change_height_pavilion) != (1, 1):
                parametrs = dict()
                parametrs["--input_file"] = "InitialModels/" + file + ".txt"
                new_file = file + "_e"
                if is_change_height_crown == 1 :
                    new_file += "c"
                    new_file += str(int(procent - 100))
                if is_change_height_pavilion == 1:
                    new_file += "p"
                    new_file += str(int(procent - 100))
                parametrs["--type_of_input_file"] = "txt"
                parametrs["--output_file"] = DIR_OUT + "/" + new_file
                parametrs["--type_of_output_file"] = "txt_obj_ply"
                parametrs["--is_change_height_crown"] = str(int(is_change_height_crown))
                parametrs["--is_change_height_pavilion"] = str(int(is_change_height_pavilion))
                parametrs["--parametr_of_height_crown"] = str(int(procent))
                parametrs["--parametr_of_height_pavilion"] = str(int(procent))
                all_parametrs_model.append(parametrs)

                parametrs_contour = dict()
                parametrs_contour["--init_file"] = DIR_OUT + "/" + new_file + ".txt"
                parametrs_contour["--directory_out"] = DIR_OUT + "/" + new_file + "_grid"
                os.makedirs(parametrs_contour["--directory_out"], exist_ok=True)
                all_parametrs_contour.append(parametrs_contour)


counter = 0
for parametrs_model, parametrs_contour in zip(all_parametrs_model, all_parametrs_contour):
    counter += 1
    command = ["./transform_model.out"] + [f"{key}={value}" for key, value in parametrs_model.items()]
    print("Running:", " ".join(command))

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stdout:
        print("stdout:", result.stdout)
    if result.stderr:
        print("stderr:", result.stderr)

    if result.returncode != 0:
        print(f"Error: Program returned non-zero exit code {result.returncode}")
        break


    command = ["./get.out"] + ["--projections=400"] + ["--init_file=" + parametrs_contour["--init_file"]]
    print("Running:", " ".join(command))

    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.stdout:
        print("stdout:", result.stdout)
    if result.stderr:
        print("stderr:", result.stderr)

    if result.returncode != 0:
        print(f"Error: Program returned non-zero exit code {result.returncode}")
        break

    command = ["./tran.out"] + ["--projections=400","--grid=1", " --gr_method=ceils", "--shift=0", " --dist_points=0.008",
            "--directory_out=" + parametrs_contour["--directory_out"]]
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

end = time.time()

print("Elapsed :", end - start, "for ", counter)



