import subprocess
import os
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


#folder_paths = ["unzipping0", "unzipping1", "unzipping2", "unzipping3", "unzipping4",  "unzipping5",  "unzipping6",  "unzipping7"]
#model_names = ["Reflect_00001", "Reflect_10001", "Reflect_20001", "Reflect_30001", "Reflect_40001", "Reflect_50001", "Reflect_60001", "Reflect_70001"] 

folder_paths = ["unzipping7"]
model_names = ["Reflect_70001"] 

max_deviation = 0

for folder_path, model_name in zip(folder_paths, model_names):
    folder_path = '/home/ilya/Downloads/' + folder_path

    angle_data = []
    slope_data = []
    azimuth_data = []

    t1 = time.time()
    for filename in os.listdir(folder_path):
        if ".txt" not in filename:
            file_path = os.path.join(folder_path, filename)
            filename = filename[0:filename.rfind("_")] + ".obj";
            file1 = "new_ChangedModels/" + filename
            if filename == "Reflect_40001_fcs2a2d2G.obj":
                continue;
            print(filename)
            print(f'./cmp_poly.out --type_of_file=obj --file1={file1} --file2={file_path}')
            p = subprocess.run(f'./cmp_poly.out --type_of_file=obj --file1={file1} --file2={file_path}', capture_output=True, text=True, shell=True)
            ps = p.stdout.split("\n\n")
            print(p.stderr)
            for ans in ps:
                if "Angle" in ans:
                    angle_data.append(ans)
                elif "Slope" in ans:
                    slope_data.append(ans)
                elif "Azimuth" in ans:
                    azimuth_data.append(ans)

    t2 = time.time()
    print("Time :", t2 - t1)

    angle_values = [float(data.split("\n")[5].split(" ")[3]) for data in angle_data]
    slope_values = [float(data.split("\n")[5].split(" ")[3]) for data in slope_data]
    azimuth_values = [float(data.split("\n")[5].split(" ")[3]) for data in azimuth_data]

    max_deviation = max(max_deviation, max(angle_values), max(slope_values), max(azimuth_values))

    edges = [int(data.split("\n")[2].split(" ")[2]) for data in angle_data]
    faces = [int(data.split("\n")[3].split(" ")[2]) for data in slope_data]


    with PdfPages(model_name + '_compare.pdf') as pdf:
        plt.figure(figsize=(10, 5))

        x_values = range(len(angle_values))

        plt.plot(x_values, angle_values, label='Angle in degrees', linewidth=1, color="green")
        plt.plot(x_values, slope_values, label='Slope in degrees', linewidth=1, color="red")
        plt.plot(x_values, azimuth_values, label='Azimuth in degrees', linewidth=1, color="blue")

        plt.title(model_name)
        plt.xlabel('number in log file')
        plt.ylabel('degrees')

        plt.xticks(range(0, len(angle_values), 5), fontsize=4, rotation=90)

        plt.legend()
        pdf.savefig(bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(10, 5))
        
        plt.plot(x_values, edges, label="Edges", linewidth=1, color="blue")
        plt.plot(x_values, faces, label="Faces", linewidth=1, color="red")
        
        plt.title(model_name)
        plt.xlabel('number in log file (angle)')
        plt.ylabel('initial - rebuild')
        
        plt.xticks(range(0, len(angle_values), 5), fontsize=4, rotation=90)

        plt.legend()
        pdf.savefig(bbox_inches='tight')  
        plt.close()

    with open(model_name + "_angle_log.txt", encoding="utf8", mode="w") as f:
        number = 0
        for line in angle_data:
            print("NUMBER =", number, file=f)
            print(line, end="\n\n", file=f)
            number += 1

    with open(model_name + "_slope_log.txt", encoding="utf8", mode="w") as f:
        number = 0
        for line in slope_data:
            print("NUMBER =", number, file=f)
            print(line, end="\n\n", file=f)
            number += 1

    with open(model_name + "_azimuth_log.txt", encoding="utf8", mode="w") as f:
        number = 0
        for line in azimuth_data:
            print("NUMBER =", number, file=f)
            print(line, end="\n\n", file=f)
            number += 1


print(f'Max deviation = {max_deviation}')
