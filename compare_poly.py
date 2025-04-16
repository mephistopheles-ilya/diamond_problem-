import subprocess
import os
import time

"""
models = ["Reflect_00001", "Reflect_10001", "Reflect_20001", "Reflect_30002"
         ,"Reflect_50001", "Reflect_60001", "Reflect_70001"]

for model in models:
    file1 = "ChangedModels/" + model + ".obj"
    file2 = "/home/ilya/Downloads/Sample21Test/" + model + "reb" + ".obj"
    p = subprocess.run(f'./cmp_poly.out --type_of_file=obj --file1={file1} --file2={file2}', capture_output=True, text=True, shell=True)
    print(p.stdout)
"""


folder_path = '/home/ilya/Downloads/unzipping'

distance_data = []
slope_data = []
azimuth_data = []

t1 = time.time()
for filename in os.listdir(folder_path):
    file1 = "ChangedModels/Reflect_20001.obj"
    file_path = os.path.join(folder_path, filename)
    p = subprocess.run(f'./cmp_poly.out --type_of_file=obj --file1={file1} --file2={file_path}', capture_output=True, text=True, shell=True)
    ps = p.stdout.split("\n\n")
    for ans in ps:
        if "Distance" in ans:
            distance_data.append(ans)
        elif "Slope" in ans:
            slope_data.append(ans)
        elif "Azimuth" in ans:
            azimuth_data.append(ans)
t2 = time.time()
print("Time :", t2 - t1)

t1 = time.time()
distance_data.sort(key=lambda line: float(line.split("\n")[3].split(" ")[3]), reverse=True)
slope_data.sort(key=lambda line: float(line.split("\n")[3].split(" ")[3]), reverse=True)
azimuth_data.sort(key=lambda line: float(line.split("\n")[3].split(" ")[3]), reverse=True)
t2 = time.time()
print("Time :", t2 - t1)

print(distance_data[0].split("\n")[3].split(" ")[3])
print(slope_data[0].split("\n")[3].split(" ")[3])
print(slope_data[0].split("\n")[3].split(" ")[3])

with open("Reflect_20001_distance_log.txt", encoding="utf8", mode="w") as f:
    for line in distance_data:
        print(line, end="\n\n", file=f)

with open("Reflect_20001_slope_log.txt", encoding="utf8", mode="w") as f:
    for line in slope_data:
        print(line, end="\n\n", file=f)

with open("Reflect_20001_azimuth_log.txt", encoding="utf8", mode="w") as f:
    for line in azimuth_data:
        print(line, end="\n\n", file=f)



