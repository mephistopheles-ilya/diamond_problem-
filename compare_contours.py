import subprocess
import os
import time


#folder_paths = ["unzipping0", "unzipping1", "unzipping2", "unzipping3", "unzipping4",  "unzipping5",  "unzipping6",  "unzipping7"]
#model_names = ["Reflect_00001", "Reflect_10001", "Reflect_20001", "Reflect_30001", "Reflect_40001", "Reflect_50001", "Reflect_60001", "Reflect_70001"] 

folder_paths = ["unzipping0"]
model_names = ["Reflect_00001"] 


t1 = time.time()
for folder_path, model_name in zip(folder_paths, model_names):

    outs = []
    folder_path = '/home/ilya/Downloads/' + folder_path

    for filename in os.listdir(folder_path):
        if "txt" not in filename :
            file_path = os.path.join(folder_path, filename)
            file_txt = os.path.join(folder_path, filename[0:filename.rfind(".")] + ".txt")
            filename = filename[0:filename.rfind("_")] + ".txt"
            file1 = "new_ChangedModels/" + filename
            if filename == "Reflect_40001_fcs2a2d2G.txt":
                continue;
            print(filename)
            print(f'./obj_to_txt.out --file_obj={file_path} --convert=obj --file_txt={file_txt}')
            p = subprocess.run(f'./obj_to_txt.out --file_obj={file_path} --convert=obj --file_txt={file_txt}'
                , capture_output=True, text=True, shell=True)
            if p.stderr :
                print(p.stderr)

            print(f'./get.out --projections=400 --init_file={file1} --directory=init_examples1 --file_with_points_numbers=0')
            p = subprocess.run(f'./get.out --projections=400 --init_file={file1} --directory=init_examples1 --file_with_points_numbers=0'
                    , capture_output=True, text=True, shell=True)
            if p.stderr :
                print(p.stderr)

            print(f'./get.out --projections=400 --init_file={file_txt} --directory=init_examples2 --file_with_points_numbers=0')
            p = subprocess.run(f'./get.out --projections=400 --init_file={file_txt} --directory=init_examples2 --file_with_points_numbers=0'
                    , capture_output=True, text=True, shell=True)
            if p.stderr :
                print(p.stderr)

            print(f'./gorizontal_compare.out init_examples1 init_examples2 400')
            p = subprocess.run(f'./gorizontal_compare.out init_examples1 init_examples2 400'
                    , capture_output=True, text=True, shell=True)
            if p.stderr :
                print(p.stderr)
                with open("error_log.txt", encoding="utf8", mode="w") as f:
                    print(p.stderr, end="\n", file=f)



            ps = p.stdout.split("\n")
            ps = [filename] + ps
            print(ps)
            outs.append(ps)
            print("\n\n")

    outs.sort(key=lambda strings: float(strings[2].split()[2]), reverse=True)
    with open(model_name + "_dist_con_log.txt", encoding="utf8", mode="w") as f:
        for strings in outs:
            for string in strings:
                print(string, end="\n", file=f)
            print("\n\n", file=f)

t2 = time.time()
print("Time :", t2 - t1)


