#!/bin/bash

files=(Reflect_00001 Reflect_10001 Reflect_20001 Reflect_30002 Reflect_40001 Reflect_50001 Reflect_60001 Reflect_70001)

for file in "${files[@]}"
do
    txtfile="${file}.txt"
    plyfile1="${file}_fcs5G.ply"
    plyfile2="${file}_fps5G.ply"
    txtfile1="${file}_fcs5G.txt"
    txtfile2="${file}_fps5G.txt"
    objfile1="${file}_fcs5G.obj"
    objfile2="${file}_fps5G.obj"

    ./obj_to_txt.out --file_ply_init=ChangedModels/${plyfile1} --convert=ply --file_txt=ChangedModels/${txtfile1}
    ./ply_to_obj.out --file_ply=ChangedModels/${plyfile1} --file_obj=ChangedModels/${objfile1}
    
    rm -rf init_examples/Contour*
    ./get.out --projections=400 --init_file=ChangedModels/${txtfile1} --directory=init_examples 
    rm -rf res_examples/Contour*
    ./tran.out --projections=400 --grid=1 --gr_method=ceils --shift=0 --dist_points=0.008 --directory_in=init_examples \
        --directory_out=res_examples
    dir="${file}_fcs5G_grid"
    mkdir -p ${dir}
    cp res_examples/Contour* ${dir}/

    ./obj_to_txt.out --file_ply_init=ChangedModels/${plyfile2} --convert=ply --file_txt=ChangedModels/${txtfile2}
    ./ply_to_obj.out --file_ply=ChangedModels/${plyfile2} --file_obj=ChangedModels/${objfile2}

    rm -rf init_examples/Contour*
    ./get.out --projections=400 --init_file=ChangedModels/${txtfile2} --directory=init_examples 
    rm -rf res_examples/Contour*
    ./tran.out --projections=400 --grid=1 --gr_method=ceils --shift=0 --dist_points=0.008 --directory_in=init_examples \
        --directory_out=res_examples
    dir="${file}_fps5G_grid"
    mkdir -p ${dir}
    cp res_examples/Contour* ${dir}/

done



