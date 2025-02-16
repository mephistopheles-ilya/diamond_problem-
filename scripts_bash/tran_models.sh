#!/bin/bash

files=(Reflect_00001 Reflect_10001 Reflect_20001 Reflect_30002 Reflect_40001 Reflect_50001 Reflect_60001 Reflect_70001)

for file in "${files[@]}"
do
    txtfile="${file}.txt"
    plyfile="${file}.ply"
    plyfile1="${file}_fcs5G.ply"
    plyfile2="${file}_fps5G.ply"
    ./splot.out --initial_model=InitialModels/${txtfile} --output=ChangedModels/${plyfile}
    ./splot.out --initial_model=InitialModels/${txtfile} --output=ChangedModels/${plyfile1}\
        --is_change_slope=1 --parametr_of_slope=0.087 --distribution=normal --sigma=0.0005 --sign_of_c=1
    ./splot.out --initial_model=InitialModels/${txtfile} --output=ChangedModels/${plyfile2}\
        --is_change_slope=1 --parametr_of_slope=0.087 --distribution=normal --sigma=0.0005 --sign_of_c=-1
done


