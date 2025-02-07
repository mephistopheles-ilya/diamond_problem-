#!/bin/bash

files=(Reflect_00001 Reflect_10001 Reflect_20001 Reflect_30002 Reflect_40001 Reflect_50001 Reflect_60001 Reflect_70001)

for file in "${files[@]}"
do
    txtfile="${file}.txt"
    rm -rf init_examples/Contour*
    ./get.out --projections=400 --init_file=InitialModels/${txtfile} --directory=init_examples 
    rm -rf res_examples/Contour*
    ./tran.out --projections=400 --grid=1 --gr_method=ceils --shift=0 --dist_points=0.008 --directory_in=init_examples \
        --directory_out=res_examples
    dir="${file}_grid"
    mkdir -p ${dir}
    cp res_examples/Contour* ${dir}/
done


