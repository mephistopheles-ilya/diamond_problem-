cmake_minimum_required(VERSION 3.11)

project (diamond_problem VERSION 1.0.0)

file(GLOB SRCS "source/*.cpp")


# static library
add_library(libgeom STATIC ${SRCS})

target_compile_options(libgeom PUBLIC 
    $<$<CONFIG:Debug>: -O0 -g -Wall -Wextra -Wpedantic -fsanitize=leak,undefined,address>
    $<$<CONFIG:Release>: -O3 -DNDEBUG -Wall -Wextra -Wpedantic>
    -std=c++20 -DBOOST_ALLOW_DEPRECATED_HEADERS 
)

target_link_options(libgeom PUBLIC
    $<$<CONFIG:Debug>: -fsanitize=leak,undefined,address>
)

target_link_libraries(libgeom PUBLIC
    boost_program_options
    gmp
)

# executable to get contours
add_executable(get.out get_contour.cpp)
target_link_libraries(get.out libgeom)
set_target_properties(get.out PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


# executable to transform contours
add_executable(tran.out transform_contour.cpp)
target_link_libraries(tran.out libgeom)
set_target_properties(tran.out PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


# executable to plot 3d model
add_executable(splot.out splot.cpp source/point3d.cpp)
target_link_libraries(splot.out PRIVATE
    boost_program_options
    gmp
)
target_compile_options(splot.out PRIVATE -O3 -DNDEBUG -DBOOST_ALLOW_DEPRECATED_HEADERS)
set_target_properties(splot.out PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# executable for gorizontal comparing of polylines
add_executable(gorizontal_compare.out gorizontal_compare.cpp source/point2d.cpp)
target_compile_options(gorizontal_compare.out PRIVATE -O3 -DNDEBUG -DBOOST_ALLOW_DEPRECATED_HEADERS)
set_target_properties(gorizontal_compare.out PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# executable for symmetric difference comparing of polylines
add_executable(sym_diff.out sym_diff.cpp source/point2d.cpp)
target_include_directories(sym_diff.out PRIVATE
    ~/libs/boost_1_86_0
)
target_compile_options(sym_diff.out PRIVATE -O3 -DNDEBUG -DBOOST_ALLOW_DEPRECATED_HEADERS)
set_target_properties(sym_diff.out PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})



