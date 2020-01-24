cmake_minimum_required(VERSION 3.5)
project(qfib)

find_package(OpenMP)

file(GLOB_RECURSE source_files sources/*.cpp)
add_executable(qfib ${source_files})


set_property(TARGET qfib PROPERTY CXX_STANDARD 11)
#target_compile_options(compressor PUBLIC -frounding-math -fsignaling-nans)

target_include_directories(qfib
    PUBLIC eigen3)
target_link_libraries(qfib
        PUBLIC OpenMP::OpenMP_CXX
        )
