cmake_minimum_required(VERSION 3.12.0)
project(lab05 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-05 src/lab-05.cpp src/Poisson3D.cpp)
deal_ii_setup_target(lab-05)

# add_executable(lab-05_parallel src/lab-05_parallel.cpp src/Poisson3D_parallel.cpp)
# deal_ii_setup_target(lab-05_parallel)
