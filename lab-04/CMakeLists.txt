cmake_minimum_required(VERSION 3.12.0)
project(lab04 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-04 src/lab-04.cpp src/Poisson2D.cpp)
deal_ii_setup_target(lab-04)
