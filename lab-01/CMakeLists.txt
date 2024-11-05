cmake_minimum_required(VERSION 3.12.0)
project(lab02 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-01 src/lab-01.cpp src/Poisson1D.cpp)
deal_ii_setup_target(lab-01)
