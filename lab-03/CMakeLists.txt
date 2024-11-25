cmake_minimum_required(VERSION 3.12.0)
project(lab03 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-03 src/lab-03.cpp src/Poisson2D.cpp)
deal_ii_setup_target(lab-03)
