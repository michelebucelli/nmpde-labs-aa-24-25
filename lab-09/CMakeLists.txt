cmake_minimum_required(VERSION 3.12.0)
project(lab09 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-09 src/lab-09.cpp src/Stokes.cpp)
deal_ii_setup_target(lab-09)