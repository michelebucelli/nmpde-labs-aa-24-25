cmake_minimum_required(VERSION 3.12.0)
project(lab08 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(lab-08-exercise1 src/lab-08-exercise1.cpp src/Heat.cpp)
deal_ii_setup_target(lab-08-exercise1)

add_executable(lab-08-exercise2 src/lab-08-exercise2.cpp src/HeatNonLinear.cpp)
deal_ii_setup_target(lab-08-exercise2)
