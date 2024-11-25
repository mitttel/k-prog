cmake_minimum_required(VERSION 3.29)
project(k_prog)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -finput-charset=UTF-8")

add_executable(k_prog main.cpp
        matrix.cpp
        matrix.h)
