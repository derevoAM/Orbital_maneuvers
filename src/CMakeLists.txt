cmake_minimum_required(VERSION 3.24)
project(Orbital_maneuvers CXX)

set(CMAKE_CXX_STANDARD 20)

file(GLOB files "*.h")

add_library(src INTERFACE ${files})