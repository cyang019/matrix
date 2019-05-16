cmake_minimum_required (VERSION 3.8)
project(matrix)

add_library(matrix SHARED)
include("include/CMakeLists.txt")
include("src/CMakeLists.txt")

