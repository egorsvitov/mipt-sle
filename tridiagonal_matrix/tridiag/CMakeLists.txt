cmake_minimum_required(VERSION 3.10)

project(tridiag)

add_library(tridiag STATIC src/tridiag.cpp)
target_include_directories(tridiag PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)