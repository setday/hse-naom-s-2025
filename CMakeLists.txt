cmake_minimum_required(VERSION 3.20)
project(HSE_NaOM_S2024)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -Wextra -Wshadow -O0 -g) # Debug
# add_compile_options(-Ofast) # Release

find_package(GSL REQUIRED)
link_libraries(GSL::gsl)

set(BLA_VENDER OpenBLAS)
find_package(BLAS REQUIRED)
link_libraries(${BLAS_LIBRARIES})

add_executable(HSE_NaOM_S2024 lab01/main.cpp)