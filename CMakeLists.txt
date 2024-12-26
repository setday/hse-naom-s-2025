cmake_minimum_required(VERSION 3.20)
project(hse_naom_s2025)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall -Wextra -Wshadow -O0 -g) # Debug
# add_compile_options(-Ofast) # Release

#find_package(GSL REQUIRED)
#link_libraries(GSL::gsl)

#set(BLA_VENDER OpenBLAS)
#find_package(BLAS REQUIRED)
#link_libraries(${BLAS_LIBRARIES})

add_executable(hse_naom_s2025 lab04/test.cpp)