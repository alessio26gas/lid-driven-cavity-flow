cmake_minimum_required(VERSION 3.10)
project(CFD VERSION 1.0 LANGUAGES C)
add_executable(cfd_solver src/main.c)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -flto -funroll-loops -ffast-math" CACHE STRING "Optimization flags for Release" FORCE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)