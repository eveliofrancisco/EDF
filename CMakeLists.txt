cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(EDF VERSION 1.0.0 LANGUAGES Fortran)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

enable_language(Fortran)

set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

find_package(MKL REQUIRED)

add_executable(edf src/edf.f)
target_link_libraries(edf MKL::MKL)
add_subdirectory(src)
target_include_directories(edf PRIVATE ${CMAKE_BINARY_DIR}/modules)