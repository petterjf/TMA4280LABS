cmake_minimum_required(VERSION 2.8)

enable_language(CXX)
enable_language(Fortran)

add_executable(poisson_pure poisson_pure.cpp)
add_executable(poisson_hybrid poisson_hybrid.cpp)

find_package(OpenMP REQUIRED)
set_source_files_properties(poisson_hybrid.cpp PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
include_directories(${MPI_CXX_INCLUDE_PATH})

add_library(common fst.f)
target_link_libraries(poisson_pure common ${MPI_CXX_LIBRARIES})
target_link_libraries(poisson_hybrid common ${MPI_CXX_LIBRARIES})

set_target_properties(poisson_pure PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(poisson_hybrid PROPERTIES LINKER_LANGUAGE CXX LINK_FLAGS ${OpenMP_CXX_FLAGS})