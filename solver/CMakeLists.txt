cmake_minimum_required(VERSION 3.0)
project(solver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")


add_executable(solver main.cpp)

add_subdirectory(properties)
add_subdirectory(utils)
add_subdirectory(semi_implicit)

# MPI
find_package(MPI REQUIRED)

target_link_libraries(solver PUBLIC MPI::MPI_CXX properties utils semi_implicit)
