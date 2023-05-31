cmake_minimum_required(VERSION 3.14)
project(solver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_BUILD_TYPE "Release")

add_executable(
        solver

        main.cpp

        material.cpp
        properties_manager.cpp
        properties_wrapper.cpp
        table_value.cpp

        explicit/explicit_solver.cpp

        semi_implicit/column_iteration_solver.cpp
        semi_implicit/row_iteration_solver.cpp
        semi_implicit/semi_implicit_solver.cpp
        semi_implicit/tridiagonal.cpp

        utils/matrix.cpp
        utils/result_saver.cpp
        utils/strings.cpp
        utils/vector.cpp
)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(solver PUBLIC MPI::MPI_CXX)

find_package(OpenMP REQUIRED)
target_link_libraries(solver PUBLIC OpenMP::OpenMP_CXX)
