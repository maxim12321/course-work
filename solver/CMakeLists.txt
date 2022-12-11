cmake_minimum_required(VERSION 3.14)
project(solver)

set(CMAKE_CXX_STANDARD 20)

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
        utils/strings.cpp
        utils/vector.cpp
)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(solver PUBLIC MPI::MPI_CXX)
