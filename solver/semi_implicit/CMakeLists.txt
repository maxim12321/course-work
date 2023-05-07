set(SEMI_IMPLICIT_SOURCES
    tridiagonal.cpp
    tridiagonal.h
    semi_implicit_solver.cpp
    semi_implicit_solver.h
)

add_library(semi_implicit ${SEMI_IMPLICIT_SOURCES})

target_link_libraries(semi_implicit properties)