# Add the source files
set(UTILS_SOURCES
  math.h
  matrix.cpp
  matrix.h
  strings.cpp
  strings.h
  computation_grid.cpp
  computation_grid.h
)

# Add the library
add_library(utils ${UTILS_SOURCES})
