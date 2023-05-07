set(PROPERTIES_SOURCES
  material.cpp
  material.h
  properties_manager.cpp
  properties_manager.h
  properties_wrapper.cpp
  properties_wrapper.h
  table_value.cpp
  table_value.h
)

add_library(properties ${PROPERTIES_SOURCES})

target_link_libraries(properties utils)
