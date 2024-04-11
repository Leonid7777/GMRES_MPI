target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE
  gmres.cpp
  mat_methods.cpp
  vec_methods.cpp
)

target_link_libraries(${CMAKE_PROJECT_NAME})