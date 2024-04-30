target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE
  gmres.cpp
  mat_methods.cpp
  vec_methods.cpp
)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "-fopenmp")
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

target_link_libraries(${CMAKE_PROJECT_NAME})