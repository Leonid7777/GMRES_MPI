target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE
  mat_methods.cpp
  vec_methods.cpp
)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "-fopenmp")
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)