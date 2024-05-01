target_sources(${CMAKE_PROJECT_NAME}
  PRIVATE
  gmres.cpp
  mat_methods.cpp
  vec_methods.cpp
  ${PROJECT_SOURCE_DIR}/class/world_class.cpp
)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "-fopenmp")
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)

target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE
  ${PROJECT_SOURCE_DIR}/class
)