add_executable(test
  algorithm_gmres.cpp
  create_mat_r.cpp
  create_mat_r.h
  ${PROJECT_SOURCE_DIR}/class/world_class.cpp

)

target_link_libraries(test
${CMAKE_PROJECT_NAME}
)

target_include_directories(test PRIVATE
  ${PROJECT_SOURCE_DIR}/class
)