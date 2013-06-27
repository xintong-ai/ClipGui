project(GreinerClip)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
link_directories(/usr/lib)

find_library(SX_LIB sx /usr/lib)

find_library(M_LIB m)

#SET(CMAKE_CXX_FLAGS "-std=c++0x")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} ${SX_LIB} ${M_LIB})
