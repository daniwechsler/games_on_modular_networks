cmake_minimum_required(VERSION 2.8.9)
project(games_on_graphs)

SET(CMAKE_INSTALL_PREFIX ./)

file(GLOB SOURCES "./src/*.cpp")

find_library(GSL_LIB gsl)
find_library(GSLCBLAS_LIB gslcblas)
find_library(IGRAPH_LIB igraph)

add_executable(games_on_graphs ${SOURCES})

target_link_libraries(games_on_graphs PRIVATE "${GSL_LIB}")
target_link_libraries(games_on_graphs PRIVATE "${GSLCBLAS_LIB}")
target_link_libraries(games_on_graphs PRIVATE "${IGRAPH_LIB}")

install(TARGETS games_on_graphs DESTINATION bin)
