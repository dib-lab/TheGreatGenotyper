cmake_minimum_required(VERSION 3.1)
project(TheGreatGenotyper)

set (CMAKE_CXX_STANDARD 17)

enable_testing()

add_executable(computeHistogram computeHistogram.cpp histogram.cpp)
target_link_libraries(computeHistogram -lKMC )
#target_link_directories(${CMAKE_BINARY_DIR}/../KMC/)
