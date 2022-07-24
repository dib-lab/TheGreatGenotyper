cmake_minimum_required(VERSION 3.1.3)

project(caches)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CACHE_INCLUDE_DIRS ${caches_SOURCE_DIR}/include CACHE INTERNAL "")

option(BUILD_TEST "Build tests for the project" ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -g -fsanitize=thread")

if (${BUILD_TEST})
  enable_testing()
endif()

add_subdirectory(deps)
add_subdirectory(src)

if (${BUILD_TEST})
  add_subdirectory(test)
endif()
