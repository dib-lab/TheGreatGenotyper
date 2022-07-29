cmake_minimum_required(VERSION 3.1)
project(TheGreatGenotyper)

set (CMAKE_CXX_STANDARD 17)




add_executable(computeHistogram computeHistogram.cpp histogram.cpp)
target_link_libraries(computeHistogram -lKMC )


file(COPY project_config.yaml  sample_table.csv subsample_table.csv Snakefile env.yaml config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
