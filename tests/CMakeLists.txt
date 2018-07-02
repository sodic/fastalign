include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

include_directories("${PROJECT_SOURCE_DIR}/src/winnowing")
add_executable(testWinnowing winnowing_test.cpp)
target_link_libraries(testWinnowing gtest_main)
target_link_libraries(testWinnowing winnowing)

#include_directories("${PROJECT_SOURCE_DIR}/src/mapper")
#add_executable(testMapper mapper_test.cpp)
#target_link_libraries(testMapper gtest_main)
#target_link_libraries(testMapper mapper)

include_directories("${PROJECT_SOURCE_DIR}/src/sweeper")
add_executable(testSweeper sweeper_test.cpp)
target_link_libraries(testSweeper gtest_main)

include_directories("${PROJECT_SOURCE_DIR}/src/statistics")
add_executable(testStats statistics_test.cpp)
target_link_libraries(testStats statistics gtest_main)
