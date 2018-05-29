add_library(statistics statistics.cpp)
install(TARGETS statistics DESTINATION bin)
target_link_libraries(statistics config)
install(FILES statistics.hpp DESTINATION include)