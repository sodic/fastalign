add_library(winnowing winnowing.cpp)
install(TARGETS winnowing DESTINATION bin)
install(FILES winnowing.hpp DESTINATION include)