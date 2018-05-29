add_library(config config.cpp)
install(TARGETS config DESTINATION bin)
install(FILES config.hpp DESTINATION include)