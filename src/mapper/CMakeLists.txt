add_library(mapper mapper.cpp)
install (TARGETS mapper DESTINATION bin)
target_link_libraries(mapper Sweeper winnowing statistics thread_pool pthread)
install(FILES mapper.hpp DESTINATION include)