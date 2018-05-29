add_library(Sweeper Sweeper.cpp)
install(TARGETS Sweeper DESTINATION bin)
target_link_libraries(Sweeper config)
install(FILES Sweeper.h DESTINATION include)