add_library(Sweeper Sweeper.cpp)
install(TARGETS Sweeper DESTINATION bin)
install(FILES Sweeper.h DESTINATION include)