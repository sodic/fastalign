add_library(types types.cpp)
install(TARGETS types DESTINATION bin)
install(FILES types.h DESTINATION include)