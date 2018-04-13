add_library(mapper mapper.cpp)
install (TARGETS mapper DESTINATION bin)
install (FILES mapper.h DESTINATION include)