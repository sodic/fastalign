add_library(FASTARead FASTARead.cpp)
install (TARGETS FASTARead DESTINATION bin)
install (FILES FASTARead.hpp DESTINATION include)