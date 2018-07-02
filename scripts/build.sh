#! /bin/bash
cmake -H. -Bbuild; # Creates CMake configuration files inside the folder "build"
cmake --build build -- -j3; # generates output programs in the folder "bin"