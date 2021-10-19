cmake_minimum_required(VERSION 3.13)

project (Shining-Pom VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(shiningpom "")

# Source Code
add_subdirectory(src)

#AMReX
find_package(AMReX REQUIRED 2D HINTS /prod/AMReX/GCC-9_2D)
find_package(RapidJSON REQUIRED)
find_package(ZLIB REQUIRED)
find_package(HDF5 1.12.0 REQUIRED)
target_link_libraries(shiningpom AMReX::amrex AMReX::Flags_CXX)
target_link_libraries(shiningpom ZLIB::ZLIB)
target_link_libraries(shiningpom stdc++fs)

target_include_directories(shiningpom PRIVATE ${RapidJSON_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-static-libstdc++")
set(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g -fsanitize=address")
#set(CMAKE_CXX_FLAGS_RELEASE  "-O2")
