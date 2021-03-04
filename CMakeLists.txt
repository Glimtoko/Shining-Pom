cmake_minimum_required(VERSION 3.13)

project (Shining-Pom VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(shiningpom "")

# Source Code
add_subdirectory(src)

#AMReX
find_package(AMReX REQUIRED 2D HINTS /prod/AMReX2D_CMAKE)
find_package(RapidJSON REQUIRED)
find_package(ZLIB REQUIRED)
target_link_libraries(shiningpom AMReX::amrex AMReX::Flags_CXX)
target_link_libraries(shiningpom ZLIB::ZLIB)

target_include_directories(shiningpom PRIVATE ${RapidJSON_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g -fsanitize=address")