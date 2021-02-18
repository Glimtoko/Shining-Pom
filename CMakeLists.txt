cmake_minimum_required(VERSION 3.13)

project (Shining-Pom VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(shiningpom "")

# Source Code
add_subdirectory(src)

#AMReX
find_package(AMReX REQUIRED 2D)
find_package(ZLIB REQUIRED)
target_link_libraries(shiningpom AMReX::amrex AMReX::Flags_CXX)
target_link_libraries(shiningpom ZLIB::ZLIB)