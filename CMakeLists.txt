cmake_minimum_required(VERSION 3.16)

project(numerical-reconstruction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(numerical-reconstruction
    src/config.cpp
    src/geometry.cpp
    src/scanner.cpp
    src/projector.cpp
    src/main.cpp
)

target_include_directories(numerical-reconstruction PRIVATE
    include
)

target_link_libraries(numerical-reconstruction PRIVATE Eigen3::Eigen)
