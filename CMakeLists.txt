cmake_minimum_required(VERSION 3.16)

project(numerical-reconstruction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(numerical-reconstruction
    src/config.cpp
    src/main.cpp
)

target_include_directories(numerical-reconstruction PRIVATE
    include
)
