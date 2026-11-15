cmake_minimum_required(VERSION 3.16)

project(numerical-reconstruction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Files extracted from downloaded archives get the extraction time as timestamp (CMake >= 3.24)
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# ════ Eigen ═══════════════════════════════════════════════════════
# Use the installed Eigen (e.g. brew install eigen), otherwise download the headers.
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen3 not found, downloading Eigen 3.4.0")
    include(FetchContent)
    FetchContent_Declare(eigen
        URL      https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
        SOURCE_SUBDIR header-only-do-not-configure  # unpack only, skip add_subdirectory
    )
    FetchContent_MakeAvailable(eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
endif()
# ══════════════════════════════════════════════════════════════════


# ════ Executable ══════════════════════════════════════════════════
add_executable(numerical-reconstruction
    src/config.cpp
    src/io.cpp
    src/main.cpp
)
target_include_directories(numerical-reconstruction PRIVATE include)
target_link_libraries(numerical-reconstruction PRIVATE Eigen3::Eigen)
target_compile_options(numerical-reconstruction PRIVATE -Wall -Wextra -Wpedantic)
# ══════════════════════════════════════════════════════════════════
