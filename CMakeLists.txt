cmake_minimum_required(VERSION 3.16)

project(numerical-reconstruction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ════ Eigen ═══════════════════════════════════════════════════════
find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(FATAL_ERROR
        "Eigen3 not found. Install it (macOS: brew install eigen, "
        "Ubuntu: sudo apt install libeigen3-dev) or pass "
        "-DEigen3_DIR=<directory containing Eigen3Config.cmake>.")
endif()
# ══════════════════════════════════════════════════════════════════

# ════ Executable ══════════════════════════════════════════════════
add_executable(numerical-reconstruction
    src/tools.cpp
    src/config.cpp
    src/geometry.cpp
    src/scanner.cpp
    src/projector.cpp
    src/main.cpp
    src/io.cpp
)

target_include_directories(numerical-reconstruction PRIVATE
    include
)

target_link_libraries(numerical-reconstruction PRIVATE Eigen3::Eigen)


