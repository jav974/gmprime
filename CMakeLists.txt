cmake_minimum_required(VERSION 3.28)
project(gmprime)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find the GMP library
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)
find_package(GMP REQUIRED)

option(USE_BOOST_UT "Use pre-installed boost.ut package" OFF)
if (USE_BOOST_UT)
    find_package(ut REQUIRED)
else ()
    include(FetchContent)
    FetchContent_Declare(
            ut
            GIT_REPOSITORY https://github.com/boost-ext/ut.git
            GIT_TAG        v2.1.0
    )
    FetchContent_MakeAvailable(ut)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
# Useful troubleshooting flags in release mode
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fsanitize=address -fno-omit-frame-pointer")

# Build gmprime static lib
add_library(gmprime STATIC
        src/gmprime.cpp
        src/gmprime.hpp
        src/gmprime.tpp
        src/gmprime_ios.hpp
        src/gmprime_reference.tpp
        src/gmprime_stream.tpp
        src/gmprime_adapter.tpp
        src/gmprime_literals.hpp)
target_link_libraries(gmprime ${GMP_LIBRARIES})

# Create 'tests' executable, linked with gmprime to run tests suites
add_executable(tests tests/main.cpp)
target_link_libraries(tests PRIVATE gmprime Boost::ut)
target_include_directories(tests PRIVATE src)

# Enable testing
enable_testing()

# Run 'tests'
add_test(NAME tests
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/tests
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})