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

# Create 'tests' executable, linked with gmprime to run tests suites
add_executable(tests main.cpp)
target_link_libraries(tests PRIVATE gmprime Boost::ut)
target_include_directories(tests PRIVATE src)

# Enable testing
enable_testing()

# Run 'tests'
add_test(NAME tests
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/tests
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

include(CTest)