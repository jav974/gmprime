cmake_minimum_required(VERSION 3.5)
project(gmprime)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find the GMP library
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)
find_package(GMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
# Useful troubleshooting flags in release mode
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fsanitize=address -fno-omit-frame-pointer")

# Code Coverage Configuration
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Add required flags (GCC & LLVM/Clang)
    target_compile_options(coverage_config INTERFACE
            -O0        # no optimization
            -g         # generate debug info
            --coverage # sets all required flags
    )
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
        target_link_options(coverage_config INTERFACE --coverage)
    else()
        target_link_libraries(coverage_config INTERFACE --coverage)
    endif()
endif()

set(GMPRIME_PUBLIC_INCLUDE_DIRECTORIES  include)
set(GMPRIME_PRIVATE_INCLUDE_DIRECTORIES src)

file(GLOB_RECURSE GMPRIME_SOURCE_FILES
        "${GMPRIME_PRIVATE_INCLUDE_DIRECTORIES}/*.cpp"
        "${GMPRIME_PRIVATE_INCLUDE_DIRECTORIES}/*.hpp*"
        "${GMPRIME_PRIVATE_INCLUDE_DIRECTORIES}/*.tpp*"
)
file(GLOB_RECURSE GMPRIME_INCLUDE_FILES
        "${GMPRIME_PUBLIC_INCLUDE_DIRECTORIES}/*.hpp*"
)

# Build gmprime static lib
add_library(${PROJECT_NAME} STATIC ${GMPRIME_INCLUDE_FILES} ${GMPRIME_SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE ${GMP_LIBRARIES} PUBLIC coverage_config)
target_include_directories(${PROJECT_NAME} PUBLIC ${GMPRIME_PUBLIC_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GMPRIME_PRIVATE_INCLUDE_DIRECTORIES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${GMPRIME_INCLUDE_FILES}")

# Install rules
install(TARGETS ${PROJECT_NAME} coverage_config
        EXPORT gmprimeConfig
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

export(EXPORT gmprimeConfig
        FILE "${CMAKE_CURRENT_BINARY_DIR}/gmprimeConfig.cmake"
        NAMESPACE gmprime::)

option (BUILD_TESTING "Build the testing tree." ON)
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
    enable_testing()
    add_subdirectory(tests)
endif()
