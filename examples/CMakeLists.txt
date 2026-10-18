# Example programs built on gmprime
add_executable(gmz gmz.cpp)
target_link_libraries(gmz PRIVATE gmprime)

if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
    set(GMZ_SAMPLE ${PROJECT_SOURCE_DIR}/CMakeLists.txt)

    # Round trips in memory, in each mode
    add_test(NAME gmz_text COMMAND gmz -C off t ${GMZ_SAMPLE})
    add_test(NAME gmz_capitals COMMAND gmz -C on t ${GMZ_SAMPLE})
    add_test(NAME gmz_auto COMMAND gmz t ${GMZ_SAMPLE})
    add_test(NAME gmz_fold COMMAND gmz -m fold t ${GMZ_SAMPLE})
    add_test(NAME gmz_fold_full COMMAND gmz -m fold-full t ${GMZ_SAMPLE})
    add_test(NAME gmz_bitmap COMMAND gmz -m bitmap t ${GMZ_SAMPLE})

    # Through files, compressing then decompressing, with a forced table size
    add_test(NAME gmz_files COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLE}
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/sample -DOPTIONS=-T\;18 -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)

    # Invalid command lines
    add_test(NAME gmz_unknown_option COMMAND gmz -q t ${GMZ_SAMPLE})
    add_test(NAME gmz_missing_command COMMAND gmz -m text)
    add_test(NAME gmz_unknown_command COMMAND gmz x ${GMZ_SAMPLE})
    set_tests_properties(gmz_unknown_option gmz_missing_command gmz_unknown_command PROPERTIES WILL_FAIL TRUE)
endif ()

if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
    # A slow test fails instead of holding up the run. ASan's stack-use-after-return detection allocates a fake
    # stack frame on most calls, which costs these runs about 200 times more with GCC on ARM64
    get_property(gmz_tests DIRECTORY PROPERTY TESTS)
    set_tests_properties(${gmz_tests} PROPERTIES TIMEOUT 60 ENVIRONMENT "ASAN_OPTIONS=detect_stack_use_after_return=0")
endif ()
