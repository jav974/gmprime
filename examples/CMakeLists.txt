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

    # Ways compared on a sample smaller than the data, then on the whole data
    add_test(NAME gmz_sampled COMMAND gmz -S 1024 t ${GMZ_SAMPLE})
    add_test(NAME gmz_sampled_fold COMMAND gmz -m fold -S 1024 t ${GMZ_SAMPLE})
    add_test(NAME gmz_exhaustive COMMAND gmz -S 0 t ${GMZ_SAMPLE})

    # Faster levels, for text and folded data
    add_test(NAME gmz_level1 COMMAND gmz -1 t ${GMZ_SAMPLE})
    add_test(NAME gmz_level2 COMMAND gmz -2 -C on t ${GMZ_SAMPLE})
    add_test(NAME gmz_level1_fold COMMAND gmz -1 -m fold t ${GMZ_SAMPLE})
    add_test(NAME gmz_level2_fold COMMAND gmz -2 -m fold-full t ${GMZ_SAMPLE})

    # Machine code models forced on the sample: binary, x86 in 32 and 64 bits, at the faster level the branch
    # transform alone
    add_test(NAME gmz_binary COMMAND gmz -x bin t ${GMZ_SAMPLE})
    add_test(NAME gmz_binary_fold COMMAND gmz -x bin -m fold t ${GMZ_SAMPLE})
    add_test(NAME gmz_x86 COMMAND gmz -x 32 t ${GMZ_SAMPLE})
    add_test(NAME gmz_x86_64 COMMAND gmz -x 64 t ${GMZ_SAMPLE})
    add_test(NAME gmz_x86_level1 COMMAND gmz -1 -x 32 t ${GMZ_SAMPLE})

    # Through files, compressing then decompressing, with a forced table size
    add_test(NAME gmz_files COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLE}
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/sample -DOPTIONS=-T\;18 -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)

    # Decompressing with another table size than the compression decodes other data: the checksum rejects it
    add_test(NAME gmz_checksum COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLE}
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/mismatch -DOPTIONS=-T\;18 -DDECOMPRESS_OPTIONS=-T\;17
             -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)

    # Records with floating point fields on lattices, and 16-bit images of either byte order, written by gmz_samples
    add_executable(gmz_samples gmz_samples.cpp)
    set(GMZ_SAMPLES ${CMAKE_CURRENT_BINARY_DIR}/samples)
    add_test(NAME gmz_samples_written COMMAND ${CMAKE_COMMAND} -E make_directory ${GMZ_SAMPLES})
    add_test(NAME gmz_samples_made COMMAND gmz_samples ${GMZ_SAMPLES})
    set_tests_properties(gmz_samples_written PROPERTIES FIXTURES_SETUP gmz_samples_directory)
    set_tests_properties(gmz_samples_made PROPERTIES FIXTURES_REQUIRED gmz_samples_directory FIXTURES_SETUP gmz_samples)
    add_test(NAME gmz_records COMMAND gmz t ${GMZ_SAMPLES}/records.bin)
    add_test(NAME gmz_records_level1 COMMAND gmz -1 t ${GMZ_SAMPLES}/records.bin)
    add_test(NAME gmz_image_little_endian COMMAND gmz -S 4096 t ${GMZ_SAMPLES}/image_le.bin)
    add_test(NAME gmz_image_big_endian COMMAND gmz -S 0 t ${GMZ_SAMPLES}/image_be.bin)
    add_test(NAME gmz_records_files COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLES}/records.bin
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/records -DOPTIONS=-T\;18 -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)
    # Decoded with another table size: garbage, which the lattice of IBM values or the checksum rejects
    add_test(NAME gmz_records_checksum COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLES}/records.bin
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/records_mismatch -DOPTIONS=-T\;18 -DDECOMPRESS_OPTIONS=-T\;17
             -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)
    # Numbers in text: a time series, JSON, a log, and a text where no byte is free to mark numbers
    add_test(NAME gmz_numbers_csv COMMAND gmz t ${GMZ_SAMPLES}/numbers.csv)
    add_test(NAME gmz_numbers_json COMMAND gmz -1 t ${GMZ_SAMPLES}/numbers.json)
    add_test(NAME gmz_numbers_log COMMAND gmz -S 2048 t ${GMZ_SAMPLES}/numbers.log)
    add_test(NAME gmz_numbers_only COMMAND gmz -2 -N on t ${GMZ_SAMPLES}/numbers.log)
    add_test(NAME gmz_numbers_off COMMAND gmz -N off t ${GMZ_SAMPLES}/numbers.csv)
    add_test(NAME gmz_numbers_crowded COMMAND gmz -N on t ${GMZ_SAMPLES}/crowded.txt)
    add_test(NAME gmz_numbers_one_line COMMAND gmz t ${GMZ_SAMPLES}/oneline.json)
    add_test(NAME gmz_numbers_files COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLES}/numbers.csv
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/numbers -DOPTIONS=-T\;18 -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)
    add_test(NAME gmz_numbers_checksum COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLES}/numbers.csv
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/numbers_mismatch -DOPTIONS=-T\;18 -DDECOMPRESS_OPTIONS=-T\;17
             -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)
    # Decimal doubles in shapes of varying length, and floats in records behind a header
    add_test(NAME gmz_scan COMMAND gmz t ${GMZ_SAMPLES}/polylines.bin)
    add_test(NAME gmz_scan_level1 COMMAND gmz -1 -S 4096 t ${GMZ_SAMPLES}/polylines.bin)
    add_test(NAME gmz_scan_files COMMAND ${CMAKE_COMMAND} -DGMZ=$<TARGET_FILE:gmz> -DINPUT=${GMZ_SAMPLES}/polylines.bin
             -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/polylines -DOPTIONS=-T\;18 -P ${CMAKE_CURRENT_SOURCE_DIR}/round_trip.cmake)
    add_test(NAME gmz_header_records COMMAND gmz t ${GMZ_SAMPLES}/header_records.bin)
    set_tests_properties(gmz_scan gmz_scan_level1 gmz_scan_files gmz_header_records PROPERTIES FIXTURES_REQUIRED gmz_samples)
    set_tests_properties(gmz_records gmz_records_level1 gmz_image_little_endian gmz_image_big_endian gmz_records_files
                         gmz_records_checksum gmz_numbers_checksum gmz_numbers_csv gmz_numbers_json gmz_numbers_log gmz_numbers_only gmz_numbers_off
                         gmz_numbers_crowded gmz_numbers_one_line gmz_numbers_files PROPERTIES FIXTURES_REQUIRED gmz_samples)

    # Invalid command lines
    add_test(NAME gmz_unknown_option COMMAND gmz -q t ${GMZ_SAMPLE})
    add_test(NAME gmz_missing_command COMMAND gmz -m text)
    add_test(NAME gmz_unknown_command COMMAND gmz x ${GMZ_SAMPLE})
    add_test(NAME gmz_unknown_code_mode COMMAND gmz -x 16 t ${GMZ_SAMPLE})
    add_test(NAME gmz_unknown_numbers_mode COMMAND gmz -N maybe t ${GMZ_SAMPLE})
    set_tests_properties(gmz_unknown_option gmz_missing_command gmz_unknown_command gmz_unknown_code_mode gmz_unknown_numbers_mode
                         PROPERTIES WILL_FAIL TRUE)
endif ()

if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
    # A slow test fails instead of holding up the run. ASan's stack-use-after-return detection allocates a fake
    # stack frame on most calls, which costs these runs about 200 times more with GCC on ARM64
    get_property(gmz_tests DIRECTORY PROPERTY TESTS)
    set_tests_properties(${gmz_tests} PROPERTIES TIMEOUT 60 ENVIRONMENT "ASAN_OPTIONS=detect_stack_use_after_return=0")
    # The structure plans compare many ways on the samples: under sanitizers and emulation (Rosetta, QEMU), minutes
    get_property(sample_tests DIRECTORY PROPERTY TESTS)
    list(FILTER sample_tests INCLUDE REGEX "gmz_(records|image|numbers|scan|header)")
    set_tests_properties(${sample_tests} PROPERTIES TIMEOUT 180)
endif ()
