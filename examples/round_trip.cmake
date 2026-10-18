# Compresses INPUT with GMZ into OUTPUT.gmz, decompresses it into OUTPUT, and checks OUTPUT equals INPUT.
# OPTIONS: gmz options, as a CMake list, given to both commands
execute_process(COMMAND ${GMZ} ${OPTIONS} c ${INPUT} ${OUTPUT}.gmz RESULT_VARIABLE compress_result)
if (NOT compress_result EQUAL 0)
    message(FATAL_ERROR "gmz c failed: ${compress_result}")
endif ()
execute_process(COMMAND ${GMZ} ${OPTIONS} d ${OUTPUT}.gmz ${OUTPUT} RESULT_VARIABLE decompress_result)
if (NOT decompress_result EQUAL 0)
    message(FATAL_ERROR "gmz d failed: ${decompress_result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${INPUT} ${OUTPUT} RESULT_VARIABLE compare_result)
if (NOT compare_result EQUAL 0)
    message(FATAL_ERROR "decompressed file differs from ${INPUT}")
endif ()
