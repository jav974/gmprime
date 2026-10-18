# Try to find the GNU Multiple Precision Arithmetic Library (GMP) and its C++ bindings (gmpxx)
# See http://gmplib.org/
#
# Defines:
#   GMP_INCLUDES, GMP_LIBRARIES (gmp + gmpxx)
#   Imported targets GMP::gmp and GMP::gmpxx

if (GMP_INCLUDES AND GMP_LIBRARY AND GMPXX_LIBRARY)
    set(GMP_FIND_QUIETLY TRUE)
endif ()

find_path(GMP_INCLUDES
        NAMES
        gmp.h
        gmpxx.h
        PATHS
        $ENV{GMPDIR}
        ${INCLUDE_INSTALL_DIR}
)

find_library(GMP_LIBRARY gmp PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})
find_library(GMPXX_LIBRARY gmpxx PATHS $ENV{GMPDIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
        GMP_INCLUDES GMP_LIBRARY GMPXX_LIBRARY)
mark_as_advanced(GMP_INCLUDES GMP_LIBRARY GMPXX_LIBRARY)

if (GMP_FOUND)
    set(GMP_LIBRARIES ${GMPXX_LIBRARY} ${GMP_LIBRARY})

    if (NOT TARGET GMP::gmp)
        add_library(GMP::gmp UNKNOWN IMPORTED)
        set_target_properties(GMP::gmp PROPERTIES
                IMPORTED_LOCATION "${GMP_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDES}")
    endif ()

    if (NOT TARGET GMP::gmpxx)
        add_library(GMP::gmpxx UNKNOWN IMPORTED)
        set_target_properties(GMP::gmpxx PROPERTIES
                IMPORTED_LOCATION "${GMPXX_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDES}"
                INTERFACE_LINK_LIBRARIES GMP::gmp)
    endif ()
endif ()
