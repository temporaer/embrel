# - Find dsdp
# Find the native DSDP headers and libraries.
#
#  DSDP_INCLUDE_DIRS - where to find dsdp/dsdp.h, etc.
#  DSDP_LIBRARIES    - List of libraries when using dsdp.
#  DSDP_FOUND        - True if dsdp found.

# Look for the header file.
FIND_PATH(DSDP_INCLUDE_DIR NAMES dsdp/dsdp5.h)
MARK_AS_ADVANCED(DSDP_INCLUDE_DIR)

# Look for the library.
FIND_LIBRARY(DSDP_LIBRARY NAMES dsdp)
MARK_AS_ADVANCED(DSDP_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set DSDP_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DSDP DEFAULT_MSG DSDP_LIBRARY DSDP_INCLUDE_DIR)

IF(DSDP_FOUND)
  SET(DSDP_LIBRARIES ${DSDP_LIBRARY})
  SET(DSDP_INCLUDE_DIRS ${DSDP_INCLUDE_DIR})
ELSE(DSDP_FOUND)
  SET(DSDP_LIBRARIES)
  SET(DSDP_INCLUDE_DIRS)
ENDIF(DSDP_FOUND)

