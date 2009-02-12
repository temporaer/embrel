# - Find csdp
# Find the native CSDP headers and libraries.
#
#  CSDP_INCLUDE_DIRS - where to find csdp/csdp.h, etc.
#  CSDP_LIBRARIES    - List of libraries when using csdp.
#  CSDP_FOUND        - True if csdp found.

# Look for the header file.
FIND_PATH(CSDP_INCLUDE_DIR NAMES declarations.h )
MARK_AS_ADVANCED(CSDP_INCLUDE_DIR)

# Look for the library.
FIND_LIBRARY(CSDP_LIBRARY NAMES sdp)
MARK_AS_ADVANCED(CSDP_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set CSDP_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CSDP DEFAULT_MSG CSDP_LIBRARY CSDP_INCLUDE_DIR)

IF(CSDP_FOUND)
  SET(CSDP_LIBRARIES ${CSDP_LIBRARY})
  SET(CSDP_INCLUDE_DIRS ${CSDP_INCLUDE_DIR})
ELSE(CSDP_FOUND)
  SET(CSDP_LIBRARIES)
  SET(CSDP_INCLUDE_DIRS)
ENDIF(CSDP_FOUND)

