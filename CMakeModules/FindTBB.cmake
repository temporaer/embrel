#
# Locate tbb include paths and libraries
# tbb can be found at http://osstbb.intel.com/index.php
# Written by Manfred Ulz, manfred.ulz_at_tugraz.at

# This module defines
# TBB_INCLUDE_DIR, where to find parallel_for.h, etc.
# TBB_LIBRARIES, the libraries to link against to use libtbb.
# TBB_FOUND, If false, don't try to use libtbb.

FIND_PATH(TBB_INCLUDE_DIR parallel_for.h
  PATH_SUFFIXES
    tbb
  PATHS
    /usr/local/include
    /usr/include
    ${SOOFEA_SOURCE_DIR}/include
)

FIND_LIBRARY(TBB_LIBRARIES
  NAMES 
    tbb
  PATHS
    /usr/local/lib
    /usr/lib
    ${SOOFEA_SOURCE_DIR}/lib
)

SET(TBB_FOUND 0)
IF(TBB_INCLUDE_DIR)
  IF(TBB_LIBRARIES)
    SET(TBB_FOUND 1)
    MESSAGE(STATUS "Found TBB")
  ENDIF(TBB_LIBRARIES)
ENDIF(TBB_INCLUDE_DIR)

MARK_AS_ADVANCED(
  TBB_INCLUDE_DIR
  TBB_LIBRARIES
) 