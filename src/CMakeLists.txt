cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )
SET(Boost_ADDITIONAL_VERSIONS "1.34" "1.34.1" "1.35" "1.35.1" "1.36" "1.36.1" )

# General stuff
PROJECT(EmbRel CXX Fortran C)
SET(CMAKE_PROJECT_NAME "EmbRel" )
SET(SERALIGN_SCRIPT_DIR ${CMAKE_BINARY_DIR}/scripts)

# find errors/eclipse auto-incl-dir discovery
#SET(CMAKE_VERBOSE_MAKEFILE ON)


# Packaging
SET(CPACK_PACKAGE_NAME                EmbRel )
SET(CPACK_PACKAGE_CONTACT             "Hannes Schulz <mail@hannes-schulz.de>" )
SET(CPACK_PACKAGE_VERSION             0.1 )
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Embedded Relational Learning" )
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386 )
SET(CPACK_DEBIAN_PACKAGE_PRIORITY     optional )
INCLUDE(CPack)

# No linebreaks in errors (gcc) (easier parsing)
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0 ")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0 ")
	MESSAGE(STATUS ${CMAKE_CXX_FLAGS})
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
	IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
		# Debug mode: Use Safe Containers
		#ADD_DEFINITIONS(-D_GLIBCXX_DEBUG)
	ELSE()
		# Other modes: Use Parallel Algorithms
		# turn off at single uses like so:
		# std::sort(v.begin(), v.end(), __gnu_parallel::sequential_tag());
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -march=native" )
		#ADD_DEFINITIONS(-D_GLIBCXX_PARALLEL) # does not work
	ENDIF()
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


# Precompiled Headers
find_package(PCHSupport QUIET)

# Boost
find_package( Boost 1.35 COMPONENTS program_options filesystem system regex serialization  REQUIRED )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} configuration )

# uBlas Backend
#FIND_PACKAGE(BLAS REQUIRED)
#FIND_PACKAGE(LAPACK REQUIRED)
include_directories ( ${LAPACK_INCLUDE_DIRS} )
link_directories    ( ${LAPACK_LIBRARY_DIRS} )
FIND_PATH(BOOST_BINDINGS_PATH boost/numeric/bindings/lapack/lapack.hpp PATHS /usr/local/include/boost-numeric-bindings )
include_directories( ${BOOST_BINDINGS_PATH} )

# Nana
find_package( Nana )
include_directories ( ${NANA_INCLUDE_DIRS} )

# general third party
include_directories ( ${CMAKE_SOURCE_DIR}/third_party )

# subdirs: sub-libs
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( actions )
ADD_SUBDIRECTORY( configuration )
ADD_SUBDIRECTORY( graphs )
ADD_SUBDIRECTORY( third_party/gnuplot-cpp )
ADD_SUBDIRECTORY( flash/abstree )
include_directories( 
  ${CMAKE_SOURCE_DIR}/tools
  ${CMAKE_SOURCE_DIR}/actions
  ${CMAKE_SOURCE_DIR}/graphs
  ${CMAKE_SOURCE_DIR}/configuration 
  )

# subdirs: optional extensions
#ADD_SUBDIRECTORY( python_bindings EXCLUDE_FROM_ALL ) 
ADD_SUBDIRECTORY( tests EXCLUDE_FROM_ALL ) 

# docs
find_package(Doxygen QUIET)
if(DOXYGEN_EXECUTABLE)
	CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/../doc/Doxyfile.in
			${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
			IMMEDIATE)
	ADD_CUSTOM_TARGET(api-docs ${DOXYGEN_EXECUTABLE} doc/Doxyfile)
endif(DOXYGEN_EXECUTABLE)

MESSAGE(STATUS "Boost Libs: ${Boost_LIBRARIES} ${Boost_LIBRARY_DIRS}" )
add_executable( erl main.cpp )
FIND_PROGRAM( CTAGS_BINARY ctags )
IF( CTAGS_BINARY )
	add_custom_command (TARGET erl
		POST_BUILD COMMAND ctags -R --c++-kinds=+p --fields=+iaS --extra=+q --exclude=third_party
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMENT "Creating tags" VERBATIM)
ENDIF( CTAGS_BINARY )
FIND_PROGRAM( MKID_BINARY mkid )
#IF( MKID_BINARY )
#    add_custom_command (TARGET erl 
#        POST_BUILD COMMAND mkid --include "C C++"
#        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} COMMENT "Creating global-tags" VERBATIM)
#ENDIF( MKID_BINARY )
target_link_libraries(erl
	${Boost_LIBRARIES}  
	${NANA_LIBRARIES} 
	${LAPACK_LIBRARIES} 
	dl
	erltools
	graphs
	erlconfig 
	)
