cmake_minimum_required( VERSION 2.6 FATAL_ERROR )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

add_library( erlconfig SHARED configuration.cpp configuration.hpp )
target_link_libraries( erlconfig ${Boost_LIBRARIES})
