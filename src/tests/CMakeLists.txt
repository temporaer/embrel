ENABLE_TESTING()

SET(Boost_USE_STATIC_LIBS ON)
find_package( Boost 1.35 COMPONENTS unit_test_framework program_options )

LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS}   )
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/configuration ${Boost_INCLUDE_DIRS} )

ADD_EXECUTABLE(test_rprop test_rprop.cpp )
TARGET_LINK_LIBRARIES ( test_rprop ${Boost_LIBRARIES} erlconfig erltools )
ADD_TEST(test_rprop ./test_rprop )
