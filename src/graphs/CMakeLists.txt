SET(SRC graph_cfg.cpp
	sdf/sdf_graph.cpp  sdf/sdf_graph.hpp
	sdf/sdf_reader.cpp sdf/sdf_reader.hpp
	)
include_directories(
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/configuration
)
add_library(graphs SHARED ${SRC})
