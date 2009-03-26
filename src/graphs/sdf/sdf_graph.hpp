#ifndef __SDF_GRPAH_HPP__
#define __SDF_GRPAH_HPP__

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/filesystem.hpp>


namespace bgl = boost;

namespace boost{
	enum edge_bondnum_t {bondnum=112};
	enum graph_classid_t {graph_classid=111};
	BOOST_INSTALL_PROPERTY(graph, classid);
	BOOST_INSTALL_PROPERTY(edge, bondnum);
};

struct SDFGraph{
		typedef bgl::adjacency_list<
			bgl::listS, // store out-edges of vertex in std::list
			bgl::vecS, // store vertex-set in std::vector
			bgl::undirectedS, // not directed
			bgl::property<bgl::vertex_name_t, std::string>,  // vertex-properties
			bgl::property<bgl::edge_weight_t, double,        // edge-properties
				bgl::property<bgl::edge_bondnum_t, int>
			>,       
			bgl::property<bgl::graph_name_t, std::string,    // graph-properties
			  bgl::property<bgl::graph_classid_t, int>
			>
		>    Graph;
		typedef bgl::property_map<Graph, bgl::graph_name_t>::type NameMap;
		typedef bgl::property_map<Graph, bgl::graph_classid_t>::type ClassIDMap;
		typedef bgl::property_map<Graph, bgl::vertex_name_t>::type VertexNameMap;
		typedef bgl::property_map<Graph, bgl::edge_weight_t>::type EdgeWeightMap;
		typedef bgl::graph_traits<Graph>::vertex_descriptor Vertex;
		typedef bgl::graph_traits<Graph>::edge_descriptor Edge;
		typedef bgl::graph_traits<Graph>::vertex_iterator VertexIterator;
		typedef bgl::graph_traits<Graph>::adjacency_iterator AdjIterator;
		typedef bgl::graph_traits<Graph>::edge_iterator EdgeIterator;
		typedef bgl::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;

		int   mClassID;
		Graph mGraph;
		std::map<std::string,std::string> mProps;

		template <class Archive>
			void serialize(Archive& ar, const unsigned int version){
				ar & mGraph;
				ar & mClassID;
				ar & mProps;
			}
};


#endif /* __SDF_GRPAH_HPP__ */

