// 2009-02-08 Hannes Schulz <mail at hannes-schulz dot de>
#include <fstream>
#include <list>
#include <numeric>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include "count.hpp"
#include <instantiate.hpp>
#include <sdf/sdf_graph.hpp>
#include <sdf/sdf_reader.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/function.hpp>
#include <boost/ref.hpp>
#include <progressbar.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <matlab_io.hpp>

using namespace std;
namespace fs = boost::filesystem;
namespace ublas = boost::numeric::ublas;

namespace bgl = boost;
using boost::format;
using boost::tie;

struct cnt_edges
{
	map<int, int>    mEdges;
	template <class E, class G>
	void operator()(const E& e, const G& g){
		int b  = bgl::get(bgl::edge_bondnum_t(),g,e);
		mEdges[b] ++;
	}
};
struct cnt_atoms
{
	map<string, int> mAtoms;
	template <class V, class G>
	void operator()(const V& v, const G& g){
		string n  = bgl::get(bgl::vertex_name_t(),g,v);
		mAtoms[n] ++;
	}
};

struct cnt_base{
	int mCnt[2];
	ublas::vector<double> mEmbPos;
	int mID;
	static int sMaxID;
	vector<int> mChildIDs;
	cnt_base(){
		memset(mCnt,0, sizeof(mCnt));
		mId = sMaxID ++;
	}
	inline double p_bar(int total_graph_num){ // p(this_feature)
		int c = accumulate(mCnt, mCnt+sizeof(mCnt)/sizeof(int), 0.0);
		return (double)c/total_graph_num;
	}
	inline double p_bar(int klass, int total_graph_num){ // p(this_feature, klass)
		int c = mCnt[klass];
		return (double)c/total_graph_num;
	}
};
int cnt_base::sMaxID=0;


struct cnt_chain : public cnt_base
{
	vector<string>  mA;
	vector<int>     mE;
	vector<SDFGraph::Vertex> mVisitedVertices;
	vector<SDFGraph::Edge>   mVisitedEdges;
	template <class Archive>
	void serialize(Archive& ar, const unsigned int version){
		ar & mA;
		ar & mE;
		ar & mCnt;
	}

	cnt_chain(){ }
	cnt_chain(const vector<string>& a, const vector<int>& e):mA(a), mE(e){ }

	template <class V, class G>
	bool check_v(const V& v, const G& g, size_t i){
		if(i>=mA.size()) return true;
		if(find(mVisitedVertices.begin(), mVisitedVertices.end(), v) != mVisitedVertices.end())
			return false;
		string n  = bgl::get(bgl::vertex_name_t(),g,v);
		if(n!=mA[i] && mA[i]!="*") return false;
		if(mE.size() <= i) return true; // don't check further
		mVisitedVertices.push_back(v);
		typedef typename bgl::graph_traits<G>::out_edge_iterator out_edge_iterator;
		typedef typename bgl::graph_traits<G>::vertex_descriptor vertex_descriptor;
		out_edge_iterator eit_end, eit;
		mVisitedVertices.push_back(v);
		for(tie(eit, eit_end) = bgl::out_edges(v,g); eit!=eit_end; ++eit){
			if(check_e(eit, g, i)) return true;
		}
		mVisitedVertices.pop_back();
		return false;
	}
	template <class E, class G>
	bool check_e(const E& e, const G& g, size_t i){
		if(find(mVisitedEdges.begin(), mVisitedEdges.end(), *e) != mVisitedEdges.end()) return false;
		if(i>=mE.size()) return true;
		int b = bgl::get(bgl::edge_bondnum_t(),g,*e);
		if(mE[i]!=b && mE[i]!=0) return false;
		typedef typename bgl::graph_traits<G>::vertex_descriptor vertex_descriptor;
		vertex_descriptor v = bgl::target(*e, g);
		mVisitedEdges.push_back(*e);
		bool res = check_v(v, g, i+1);
		mVisitedEdges.pop_back();
		return res;
	}
	template <class V, class G>
	inline bool operator()(V v, const G& g, int klass){
		if(check_v(v,g, 0)) {
			if(klass == 0)
				mCnt[0]++;
			else
				mCnt[1]++;
			return true;
		}
		mVisitedEdges.clear();
		mVisitedVertices.clear();
		return false;
	}
};



ostream& operator<<(ostream& o, const cnt_chain& c){
	for(unsigned int i=0;i<c.mA.size();i++){
		o << c.mA[i] << "-";
		if(i<c.mE.size())
			o << c.mE[i] <<"-";
	}
	o << ":  "<<c.mCnt[0] <<" / " << c.mCnt[1];
	return o;
}

vector<cnt_chain> 
getFreqs(const SDFReader& sdf_read, bool verbose, bool force){
	vector<cnt_chain> ccs;
	fs::path freq_fn = gCfg().getString("output-dir");
	freq_fn /= "freqs.ser";
	if(fs::exists(freq_fn) && !force){
		ifstream ifs(freq_fn.string().c_str());
		boost::archive::binary_iarchive ar(ifs);
		ar >> ccs;
		return ccs;
	}
	cnt_atoms ca;
	cnt_edges ce;
	
	// find all atoms and all edge types
	BOOST_FOREACH(const SDFGraph& gobj, sdf_read){
		const SDFGraph::Graph& g = gobj.mGraph;
		SDFGraph::VertexIterator vit, vit_end;
		for(tie(vit, vit_end) = bgl::vertices(g); vit!=vit_end; ++vit){
			ca(*vit, g);
		}

		SDFGraph::EdgeIterator eit, eit_end;
		for(tie(eit, eit_end) = bgl::edges(g); eit!=eit_end; ++eit){
			ce(*eit, g);
		}
	}

	list<cnt_chain> lastround;
	lastround.push_back(cnt_chain(vector<string>(), vector<int>()));

	// create counters for atom-edge pair
	typedef pair<string, int> sipair;
	typedef pair<int, int>    iipair;
	for(int i=0; i<20; i++){
		list<cnt_chain> tmp;
		BOOST_FOREACH(const cnt_chain& c_old, lastround){
			if(c_old.mA.size()>c_old.mE.size()){
				BOOST_FOREACH(const iipair& pe, ce.mEdges){ // insert edge
					cnt_chain c(c_old.mA, c_old.mE);
					c.mE.push_back(pe.first);
					tmp.push_back(c);
					c.mA.push_back("*");                    // and arbitrary atom
					tmp.push_back(c);
				}
			}
			else{
				BOOST_FOREACH(const sipair& pa, ca.mAtoms){
					cnt_chain c(c_old.mA, c_old.mE);        // insert atom
					c.mA.push_back(pa.first);
					tmp.push_back(c);
					c.mE.push_back(0);                      // and arbitrary edge
					tmp.push_back(c);
				}
			}
		}

		// find all occurrences 
		ProgressBar pb(sdf_read.size(), (format("CountLevel %d, t:%d")%i%tmp.size()).str());
		BOOST_FOREACH(const SDFGraph& gobj, sdf_read){
			if(verbose) pb.inc();
			const SDFGraph::Graph& g = gobj.mGraph;
			SDFGraph::VertexIterator vit, vit_end;
			for(tie(vit, vit_end) = bgl::vertices(g); vit!=vit_end; ++vit){
				BOOST_FOREACH(cnt_chain& f, tmp){
					if(f(*vit, g,gobj.mClassID)) break;
				}
			}
		}
		if(verbose) pb.finish();
		lastround.clear();
		// remove non-occurring items from tmp
		BOOST_FOREACH(cnt_chain& c, tmp){
			if(c.mCnt[0]+c.mCnt[1]>=6){
				ccs.push_back(c);
				lastround.push_back(c);
			}
		}
		if(verbose) cout << "Found "<<lastround.size() <<" new entries on level "<<i<<endl;
	}

	// now try to find frequent conjunctions
	if(0){
		list<cnt_chain> tmp;
		ProgressBar pb(ccs.size()*ccs.size(), "Counting Conjunctions");
		BOOST_FOREACH(cnt_chain& c1, ccs){
			BOOST_FOREACH(cnt_chain& c2, ccs){
				pb.inc();
				cnt_chain conj;
				BOOST_FOREACH(const SDFGraph& gobj, sdf_read){
					const SDFGraph::Graph& g = gobj.mGraph;
					SDFGraph::VertexIterator vit, vit_end;
					bool b1=false, b2=false;
					for(tie(vit, vit_end) = bgl::vertices(g); vit!=vit_end; ++vit){
						if(!b1 && c1(*vit, g,gobj.mClassID)) b1=true;
						if(!b2 && c2(*vit, g,gobj.mClassID)) b2=true;
						if(b1&&b2) break;
					}
					if(b1 && b2)  conj.mCnt[gobj.mClassID]++;
				}
				if(conj.mCnt[0]+conj.mCnt[1]>10)
					tmp.push_back(conj);
			}
		}
		pb.finish();
		copy(tmp.begin(),tmp.end(),back_inserter(ccs));
	}
#if 0
		ProgressBar pb(ccs.size()*ccs.size(), "Checking for Instantiations");
		BOOST_FOREACH(cnt_chain& c1, ccs){
			stringstream ss;
			for(int i=0;i<c1.mA.size();i++){
				ss << c1.mA[i];
				if(i<c1.mE.size()) ss << c1.mE[i];
			}
			BOOST_FOREACH(cnt_chain& c2, ccs){
			}
		}
#endif

	ofstream os(freq_fn.string().c_str());
	boost::archive::binary_oarchive ar(os);
	ar << const_cast<const vector<cnt_chain>&>(ccs);

	return ccs;
}


void Count::operator()()
{
	SDFReader sdf_read;
	sdf_read.configure();
	bool force_reread = gCfg().getBool("count.recount");
	vector<cnt_chain> freqs = getFreqs(sdf_read, mVerbose, force_reread);

	fs::path pxy_fn = gCfg().getString("output-dir");
	pxy_fn /= "pxy.m";
	ofstream pxy_os(pxy_fn.string().c_str());
	
	ublas::matrix<double> pxy(2,freqs.size());
	int i=0;
	BOOST_FOREACH(const cnt_chain& c, freqs){
		pxy(0,i) = c.mCnt[0];
		pxy(1,i) = c.mCnt[1];
		i++;
	}
	matlab_matrix_out(pxy_os, "pxy_data", pxy);
}

void Count::configure()
{
	Action::configure();
	mVerbose = gCfg().getBool("verbose");
}

Count::Count()
{
}

Count::~Count()
{
  // cleanup
}

namespace{ registerInFactory<Action, Count> registerBase("count"); }
