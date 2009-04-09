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
#include <boost/array.hpp>
#include <progressbar.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <matlab_io.hpp>
#include <stats.hpp>
#include <normalize.hpp>

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
	vector<int> mCnt;
	ublas::vector<double> mEmbPos;
	int mID;
	float mCoor[2], mColor;
	float mSize;
	static int sMaxID;
	vector<int> mChildIDs;
	int         mParent;
	cnt_base()
		:mCnt(100,0)
	{
		mID = sMaxID ++;
	}
	template <class Archive>
	void serialize(Archive& ar, const unsigned int version){
		ar & mID;
		ar & mParent;
		ar & mChildIDs;
	}
	inline double p_bar(int total_graph_num){ // p(this_feature)
		int c = accumulate(mCnt.begin(), mCnt.end(), 0.0);
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
	mutable vector<SDFGraph::Vertex> mVisitedVertices;
	mutable vector<SDFGraph::Edge>   mVisitedEdges;
	template <class Archive>
	void serialize(Archive& ar, const unsigned int version){
		ar & boost::serialization::base_object<cnt_base>(*this);
		ar & mA;
		ar & mE;
		ar & mCnt;
	}

	cnt_chain(){ }
	cnt_chain(const vector<string>& a, const vector<int>& e):mA(a), mE(e){ }

	string str()const{
		int i = 0;
		stringstream ss;
		while(i<mA.size()){
			ss << mA[i];
			if(i<mE.size()){
				switch(mE[i]){
					case 0: ss << "."; break;
					case 1: ss << "-"; break;
					case 2: ss << "="; break;
					default:ss << "?"; break;
				}
			}
			i++;
		}
		int cnt = accumulate(mCnt.begin(), mCnt.end(),0.0);
		ss << "\t" << cnt;
		return ss.str();
	}

	template <class V, class G>
	bool check_v(const V& v, const G& g, size_t i)const{
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
	bool check_e(const E& e, const G& g, size_t i)const{
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
	inline bool operator()(V v, const G& g, int klass)const{
		if(check_v(v,g, 0)) {
			return true;
		}
		mVisitedEdges.clear();
		mVisitedVertices.clear();
		return false;
	}
	template <class V, class G>
	inline bool operator()(V v, const G& g, int klass){
		if(check_v(v,g, 0)) {
			mCnt[klass] ++;
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

void printCountsToCSV(const char* fn,const SDFReader& sdf_read, const vector<cnt_chain>& css,bool verbose=true){
	ofstream os(fn);

	ProgressBar pb(sdf_read.size(), "Writing CSV...");
	int graphidx=0;
	BOOST_FOREACH(const SDFGraph& gobj, sdf_read){
		if(verbose) pb.inc();
		ublas::vector<int> feat_matches(css.size(),0);
		const SDFGraph::Graph& g = gobj.mGraph;
		SDFGraph::VertexIterator vit, vit_end;
		int feat = 0;
		BOOST_FOREACH(const cnt_chain& f, css){
			for(tie(vit, vit_end) = bgl::vertices(g); vit!=vit_end; ++vit){
				if(f(*vit, g,gobj.mClassID)){
					feat_matches[feat] = 1;
					break;
				}
			}
			feat++;
		}
		os << "graph("<<graphidx++ <<"),";
		copy(feat_matches.begin(), feat_matches.end(), ostream_iterator<int>(os,","));
		os << gobj.mClassID << endl;
	}
	ofstream featn("/tmp/erl/chains-names.csv");
	BOOST_FOREACH(const cnt_chain& f, css){
		featn << f.str() << ":" <<f.mID << ":"<<f.mParent <<endl;
	}
	pb.finish();
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
#if 0
	typedef pair<string, int> si_pair;
	BOOST_FOREACH(si_pair p, ca.mAtoms){
		cout << p.first << " " << p.second<<endl;
	}
	typedef pair<int, int> ii_pair;
	BOOST_FOREACH(ii_pair p, ce.mEdges){
		cout << p.first << " " << p.second<<endl;
	}
#endif

	vector<cnt_chain*> lastround;
	cnt_chain dummy_chain = cnt_chain(vector<string>(),vector<int>());
	lastround.push_back(&dummy_chain);
	int minFreq = gCfg().getInt("count.min_freq");
	int maxLevel = gCfg().getInt("count.max_level");

	// create counters for atom-edge pair
	typedef pair<string, int> sipair;
	typedef pair<int, int>    iipair;
	for(int i=0; i<maxLevel; i++){
		list<cnt_chain> tmp;
		BOOST_FOREACH(cnt_chain* c_old, lastround){
			if(c_old->mA.size() -1 == c_old->mE.size()){
				BOOST_FOREACH(const iipair& pe, ce.mEdges){ // insert edge
					cnt_chain c(c_old->mA, c_old->mE);
					c.mE.push_back(pe.first);
					c.mParent = c_old->mID;
					c_old->mChildIDs.push_back(c.mID);

					tmp.push_back(c);

					//c = cnt_chain(c_old->mA, c_old->mE);
					//c.mE.push_back(pe.first);
					//c.mA.push_back("*");                    // and arbitrary atom
					//c_old->mChildIDs.push_back(c.mID);
					//tmp.push_back(c);
				}
			}
			else if(c_old->mA.size() == c_old->mE.size()){
				BOOST_FOREACH(const sipair& pa, ca.mAtoms){
					cnt_chain c(c_old->mA, c_old->mE);        // insert atom
					c.mA.push_back(pa.first);
					c.mParent = c_old->mID;
					c_old->mChildIDs.push_back(c.mID);

					// check whether the chain exists (mirrored) already
					bool found_mirr = false;
					BOOST_FOREACH (cnt_chain mirr, tmp){
						bool found_diff = false;
						int as = mirr.mA.size();
						int es = mirr.mE.size();
						for(int i=0;!found_diff && i<as;++i)
							if(mirr.mA[as-i-1] != c.mA[i])
								found_diff=true;
						for(int i=0;!found_diff && i<es;++i){
							if(mirr.mE[es-i-1] != c.mE[i])
								found_diff=true;
						}
						if(!found_diff){
							found_mirr = true;
							//cout << "will not add "<<c.str()<<" bc it mirrors "<<mirr.str()<<endl;
							break;
						}
					}
					if(found_mirr)
						continue;
					tmp.push_back(c);

					//c = cnt_chain(c_old->mA, c_old->mE);
					//c.mA.push_back(pa.first);
					//c.mE.push_back(0);                      // and arbitrary edge
					//c_old->mChildIDs.push_back(c.mID);
					//tmp.push_back(c);
				}
			}
		}
		if(0){
		BOOST_FOREACH(cnt_chain& f, tmp){
			BOOST_FOREACH(string& n, f.mA){
				cout <<"a:"<< n << " ";
			}
			cout <<endl;
			BOOST_FOREACH(int& n, f.mE){
				cout <<"e:"<< n << " ";
			}
			cout <<endl;
			cout <<endl;
		}}


		// find all occurrences 
		ProgressBar pb(sdf_read.size(), (format("CountLevel %d, t:%d")%i%tmp.size()).str());
		BOOST_FOREACH(const SDFGraph& gobj, sdf_read){
			if(verbose) pb.inc();
			const SDFGraph::Graph& g = gobj.mGraph;
			SDFGraph::VertexIterator vit, vit_end;
			BOOST_FOREACH(cnt_chain& f, tmp){
				for(tie(vit, vit_end) = bgl::vertices(g); vit!=vit_end; ++vit){
					if(f(*vit, g,gobj.mClassID)){
						break;
					}
				}
			}
		}
		if(verbose) pb.finish();
		lastround.clear();
		ccs.reserve(ccs.size() + tmp.size());
		// remove non-occurring items from tmp
		BOOST_FOREACH(cnt_chain& c, tmp){
			int cnt = accumulate(c.mCnt.begin(), c.mCnt.end(),0.0);
			if(cnt>=minFreq){
				ccs.push_back(c);
				lastround.push_back(&ccs.back());
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
				if(accumulate(conj.mCnt.begin(),conj.mCnt.end(),0.0)>10)
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


template<class T>
void
determineChildren(T& cont, cnt_base& b){
	bool addedSth = true;
	while(addedSth){
		addedSth = false;
		BOOST_FOREACH(int id, b.mChildIDs){
			BOOST_FOREACH(const cnt_base& c, cont){
				if(c.mID != id) continue;
				vector<int> add;
				BOOST_FOREACH(int nid, c.mChildIDs){
					if(find(b.mChildIDs.begin(), b.mChildIDs.end(), nid) != b.mChildIDs.end())
						continue;
					if(find(add.begin(), add.end(), nid) != add.end())
						continue;
					add.push_back(nid);
					addedSth = true;
				}
				copy(add.begin(),add.end(),back_inserter(b.mChildIDs));
			}
		}
	}
	vector<int> children;
	BOOST_FOREACH(int id, b.mChildIDs){
		BOOST_FOREACH(const cnt_base& c, cont){
			if(c.mID != id) continue;
			children.push_back(id);
			break;
		}
	}
	b.mChildIDs = children;
}

template<class T>
void
makeXML(T& cont){
	ifstream is("/tmp/erl/phix.txt");
	ofstream os("/tmp/erl/points.xml");
	ExactDescriptiveStatistics xstats("xstats"), ystats("ystats"), cstats("cstats"), sstats("sstats");
	double x, y, c;
	BOOST_FOREACH(cnt_base& b, cont){
		is >>x>>y>>c;
		b.mCoor[0]= x;  xstats += x;
		b.mCoor[1]= y;  ystats += y;
		b.mColor  = c;  cstats += c;
		b.mSize   = log(log(accumulate(b.mCnt.begin(), b.mCnt.end(), 0.0)));
		sstats += b.mSize;
	}
	cout << xstats<<ystats<<cstats<<endl;
	
	// normalize
	BOOST_FOREACH(cnt_base& b, cont){
		b.mCoor[0] =       normalize_minmax(b.mCoor[0], 0.0, 800.0, xstats);
		b.mCoor[1] = 600 - normalize_minmax(b.mCoor[1], 0.0, 600.0, ystats);
		b.mColor   =       normalize_minmax(b.mColor,   0.0, 255.0, cstats);
		b.mSize    =       normalize_minmax(b.mSize,    1.0, 7.0,   sstats);
	}

	// print
	os << "<?xml version='1.0' encoding='latin1' ?>"<<endl;
	os << "<data>"<<endl;
	BOOST_FOREACH(cnt_base& b, cont){
		os  << format("<node id='%d' color='0x10%02X10' size='%.0d' x='%.0f' y='%.0f'>")
			                     %b.mID     %((int)b.mColor)   %max(1.0,round(b.mSize))   %b.mCoor[0] %b.mCoor[1]
			<<endl;
		BOOST_FOREACH(int cid, b.mChildIDs){
			os << format("<child id='%d' />") % cid<<endl;
		}
		int maxc = distance(b.mCnt.begin(), max_element(b.mCnt.begin(), b.mCnt.end()));
		os << format("<child id='class%d' />") % maxc<<endl;
		os << "</node>"<<endl;
	}

	ifstream classis("/tmp/erl/psiy.txt");
	int i=0;
	while(!classis.eof()){
		classis >> x >> y;
		if(classis.eof()) break;
		x =       normalize_minmax(x, 0.0, 800.0, xstats);
		y = 600 - normalize_minmax(y, 0.0, 600.0, ystats);
		os  << format("<node id='class%1$d' color='0x%2$02X%2$02X%2$02X' size='%3$d' x='%4$.0f' y='%5$.0f'>")
			                    %i            %(i*40)                        %15      %x         %y
			<<endl;
		os << "</node>"<<endl;
		i++;
	}
	
	os << "</data>"<<endl;
}

void Count::operator()()
{
	SDFReader sdf_read;
	sdf_read.configure();
if(1){
	// assign classes based on properties
	BOOST_FOREACH(SDFGraph& gobj, sdf_read){
		if(gobj.mProps["ActivityCategory_ER_RBA"] == "active strong")
			gobj.mClassID = 4;
		if(gobj.mProps["ActivityCategory_ER_RBA"] == "active medium")
			gobj.mClassID = 3;
		if(gobj.mProps["ActivityCategory_ER_RBA"] == "active weak")
			gobj.mClassID = 2;
		if(gobj.mProps["ActivityCategory_ER_RBA"] == "slight binder")
			gobj.mClassID = 1;
		if(gobj.mProps["ActivityCategory_ER_RBA"] == "inactive")
			gobj.mClassID = 0;
		if(gobj.mProps.find("ActivityScore_NCTRER") != gobj.mProps.end()){
			gobj.mClassID = boost::lexical_cast<int>(gobj.mProps["ActivityScore_NCTRER"]);
		}
	}
}
	bool force_reread = gCfg().getBool("count.recount");
	fs::path pxy_fn   = gCfg().getString("output-dir");
	vector<cnt_chain> freqs = getFreqs(sdf_read, mVerbose, force_reread);
	printCountsToCSV((pxy_fn/"chains.csv").string().c_str(),sdf_read,freqs);

	pxy_fn /= "pxy.m";
	ofstream pxy_os(pxy_fn.string().c_str());
	
	ublas::matrix<double> pxy(freqs.front().mCnt.size(),freqs.size());
	int i=0;
	BOOST_FOREACH(const cnt_chain& c, freqs){
		for(unsigned int j=0;j<c.mCnt.size();j++)
			pxy(j,i) = c.mCnt[j];
		i++;
	}
	matlab_matrix_out(pxy_os, "pxy_data", pxy);
	if(!force_reread){
		BOOST_FOREACH(cnt_chain& c, freqs){
			determineChildren(freqs, c);
		}
		makeXML(freqs);
	}
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
