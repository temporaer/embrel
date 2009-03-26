

#include <sstream>
#include <fstream>
#include <numeric>
#include <configuration.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/graphviz.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

#include <factory/factory.h>
#include <matlab_io.hpp>
#include "sdf_reader.hpp"
#include <nana.h>

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;
namespace fs = boost::filesystem;
namespace ll = boost::lambda;

bool nextLineMatches(istream& is, string re){
	string s;
	getline(is, s);
	regex e(re);
	bool b = regex_match(s, e);
	if(!b){
		cerr << "Warning: Line `"<<s<<"' does not match expected pattern `"<<re<<"'"<<endl;
		return false;
	}
	return true;
}

void SDFReader::configure()
{
	string files = gCfg().getString("SDFReader.files");
	mFixedSize   = gCfg().getInt("SDFReader.fixed_size");
	setInputFiles(files);
	if(mGraphs.size()==0){
		if(fs::exists("/tmp/erl/graphs.ser") && gCfg().getBool("SDFReader.use_cache")){
			ifstream ifs("/tmp/erl/graphs.ser");
			archive::binary_iarchive ar(ifs);
			ar >> (*this);
		}else
			readInputFiles();
		mOutputCounter = 0;
	}
}

SDFReader::SDFReader()
	:mFixedSize(0)
{
}
SDFReader::~SDFReader()
{
}

void SDFReader::setInputFiles(const std::string& inputfileline)
{
	string line = boost::trim_copy(inputfileline);

	vector<string> file_class_pairs;
	boost::split( file_class_pairs, line, boost::is_any_of(",") );

	BOOST_FOREACH(const string& pair, file_class_pairs){
		vector<string> tmp;
		boost::split( tmp, pair, boost::is_any_of(":") );
		if(!tmp.size()==2)
		{
			cerr << "Warning: SDFReader: InputFormatError: "<<pair<< " in " << line<< " should have format filename:classid"<<endl;
			continue;
		}
		FileDescriptor fd;
		fd.name = tmp[0];
		fd.classid = boost::lexical_cast<int>(tmp[1]);
		mInputFiles.push_back(fd);
	}
}

void SDFReader::readInputFiles()
{
	int graphCount = 0;
	BOOST_FOREACH(const FileDescriptor& fd, mInputFiles){
		L("Reading SDF input file `%s'...", fd.name.c_str());
		ifstream is(fd.name.c_str());
		int fcnt = 0;
		while(readMolekule(is, fd.classid, graphCount++) && fcnt++<400);
		L(". %d molecules read.\n", graphCount);
	}
	ofstream os("/tmp/erl/graphs.ser");
	archive::binary_oarchive ar(os);
	random_shuffle(mGraphs.begin(), mGraphs.end());
	ar << const_cast<const SDFReader&>(*this);
}

bool SDFReader::readMolekule(std::istream& is, int klass, int graphCount)
{
	boost::regex whitespace("\\s+");
	SDFGraph desc;
	desc.mClassID = klass;

	string line, id;
	if(!nextLineMatches(is, "\\s*"))            return false;
	if(!nextLineMatches(is, ".*OpenBabel.*"))   return false;
	if(!nextLineMatches(is, "\\s*"))            return false;
	if(!nextLineMatches(is, ".*V3000"))         return false;
	if(!nextLineMatches(is, ".*BEGIN\\s+CTAB")) return false;
	getline(is, line); // M V30 COUNTS numberAtoms numberBonds 0 0 0
	boost::trim(line);
	vector<string> strvec;
	boost::sregex_token_iterator tokit(line.begin(), line.end(), whitespace, -1), tokit_end;
	copy(tokit, tokit_end, back_inserter(strvec));
	if(!strvec.size()==8){
		cerr << "Warning: SDFReader: InputFileFormat Error: input description line `"<<line<<"' has unexpected format."<<endl;
		return false;
	}
	if(!nextLineMatches(is, ".*BEGIN\\s+ATOM")) return false;

	// create a matrix with numberAtoms size
	int numberAtoms = boost::lexical_cast<int>(strvec[3]);
	int numberBonds = boost::lexical_cast<int>(strvec[4]);
	int n = (mFixedSize==0) ? numberAtoms : mFixedSize;

	desc.mGraph            = SDFGraph::Graph(n);
	SDFGraph::Graph& graph = desc.mGraph;

	// read all atoms
	SDFGraph::VertexIterator vi, vi_end;
	int atomCnt=0;
	for(tie(vi,vi_end) = bgl::vertices(graph); vi!=vi_end; ++vi, ++atomCnt){
		if(atomCnt == numberAtoms) break;
		getline(is, line);
		boost::trim(line);
		strvec.clear();
		boost::sregex_token_iterator tokit(line.begin(), line.end(), whitespace, -1), tokit_end;
		copy(tokit, tokit_end, back_inserter(strvec));

		if(strvec.size() < 8){
			cerr << "Warning: SDFReader: InputFileFormat Error: atom "<<graphCount<<" description line `"<<line<<"' has unexpected format."<<endl;
			return false;
		}
		string atomType = strvec[3]; // C or N or...
		bgl::get(bgl::vertex_name_t(),graph)[*vi] = atomType;
	}
	if(!nextLineMatches(is,".*END\\s+ATOM"))return false;
	if(!nextLineMatches(is,".*BEGIN\\s+BOND"))return false;
	
	// read all bonds
	for(int i=0;i<numberBonds; i++){
		getline(is, line);
		boost::trim(line);
		strvec.clear();
		boost::sregex_token_iterator tokit(line.begin(), line.end(), whitespace, -1), tokit_end;
		copy(tokit, tokit_end, back_inserter(strvec));

		if(strvec.size() < 6){
			cerr << "Warning: SDFReader: InputFileFormat Error: bond description line `"<<line<<"' has unexpected format."<<endl;
			return false;
		}
		int atomFrom = boost::lexical_cast<int>(strvec[4])-1;
		int atomTo   = boost::lexical_cast<int>(strvec[5])-1;
		int bondType = boost::lexical_cast<int>(strvec[3]);
		if(atomFrom >= numberAtoms || atomTo >= numberAtoms){
			cerr << "Warning: SDFReader: InputFileFormat Error: bond description line `"<<line<<"' has unexpected format."<<endl;
			return false;
		}
		bgl::add_edge(atomFrom, atomTo, graph);
		bgl::get(bgl::get(bgl::bondnum, graph),bgl::edge(atomFrom, atomTo, graph).first) = bondType;
	}
	if(!nextLineMatches(is,".*END\\s+BOND"))return false;
	if(!nextLineMatches(is,".*END\\s+CTAB"))return false;
	if(!nextLineMatches(is,".*END"))return false;
	
	getline(is,line);
	regex prop_match(">\\s+<(\\w+)>");
	boost::cmatch matches;
	while(line[0] == '>'){ // additional properties
		boost::trim(line);
		if(line.size() == 0) break;
		if(regex_match(line.c_str(),matches,prop_match)){
			string       propname = matches[1];
			stringstream propval;
			while(line.size()>0){
				getline(is,line);
				boost::trim(line);
				propval << line;
			}
			desc.mProps[propname] = propval.str();
		}else{
			cerr << "Warning: Line `"<<line<<"' does not match expected pattern `> <feature>'"<<endl;
			return false;
		}
		getline(is, line);
	}
	if(line != "$$$$"){
		cerr << "Warning: Line `"<<line<<"' does not match expected pattern `$$$$'"<<endl;
		return false;
	}

	stringstream sname;
	sname << "SDF_class"<<klass<<"_"<<graphCount;
	bgl::get_property(graph, bgl::graph_classid) = klass;

	mGraphs.push_back(desc);

	return true;
}

