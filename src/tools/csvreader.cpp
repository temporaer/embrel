#include <vector>

#include <boost/algorithm/string.hpp>

#include <boost/foreach.hpp>

#include "csvreader.hpp"
#define foreach BOOST_FOREACH 

CSVReader::CSVReader(bool readHeader, const char* div)
	:mReadHeader(readHeader)
	,mDivider(div)
{
}

void CSVReader::read(std::istream& is){
	std::string line;
	int idx=0;
	while(is){
		getline(is,line);
		boost::trim(line);
		read_line(line,idx++);
	}
}
void CSVReader::read_line(const std::string& s, int l_idx){
	std::vector<std::string> fields;
	boost::split( fields, s, boost::is_any_of(mDivider));
	int f_idx=0;
	fields.erase(remove(fields.begin(), fields.end(), ""),fields.end());
	int numfields=fields.size();

	foreach(std::string& s, fields){
		boost::trim(s);
		if(l_idx==0 && mReadHeader){
			read_header(s, l_idx, f_idx, numfields);
		}else{
			read_field(s, l_idx, f_idx, numfields);
		}
		f_idx++;
	}

}

void CSVReader::read_field(const std::string& s, int lidx, int idx, int numfields){
}
void CSVReader::read_header(const std::string& s, int lidx, int idx, int numfields){
}

