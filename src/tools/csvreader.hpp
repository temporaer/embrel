#ifndef __CSVREADER_HPP__
#define __CSVREADER_HPP__

#include <string>
#include <iostream>

class CSVReader {
	private:
		bool mReadHeader;
		const char* mDivider;
	public:
		CSVReader(bool readHeader=false, const char* div = ",");
		virtual void read(std::istream& is);
		virtual void read_line(const std::string& s, int l_idx);
		virtual void read_header(const std::string& s, int lidx, int idx, int numfields);
		virtual void read_field (const std::string& s, int lidx, int idx, int numfields);
};
#endif /* __CSVREADER_HPP__ */
