#ifndef __SDF_READER_HPP__
#define __SDF_READER_HPP__
#include "sdf_graph.hpp"

class SDFReader{
	private:
		int          mFixedSize;
	public:

		SDFReader();
		virtual void configure();
		virtual ~SDFReader();
		
		typedef      std::vector<SDFGraph>::iterator iterator;
		typedef      std::vector<SDFGraph>::const_iterator const_iterator;
		inline iterator       begin()     {return mGraphs.begin();}
		inline const_iterator begin()const{return mGraphs.begin();}
		inline iterator       end  ()     {return mGraphs.end();}
		inline const_iterator end  ()const{return mGraphs.end();}
		inline size_t         size ()const{return mGraphs.size();}
	private:
		int mOutputCounter;

		struct FileDescriptor{
			std::string name;
			int         classid;
		};
		std::vector<FileDescriptor> mInputFiles; 

		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive& ar, const unsigned int version){
			ar & mGraphs;
		}

		std::vector<SDFGraph>  mGraphs;

		void readInputFiles();
		bool readMolekule(std::istream& is,int klass, int cnt);
		void setInputFiles(const std::string& s);
};
#endif /* __SDF_READER_HPP__ */

