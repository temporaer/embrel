// vim:fdm=syntax
#include <numeric>
#include <boost/shared_ptr.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <boost/assign.hpp>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include <csvreader.hpp>
#include <arg_max.hpp>

#include <signal.h>

#include <factory/factory.h>
#include <matlab_io.hpp>
#include <matlab_matio.hpp>
#include <stats.hpp>
#include <normalize.hpp>
#include <progressbar.hpp>
#include <configuration.hpp>

#include <rcode.hpp>

#include <unistd.h>

#include "code_data_gen.hpp"

#define V(X) #X << "="<< (X) << " "

using namespace std;
namespace fs = boost::filesystem;
namespace ublas = boost::numeric::ublas;
namespace ll = boost::lambda;
using namespace boost::assign;
#define foreach BOOST_FOREACH 


CoocReader::~CoocReader(){}
CoocReader::CoocReader(int lines)
	: CSVReader(false)
	, mLines(lines)
	, mProgress(lines,"reading")
{
}
CoocReader::CoocReader()
	:CSVReader(false)
	 ,mLines(0)
{
}
void CoocReader::read_field(const std::string& s, int lidx, int idx, int numfields){
	if(lidx==0 && idx==0){
		cout << "mObsFeatMat size: " << mLines << " x " << (numfields-2) << endl;
		mObsFeatMat.reset(new matrix_itype(mLines,numfields-2));
		*mObsFeatMat = ublas::zero_matrix<int>(mLines,numfields-2);
	}
	if(0);
	else if(idx==0){
		mCurrObsDesc.mDesc = s;
	}
	else if(idx==numfields-1){
		if(s=="act")
			mCurrObsDesc.mKlass = 1;
		else if (s=="inact")
			mCurrObsDesc.mKlass = 0;
		else
			mCurrObsDesc.mKlass = boost::lexical_cast<klass_type>(s);
		//if(mCurrObsDesc.mKlass < 2) mCurrObsDesc.mKlass = 0;
		//else                        mCurrObsDesc.mKlass = 1;
		mKlasses[mCurrObsDesc.mKlass] ++;
		mObsDesc += mCurrObsDesc;

		mProgress.inc();
		if(lidx+1==mLines)
			mProgress.finish(false);
	}
	else{
		int i= (int)(0.5 + boost::lexical_cast<double>(s));
		if(i>0){
			(*mObsFeatMat)(lidx,idx-1) = i;
		//	mObsFeatMat->push_back(lidx,idx-1,i);
		}
	}
}

void CoocReader::init_features(int aleph_queries, vector<string>& aleph_names)
{
	map<observation::klass_type, int> klass_cnt;
	map<observation::klass_type, double> klass_prior;
	if(mObsDesc.size() != mObsFeatMat->size1())
		throw runtime_error(string("number of observation inconsistency"));
	foreach(const observation& o, mObsDesc){
		klass_cnt[o.mKlass] ++;
	}
	cout << "Number of Classes: " << klass_cnt.size()<<endl;
	typedef pair<observation::klass_type, double> kd_p;
	typedef pair<observation::klass_type, int>    ki_p;
	foreach(const kd_p& p, klass_cnt){
		klass_prior[p.first] = (double)klass_cnt[p.first] / mObsFeatMat->size2();
	}
	fs::path featnamfile(gCfg().getString("code.input_file") + "-names.csv");
	fs::ifstream featnamstream;
	featnamstream.open(featnamfile);
	int running_feat_num = 0;

	ProgressBar pb1(mObsFeatMat->size2(),"feat props");
	ExactDescriptiveStatistics estat("Entropy");
	boost::regex key_re1("\\bkey\\(A\\),\\s*");
	boost::regex key_re2("\\bA,\\s*");
	//boost::regex muta_soln(
		//"(?:attyp\(., ., 195\))"
		//,boost::regex::mod_x);
	bool regression = klass_cnt.size() >= 5;
	reverse(aleph_names.begin(), aleph_names.end());
	for(unsigned int i=0;i<mObsFeatMat->size2();i++){
		pb1.inc();
		feature f;
		// id
		string str;
		f.mAlephFeat = ( i >= mObsFeatMat->size2() - aleph_queries);
		if(f.mAlephFeat){
			f.mId = aleph_names.back();
			aleph_names.pop_back();
		}else{
			assert(!featnamstream.eof());
			getline(featnamstream, str);
			boost::trim(str);
			vector<string> tmp;
			boost::split( tmp, str, boost::is_any_of(":"));
			f.mId       = tmp[0];
			f.mFromId   = tmp.size()>1 ? tmp[1] : "_";
			f.mTargetId = tmp.size()>2 ? tmp[2] : "_";
			f.mId = boost::regex_replace(f.mId, key_re1, "");
			f.mId = boost::regex_replace(f.mId, key_re2, "");
		}

		f.mComplexity = f.mId.length();

		// class count
		ublas::matrix_column<matrix_itype> mc(*mObsFeatMat,i);
		for(unsigned int j=0;j<mc.size();j++){
			f.mKlassCount[mObsDesc[j].mKlass] += mc[j];
		}

		// frequency
		f.mFreq = 0;
		foreach(const ki_p& p,f.mKlassCount){
			f.mFreq += p.second;
		}

		// entropy
		f.mEntropy = 0;
		foreach(const ki_p& p,f.mKlassCount){
			double x = ((double)p.second / f.mFreq) * klass_prior[p.first];
			if(x<=0) continue;
			f.mEntropy += x * log(x);
		}
		estat += f.mEntropy;

		// bestklass
		double maxval;
		f.mBestKlass = util::arg_max(
				f.mKlassCount.begin(), 
				f.mKlassCount.end(),maxval,
				ll::bind<int>(&kd_p::second,ll::_1))->first;

		// precision & recall
		float precision,recall;
		ExactDescriptiveStatistics regress_stats;
		if(regression){
			for(map<observation::klass_type,int>::iterator p=f.mKlassCount.begin();p!=f.mKlassCount.end();p++){
				for(int i=0;i<p->second;i++)
					regress_stats.notify(p->first);
			}
			precision = regress_stats.getMean();
			recall    = regress_stats.getN();
		}else{
			precision = (float) f.mKlassCount[f.mBestKlass] / f.mFreq;
			recall    = (float) f.mKlassCount[f.mBestKlass] / klass_cnt[f.mBestKlass];
		}
		const float beta= 0.2f;
		f.mFMeasure     = (1+beta*beta) * (precision * recall) / (beta*beta*precision + recall);

		//f.mId           = (boost::format("%s p:%1.2f r:%1.2f") % f.mId % precision % recall).str();
		f.mId           = (boost::format("%s") % f.mId ).str();

		// color factor
		f.mColorFact = precision;

		// size 
		f.mSize        = recall;
		//f.mSize      = (double) log((double)f.mFreq);
		//f.mSize      = (double)f.mEntropy * log((double)f.mFreq);
		//f.mSize      = f.mFMeasure;

		// selection criterion
		//f.mSelectCrit = f.mFreq / log(f.mComplexity);
		f.mSelectCrit = f.mFMeasure / f.mComplexity;

		f.mIgnore = false;
		f.mRunningNumber = running_feat_num++;

		mFeaDesc += f;
	}
	pb1.finish();

if(gCfg().getBool("code.remove_sim")){
	cout << "remove very similar features..."<<endl;
	ProgressBar pbpf(mObsFeatMat->size2() * mObsFeatMat->size2() / 2 - mObsFeatMat->size2()/2,"Prefilt");
	vector<int> colsToRetain;
	for(unsigned long int i=0;i<mObsFeatMat->size2();i++){
		bool retain=true;
		for(unsigned long int j=i+1;j<mObsFeatMat->size2();j++)
		{
			const ublas::matrix_column<matrix_itype> a(*mObsFeatMat,i);
			const ublas::matrix_column<matrix_itype> b(*mObsFeatMat,j);
			if(equal(a.begin(),a.end(),b.begin())){
				assert(i<mFeaDesc.size());
				assert(j<mFeaDesc.size());
				if(mFeaDesc[i].mComplexity <= mFeaDesc[j].mComplexity){
					if(mFeaDesc[i].mAlephFeat)
						mFeaDesc[j].mAlephFeat=true; // conserve alephfeat flag
					swap(mFeaDesc[i],mFeaDesc[j]);
				}
				mFeaDesc[i].mRunningNumber=-1; // mark for deletion
				retain=false;
				break;
			}
		}
		if(retain) colsToRetain += i;
		pbpf.inc(mObsFeatMat->size2()-i);
	}
	pbpf.finish();
	cout << "Need to delete "<< (mObsFeatMat->size2()-colsToRetain.size()) << " features."<<endl;
	cout << "New feature number: "<< colsToRetain.size() << " features."<<endl;
	cout << "deleting..."<<flush;
	mFeaDesc.erase(
		std::remove_if(mFeaDesc.begin(),mFeaDesc.end(),ll::bind(&feature::mRunningNumber,ll::_1)<0),mFeaDesc.end());
	cout << "..."<<flush;
	matrix_pitype tmp(new matrix_itype(mObsFeatMat->size1(),mFeaDesc.size()));
	int newidx=0;
	foreach(int i, colsToRetain){
		const ublas::matrix_column<matrix_itype> source(*mObsFeatMat,i);
		      ublas::matrix_column<matrix_itype> target(*tmp,newidx);
		ublas::noalias(target) = source;
		newidx++;
	}
	mObsFeatMat = tmp;
	// renumber everything
	newidx=0;
	foreach(feature& f, mFeaDesc){ f.mRunningNumber=newidx++; }
	cout << "done."<<endl;
	
	//exit(0);
}


	unsigned int fea_num = mFeaDesc.size();
	cout << "creating f x f matrix..."<<flush;
	mFeatFeatMat.reset(new matrix_dtype(fea_num,fea_num)); 
	cout << "zeroing..."<<flush;
	*mFeatFeatMat = ublas::zero_matrix<double>(fea_num,fea_num);
	cout << "done."<<endl;
	
	if(!gCfg().getBool("code.dont_run_code")){
#if 1
		cout << "f x f... multiplication" <<flush;
		ublas::noalias(*mFeatFeatMat) = ublas::prod(ublas::trans(*mObsFeatMat), *mObsFeatMat);
		//ublas::axpy_prod(ublas::trans(*mObsFeatMat), *mObsFeatMat, *mFeatFeatMat, true);
		cout << "done."<<endl;
#else
		ProgressBar pb(fea_num*fea_num/2+fea_num/2, "f x f");
		matrix_dtype& featfeat = *mFeatFeatMat;
		matrix_itype& obsfeat = *mObsFeatMat;
		for(unsigned int i=0; i<featfeat.size1();i++){
			for(unsigned int j=i+1; j<featfeat.size2();j++){
				int f=0;
				for(unsigned int k=0;k<obsfeat.size1();k++)
					if( obsfeat(k,i)>0 && obsfeat(k,j)>0 )
						f++;
				if(f==0) continue;
				featfeat.insert_element(i,j,f);
				featfeat.insert_element(j,i,f);
			}
			featfeat.insert_element(i,i,mFeaDesc[i].mFreq); // diagonal
			pb.inc(fea_num-i);
		}
		pb.finish();
#endif
	}
}

template<class T>
struct PosReader : CSVReader{
	T& mVec;
	PosReader(T& v)
		:CSVReader(false," ")
		,mVec(v)
	{
	}
	void read_field(const std::string& s, int lidx, int idx, int numfields){
		if(idx==0)
			mVec[lidx].mPos = ublas::vector<double>(numfields);
		mVec[lidx].mPos[idx] = boost::lexical_cast<double>(s);
	}
};

void
CoocReader::load_obs_pos(const string& fn){
	PosReader<vector<observation> > pr(mObsDesc);
	ifstream obs_stream(fn.c_str());
	pr.read(obs_stream);
	obs_stream.close();
}
void
CoocReader::load_feat_pos(const string& fn){
	ifstream fea_stream(fn.c_str());
	PosReader<vector<feature> >     pr(mFeaDesc);
	pr.read(fea_stream);
	fea_stream.close();
}


CODE_data_gen::~CODE_data_gen(){
}

void CODE_data_gen::configure(){
}


void CODE_data_gen::run(){
	// read input file
	boost::filesystem::path infile(gCfg().getString("code.input_file") + ".csv");
	boost::filesystem::path serfile(gCfg().getString("code.input_file") + ".ser");
	cout << "Reading CSV: " << infile<<endl;
	fs::ifstream ifs;

	ifs.open(infile);

	int count=0;
	string str;
	for(;getline(ifs,str);count++);
	ifs.close();
	CoocReader cr(count);

	bool read_cr = false;
	if(fs::exists(serfile)){
		ifstream ifs(serfile.string().c_str());
		boost::archive::binary_iarchive ar(ifs);
		int saved_count=0;
		ar >> saved_count;
		if(count==saved_count){
			ar >> cr;
			read_cr = true;
		}
		cout <<endl<< "Read data from cache; ObsFeatMat:" <<endl
		     <<cr.getObsFeatMat()->size1() << "x"<<cr.getObsFeatMat()->size2()<<endl
				 <<"FeatureNum: " << cr.mFeaDesc.size()<<endl
				 <<"ObsNum:     " << cr.mObsDesc.size()<<endl;
	}
	if(!read_cr){
		ifs.open(infile);
		cr.read(ifs);
		ifs.close();
		int running_obs_num=0;
		foreach(observation& o, cr.mObsDesc){ o.mRunningNumber=running_obs_num++;}

		int aleph_feat = 0;
		vector<string> aleph_names;
		string alephfn = gCfg().getString("code.aleph_queries");
		if(alephfn.length()>0){
			cout << "Reading ALEPH queries from "<<alephfn<< "..."<<flush;
			ifstream ais(alephfn.c_str());
			string descr,line1, line2;
			int count=0;
			for(;getline(ais,str);count++);
			ais.close();
			count /=3;

			// adjust size of  matrix
			unsigned int s1 = cr.getObsFeatMat()->size1();
			unsigned int s2 = cr.getObsFeatMat()->size2();
			CoocReader::matrix_itype m(s1, s2+count);
			m.clear();
			ublas::matrix_range<CoocReader::matrix_itype> mr (m, ublas::range (0, s1), ublas::range (0, s2));
			mr = *cr.getObsFeatMat();

			ais.open(alephfn.c_str());
			while(!ais.eof()){
				getline(ais, descr);
				boost::trim(descr);
				if(!descr.length())
					break;
				getline(ais, line1); // 1st line: plain query
				getline(ais, line2); // 2nd line: w/o previous queries
				//cout << V(descr) << V(line1) <<endl;
				boost::trim(line1); boost::trim(line2);
				vector<string> tmp;
				boost::split( tmp, line1, boost::is_any_of(","));
				foreach(const string& s, tmp){
					unsigned int obs = boost::lexical_cast<int>(s)-1;
					m(obs,s2+aleph_feat) = 1;
				}
				aleph_names.push_back(descr);
				aleph_feat ++;
			}
			* cr.getObsFeatMat()  = m;
			cout << "done, "<< aleph_feat << " read."<<endl;
		}

		cr.init_features(aleph_feat, aleph_names);

		ofstream ofs(serfile.string().c_str());
		boost::archive::binary_oarchive ar(ofs);
		ar << count;
		ar << cr;
	}



	// this block only for NIPS
#if 0
	ifstream fnis("/home/schulzha/data/nips_counts-names.csv");
	foreach(feature& f, cr.mFeaDesc){
		string str;
		getline(fnis, str);
		boost::trim(str);
		f.mId   = str;
		f.mSize = f.mFreq;
		f.mEntropy = 1.0;
		f.mColorFact = f.mFreq;
		f.mSelectCrit = f.mFreq;
	}
#endif
	

	// call matlab.
	RCode rc(2);
	rc.configure();
	rc.setPxy(ublas::trans(*cr.getObsFeatMat()));
	rc.setPxx(*cr.getFeatFeatMat());
	rc.init_positions(); // randomly
	if(!gCfg().getBool("code.dont_run_code")){
		if(!gCfg().getBool("code.entropy_emb")){
			// straight forward embedding
			rc.mFixXpos = false;
			rc.mFixYpos = false;
			run_code(rc,cr,false);
		}else{
			// entropybased 2nd embedding after 1st
			rc.mFixXpos = false;
			rc.mFixYpos = false;
			run_code(rc,cr,false); // use only feat pos

			cout << "weighing cooccurrence by squared normalized entropy..."<<flush;
			ExactDescriptiveStatistics estat;
			estat.notify(cr.mFeaDesc.begin(), cr.mFeaDesc.end(), ll::bind(&feature::mEntropy,ll::_1));
			for(unsigned int i=0;i<cr.getObsFeatMat()->size2();i++){
				ublas::matrix_column<CoocReader::matrix_itype> col(*cr.getObsFeatMat(),i);
				col *= pow(normalize_minmax(cr.mFeaDesc[i].mEntropy,0.0,1.0,estat),2);
			}
			for(unsigned int i=0;i<cr.getObsFeatMat()->size1();i++){
				ublas::matrix_row<CoocReader::matrix_itype> row(*cr.getObsFeatMat(),i);
				row /= accumulate(row.begin(),row.end(),0.0);
			}
			rc.setPxy(ublas::trans(*cr.getObsFeatMat()));
			rc.setPxx(*cr.getFeatFeatMat());
			cout <<"done."<<endl;

			rc.mUse_Pxx = false;
			rc.mUse_Pxy = true;
			rc.mFixXpos = true;
			rc.mFixYpos = false;
			run_code(rc,cr,true); // load previously determined feat pos
		}
	}

	if(gCfg().getBool("code.linlog_write")){
		cout << "Writing Linlog Graph to file..."<<flush;
		ofstream linlog(gCfg().getOutputFile("code.linlog_graph_out").c_str());
		CoocReader::matrix_itype& ofm = *cr.getObsFeatMat();
		CoocReader::matrix_dtype& ffm = *cr.getFeatFeatMat();
		for(unsigned int i=0;i<ofm.size1(); i++){
			for(unsigned int j=0;j<ofm.size2(); j++){
			  if(ofm(i,j)> .001)
					linlog << j << " " << 100000+i <<" "<< 1.0 <<endl;
			}
		}
		for(unsigned int i=0;i<ffm.size1(); i++){
			for(unsigned int j=0;j<ffm.size2(); j++){
			  if(ffm(i,j)> .001)
					linlog << j << " " << i <<" "<< ffm(i,j) <<endl;
			}
		}
		linlog.close();
		cout << "done. Exiting. "<<endl<<endl
		     << "You should now run LinLogLayout from http://code.google.com/p/linloglayout/"<<endl
				 << "  $  java -Xmx2000m -Xms1000m -cp ../../LinLogLayout/bin LinLogLayout 2 " 
				      << gCfg().getOutputFile("code.linlog_graph_out") << " "
							<< gCfg().getOutputFile("code.linlog_pos_in") << endl
				 << "Afterwards, run erl with --code.linlog_read --code.linlog_pos_in="<<gCfg().getOutputFile("code.linlog_pos_in")<<endl;
		exit(0);
	}

	if(gCfg().getBool("code.linlog_read")){
		cout << "Reading linlog positions..." <<flush;
	  ifstream linlog(gCfg().getOutputFile("code.linlog_pos_in").c_str());
		int id;
		float x,y,z, clus;
		while(!linlog.eof()){
			linlog >> id >> x >> y >> z >> clus;
			if(id >= 100000){
			  id -= 100000;
				cout << id << " " << x<< " " << y <<endl;
				cr.mObsDesc[id].mPos[0]=x;
				cr.mObsDesc[id].mPos[1]=y;
			}else{
				cout << id << " " << x<< " " << y <<endl;
				cr.mFeaDesc[id].mPos[0]=x;
				cr.mFeaDesc[id].mPos[1]=y;
			}
		}
		cout << "done."<<endl;
	} else{
		cout << "Reading CODE positions..." << flush;
		cr.load_feat_pos(gCfg().getOutputFile("fea.txt").c_str());
		cr.load_obs_pos(gCfg().getOutputFile("cla.txt").c_str());
		cout << "done."<<endl;
	}

	CoocReader::featpos fp;
	foreach(feature& f, cr.mFeaDesc){
		fp.mFeat = &f;
		cr.mFeatTree.insert(fp);
	}
	cr.mFeatTree.optimize();

	ProgressBar pb2(cr.mFeaDesc.size(), "registering");
	unsigned int num_neigh = gCfg().getInt("code.contrast_eq_neigh");
	foreach(feature& f, cr.mFeaDesc){
		pb2.inc();
		fp.mFeat = &f;
		vector<CoocReader::featpos> res;
		float range = 0.1;
		do{
			res.clear();
			cr.mFeatTree.find_within_range(fp,range,back_inserter(res));
			range *= 1.1;
		}while(res.size()<num_neigh);

		float sum=0;
		foreach(const CoocReader::featpos& p, res){
			sum += p.mFeat->mSelectCrit;
		}
		f.mSelectCrit2 = f.mSelectCrit - sum/res.size();
	}
	pb2.finish();

	int maxiter = gCfg().getInt("code.hebb_iter");

	ProgressBar pb3(maxiter,"hebb");
	ExactDescriptiveStatistics sc("selectcrit");
	foreach(feature& f, cr.mFeaDesc){ sc += f.mSelectCrit2; }
	foreach(feature& f, cr.mFeaDesc){ f.mSelectCrit2 = normalize_minmax(f.mSelectCrit2, 0.0,1.0, sc); }
	for(int iter=0;iter<maxiter;iter++){
		pb3.inc();
		foreach(feature& query, cr.mFeaDesc){
			CoocReader::featpos fp;
			fp.mFeat = &query;
			vector<CoocReader::featpos> res;
			float range = 0.01;
			float tmp = -1E6;
			do{
				res.clear();
				cr.mFeatTree.find_within_range(fp,range,back_inserter(res));
				range *= 1.1;
			}while(res.size()<num_neigh);
			fp = *util::arg_max(res.begin(), res.end(), tmp,
					ll::bind<float&>(&feature::mSelectCrit2,*ll::bind<feature*>(&CoocReader::featpos::mFeat,ll::_1)));
			foreach(CoocReader::featpos& p, res){
				if(fp.mFeat == p.mFeat){
					p.mFeat->mSelectCrit2 *= 1.01;
					p.mFeat->mSelectCrit2 = min(p.mFeat->mSelectCrit2, 1000.f);
				}
				else{
					p.mFeat->mSelectCrit2 *= 0.998;
				}
			}
		}
	}
	pb3.finish();

	vector<feature*> fvec(cr.mFeaDesc.size());
	transform(cr.mFeaDesc.begin(),cr.mFeaDesc.end(), fvec.begin(), &ll::_1);
	sort(fvec.begin(),fvec.end(), 
			ll::bind(&feature::mSelectCrit2,ll::_1) <
			ll::bind(&feature::mSelectCrit2,ll::_2));
	int threshidx = (int)(gCfg().getFloat("code.delete_perc") * cr.mFeaDesc.size());
	float thresh = fvec[threshidx]->mSelectCrit2;
	int ignoredcnt=0;
	foreach(feature& f, cr.mFeaDesc){ 
		if(f.mSelectCrit2 >= thresh)
			continue;
		f.mIgnore = true;
		ignoredcnt++;
	}
	
	ExactDescriptiveStatistics xstat("x"), ystat("y"), cstat("color"), sstat("size");
	foreach(feature& f, cr.mFeaDesc){
		//if(f.mIgnore) continue;
		xstat += f.mPos[0];
		ystat += f.mPos[1];
		cstat += f.mColorFact;
		sstat += f.mSize;
	}
	ExactDescriptiveStatistics klass_stats;
	foreach(observation& o, cr.mObsDesc){
		xstat += o.mPos[0];
		ystat += o.mPos[1];
		klass_stats += o.mKlass;
	}
	float img_width  = gCfg().getFloat("code.img_width");
	float img_height = gCfg().getFloat("code.img_height");
	float size_fact  = gCfg().getFloat("code.size_fact");
	float pos_rand   = gCfg().getFloat("code.pos_rand");
	ExactDescriptiveStatistics posnormstat = xstat.getRange() > ystat.getRange() ? xstat : ystat; 
	posnormstat += xstat.getMin();
	posnormstat += xstat.getMax();
	posnormstat += ystat.getMin();
	posnormstat += ystat.getMax();
	foreach(feature& f, cr.mFeaDesc){
		//if(f.mIgnore) continue;
		f.mPos[0] = normalize_minmax(f.mPos[0], 0.0f, img_width,  xstat) + (2*(drand48()-0.5))*pos_rand;
		f.mPos[1] = normalize_minmax(f.mPos[1], 0.0f, img_height, ystat) + (2*(drand48()-0.5))*pos_rand;
		//f.mPos[0] = normalize_sd(f.mPos[0], 0.0f, img_width,0.3f,  xstat) + (2*(drand48()-0.5))*pos_rand;
		//f.mPos[1] = normalize_sd(f.mPos[1], 0.0f, img_height,0.3f, ystat) + (2*(drand48()-0.5))*pos_rand;
		if(cr.mKlasses.size()<=4)
			f.mColorFact = normalize_minmax(f.mColorFact,0.0f,255.0f,cstat);
		else // regression
		 ;
		f.mSize = normalize_minmax(f.mSize,size_fact*4.0f,size_fact*10.0f,sstat);
	}

	fs::path pointsxml(gCfg().getOutputFile("points.xml").c_str());
	fs::ofstream xml(pointsxml);
	xml << "<?xml version='1.0' encoding='utf-8' ?>"<<endl
		<< "<data>"<<endl;
	string color;
	cout << "Writing Features: "<<cr.mFeaDesc.size() <<endl;
	foreach(feature& f, cr.mFeaDesc){
		//if(f.mIgnore) continue;
		string color;
		if(cr.mKlasses.size()<4){
			switch(f.mBestKlass){
				case 0: color="0x%02X0101"; break;
				case 1: color="0x01%02X01"; break;
				case 2: color="0x0101%02X"; break;
				default: throw runtime_error(string("too few colors available!"));
			}
		}else{ 
			//regression
			const float splitpoint = 5.f;
			if(f.mColorFact <= splitpoint){
				f.mColorFact = 255-normalize_minmax(f.mColorFact,150.0f,255.0f,cstat);
				color="0x%02X0101";
			}else{
				f.mColorFact = normalize_minmax(f.mColorFact,0.0f,255.0f,cstat);
				color="0x01%02X01";
			}
		}
		color = (boost::format(color)%(int)f.mColorFact).str();
		xml << boost::format(
				"<node id='%s' ignore='%d' color='%s' size='%1.0f' x='%03.0f' y='%03.0f' objtype='%d' alephfeat='%d' >")
			%  f.mId  % f.mIgnore % color   % f.mSize   % f.mPos[0] % f.mPos[1] % 1 % f.mAlephFeat
			<< endl;
		writeFeaChildren(xml,cr,f.mRunningNumber);
		xml << "</node>"<<endl;
	}
	vector<string> obscolor;
	bool regression=false;
	if(cr.mKlasses.size() <= 4)
		obscolor += "0xFFAAAA","0x0000FF", "0xAA0808", "0x0000AA";
	else if(cr.mKlasses.size() <= 5)
		obscolor += "0x000000", "0x00003E", "0x00007C", "0x0000B2", "0x0000FF";
	else{
		regression=true;
	}
	foreach(observation& o, cr.mObsDesc){
		o.mPos[0] = normalize_minmax(o.mPos[0], 0.0f, img_width,  xstat);
		o.mPos[1] = normalize_minmax(o.mPos[1], 0.0f, img_height, ystat);
		//o.mPos[0] = normalize_sd(o.mPos[0], 0.0f, img_width, 0.3f,  xstat);
		//o.mPos[1] = normalize_sd(o.mPos[1], 0.0f, img_height, 0.3f, ystat);
		string color;
		if(regression){
			int c;
			if(o.mKlass<5)
				c = normalize_minmax(o.mKlass,0.0,255.0, klass_stats);
			else
				c = normalize_minmax(o.mKlass,150.0,255.0, klass_stats);
			color = (boost::format("0x%02X0101") %c).str();
		}else{
			color = obscolor[o.mKlass];
		}
		xml << boost::format(
				"<node id='obs-%d' color='%s' size='%1.0f' x='%03.0f' y='%03.0f' objtype='%d'>")
			%      o.mRunningNumber     % color       % 4    % o.mPos[0]  % o.mPos[1] % 0 
			<< endl;
		writeObsChildren(xml,cr,o.mRunningNumber);
		xml << "</node>"<<endl;
	}

	xml << "</data>"<<endl;
	xml.close();

	if(gCfg().getBool("code.dvc")){
		ofstream dist_vs_norm(gCfg().getString("code.dvc_file").c_str());
		ProgressBar pb2a(cr.mFeaDesc.size(),"dist_vs_cooc");
		foreach(feature& f1, cr.mFeaDesc){
			foreach(feature& f2, cr.mFeaDesc){
				dist_vs_norm << pow(ublas::norm_2(f1.mPos-f2.mPos),2.0);
				dist_vs_norm << " ";
				ublas::matrix_column<CoocReader::matrix_itype> c1(*cr.getObsFeatMat(),f1.mRunningNumber);
				ublas::matrix_column<CoocReader::matrix_itype> c2(*cr.getObsFeatMat(),f2.mRunningNumber);
				double c=0;
				for(unsigned int i=0;i<c1.size();i++) {
					if(c1[i]>0 && c2[i]>0)
						c++;
				}
				dist_vs_norm << c;
				dist_vs_norm << endl;
			}
			pb2a.inc();
		}
		pb2a.finish();
		dist_vs_norm.close();
	}


	int res=0;
	res |= system((boost::format("cd ../.. ; perl pointsxml2svg.pl %s")%gCfg().getOutputFile("points.xml")).str().c_str());
	res |= system((boost::format("convert %s %s")%gCfg().getOutputFile("points-static.svg")%gCfg().getOutputFile("points.png")).str().c_str());
	//res |= system("/usr/bin/play -q /usr/lib/xcdroast/sound/test.wav 2>&1>/dev/null");
	if(res!=0){
		cout<<"error during system(...)"<<endl;
	}
}

void CODE_data_gen::run_code(RCode& rc, CoocReader& cr, bool load_fea_pos){
#define USE_RCODE 1
#if USE_RCODE
		double lastLogLik=-1E9;
		int n_restarts = gCfg().getInt("code.nrestarts");
		for(int i=0;i<n_restarts;i++){
			rc.init_positions(); // randomly
			if(load_fea_pos){
				cr.load_feat_pos(gCfg().getOutputFile("fea.txt").c_str());
				for(unsigned int i=0;i<cr.mFeaDesc.size();i++){
					ublas::row(rc.mXpos,i) = cr.mFeaDesc[i].mPos;
				}
			}
			double loglik = rc.run(i);
			if(loglik<lastLogLik) continue;
			lastLogLik = loglik;
			ofstream os1(gCfg().getOutputFile("fea.txt").c_str());
			for(unsigned int i=0;i<rc.mXpos.size1();i++) {
				ublas::matrix_row<RCode::mat_t> r(rc.mXpos,i);
				copy(r.begin(),r.end(),ostream_iterator<double>(os1," "));
				os1<<endl;
			}
			ofstream os2(gCfg().getOutputFile("cla.txt").c_str());
			for(unsigned int i=0;i<rc.mYpos.size1();i++) {
				ublas::matrix_row<RCode::mat_t> r(rc.mYpos,i);
				copy(r.begin(),r.end(),ostream_iterator<double>(os2," "));
				os2<<endl;
			}
			os1.close(), os2.close();
		}
			
#else
		CoocReader::matrix_itype& obsfea = *cr.getObsFeatMat();
		CoocReader::matrix_dtype& feafea = *cr.getFeatFeatMat();
		cout << "writing matrices to matlab file..."<<flush;
		if(matlab_matrix_out("/tmp/code_data.mat","feat_feat",feafea))
			throw runtime_error(string("could not write feat_feat"));
		CoocReader::matrix_itype tmpmat(ublas::trans(obsfea));
		if(matlab_matrix_out("/tmp/code_data.mat","feat_klass",tmpmat))
			throw runtime_error(string("could not write feat_klass"));
		cout <<"done."<<endl;

		chdir("../../src/matlab");
		const char* matlab_out = "/tmp/matlab.out";
		int res = system(
				//(boost::format("MALLOC_CHECK_=1 matlab -glnxa64 -nosplash -nodisplay -nojvm -r eval_codtest -logfile %s") % matlab_out).str().c_str()
				(boost::format("matlab -glnx86 -nosplash -nodisplay -nojvm -r eval_codtest -logfile %s") % matlab_out).str().c_str()
				);
		if(res == -1)
			throw runtime_error(std::string("Matlab execution failed!"));
		if(WIFSIGNALED(res) && (WTERMSIG(res) == SIGINT || WTERMSIG(res) == SIGQUIT))
			throw runtime_error(std::string("Got interrupted."));
#endif
}

double CODE_data_gen::getLogLik(RCode& rc, CoocReader& cr){
	rc.init_positions(); // make sure we set mPosInitialized
	for(unsigned int i=0;i<cr.mFeaDesc.size();i++){
		ublas::row(rc.mXpos,i) = cr.mFeaDesc[i].mPos;
	}
	for(unsigned int i=0;i<cr.mObsDesc.size();i++){
		ublas::row(rc.mYpos,i) = cr.mObsDesc[i].mPos;
	}
	double loglik = 0;
	int  maxRPropIters = 1;
	swap(rc.mRPropMaxIter, maxRPropIters);
	loglik = rc.run(0);
	swap(rc.mRPropMaxIter, maxRPropIters);
	return loglik;
}

void CODE_data_gen::writeObsChildren(ostream& os, CoocReader&cr, int obsnr)
{
	CoocReader::matrix_itype& ofmat = *cr.getObsFeatMat();
	for(unsigned int f=0;f<ofmat.size2();f++){
		if(ofmat(obsnr,f) == 0)
			continue;
		//if(cr.mFeaDesc[f].mIgnore)
			//continue;
		os << boost::format("<child id='%s' color='0xFF0000'/>")%cr.mFeaDesc[f].mId;
	}
}
void CODE_data_gen::writeFeaChildren(std::ostream& os, CoocReader&cr, int feanr)
{
	CoocReader::matrix_itype& ofmat = *cr.getObsFeatMat();
	for(unsigned int o=0;o<ofmat.size1();o++){
		if(ofmat(o,feanr)==0)
			continue;
		os << boost::format("<child id='obs-%d' color='0xFF0000'/>")%cr.mObsDesc[o].mRunningNumber;
	}
}


namespace{ registerInFactory<Action, CODE_data_gen> registerBase("code"); }
