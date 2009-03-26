// vim:fdm=syntax
#include <numeric>
#include <boost/shared_ptr.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/assign.hpp>

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
#include <stats.hpp>
#include <normalize.hpp>
#include <progressbar.hpp>
#include <configuration.hpp>

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
{
}
CoocReader::CoocReader()
	:CSVReader(false)
	 ,mLines(0)
{
}
void CoocReader::read_field(const std::string& s, int lidx, int idx, int numfields){
	if(lidx==0 && idx==0){
		mObsFeatMat.reset(new matrix_itype(mLines,numfields-2));
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
	}
	else{
		(*mObsFeatMat)(lidx,idx-1) = boost::lexical_cast<int>(s);
	}
}

void CoocReader::init_features()
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
	for(unsigned int i=0;i<mObsFeatMat->size2();i++){
		pb1.inc();
		feature f;
		// id
		string str;
		assert(!featnamstream.eof());
		getline(featnamstream, str);
		boost::trim(str);
		vector<string> tmp;
		boost::split( tmp, str, boost::is_any_of(":"));
		f.mId       = tmp[0];
		f.mFromId   = tmp.size()>1 ? tmp[1] : "_";
		f.mTargetId = tmp.size()>2 ? tmp[2] : "_";


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
		float precision = (float) f.mKlassCount[f.mBestKlass] / f.mFreq;
		float recall    = (float) f.mKlassCount[f.mBestKlass] / klass_cnt[f.mBestKlass];
		float beta      = 0.2f;
		f.mFMeasure     = (1+beta*beta) * (precision * recall) / (beta*beta*precision + recall);

		// color factor
		f.mColorFact = (float)f.mKlassCount[f.mBestKlass] / (float)f.mFreq;

		// size 
		//f.mSize      = (double)f.mEntropy * log((double)f.mFreq);
		f.mSize      = f.mFMeasure;

		// selection criterion
		//f.mSelectCrit = f.mFreq / log(f.mId.length());
		f.mSelectCrit = f.mFMeasure / f.mId.length();

		f.mIgnore = false;
		f.mRunningNumber = running_feat_num++;

		mFeaDesc += f;
	}
	pb1.finish();

#if 0
	// weight cooccurrence by entropy
	for(unsigned int i=0;i<mObsFeatMat->size2();i++){
		ublas::matrix_column<matrix_itype> col(*mObsFeatMat,i);
		col *= normalize_minmax(mFeaDesc[i].mEntropy,0.0,1.0,estat);
		col /= accumulate(col.begin(),col.end(),0.0);
	}
#endif

	unsigned int fea_num = mFeaDesc.size();
	mFeatFeatMat.reset(new matrix_dtype(fea_num,fea_num,2*fea_num/3)); // last param relevant only for sparse matrices, otherwise overwritten by zero_matrix!
	*mFeatFeatMat = ublas::zero_matrix<double>(fea_num,fea_num);
	
	if(!gCfg().getBool("code.dont_run_code")){
#if 1
		cout << "f x f..." ;
		ublas::noalias(*mFeatFeatMat) = ublas::prod(ublas::trans(*mObsFeatMat), *mObsFeatMat);
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


CODE_data_gen::~CODE_data_gen(){
}

void CODE_data_gen::configure(){
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

void CODE_data_gen::run(){
	// read input file
	boost::filesystem::path infile(gCfg().getString("code.input_file") + ".csv");
	cout << "Reading CSV: " << infile<<endl;
	fs::ifstream ifs;

	ifs.open(infile);

	int count=0;
	string str;
	for(;getline(ifs,str);count++);
	ifs.close();

	ifs.open(infile);
	CoocReader cr(count);
	cr.read(ifs);
	ifs.close();
	int running_obs_num=0;
	foreach(observation& o, cr.mObsDesc){ o.mRunningNumber=running_obs_num++;}

#if 0
	cout << "remove very similar features..."<<endl;
	ProgressBar pbpf(cr.getObsFeatMat()->size2() * cr.getObsFeatMat()->size2() / 2 - cr.getObsFeatMat()->size2()/2,"Prefilt");
	vector<int> colsToDelete;
	for(int i=0;i<cr.getObsFeatMat()->size2();i++){
		for(int j=i+1;j<cr.getObsFeatMat()->size2();j++)
		{
			int f = ublas::norm_1(
					ublas::column(*cr.getObsFeatMat(),i)-
					ublas::column(*cr.getObsFeatMat(),j));
			if(f>5)
				colsToDelete += j;
		}
		pbpf.inc(cr.getObsFeatMat()->size2()-i);
	}
	pbpf.finish();

	cout << "Need to delete "<< colsToDelete.size() << " features."<<endl;
	//exit(0);
#endif

	cr.init_features();
	CoocReader::matrix_itype& obsfea = *cr.getObsFeatMat();
	CoocReader::matrix_dtype& feafea = *cr.getFeatFeatMat();

	fs::path code_data("/tmp/code_data.m");
	fs::ofstream code_data_stream(code_data);
	matlab_matrix_out(code_data_stream, "feat_feat", feafea);
	matlab_matrix_out(code_data_stream, "feat_klass", ublas::trans(obsfea));

	// call matlab.
	if(!gCfg().getBool("code.dont_run_code")){
		chdir("../../src/matlab");
		const char* matlab_out = "/tmp/matlab.out";
		int res = system(
				(boost::format("matlab -nodisplay -nojvm -r eval_codtest -logfile %s") % matlab_out).str().c_str());
		if(res == -1)
			throw runtime_error(std::string("Matlab execution failed!"));
		if(WIFSIGNALED(res) && (WTERMSIG(res) == SIGINT || WTERMSIG(res) == SIGQUIT))
			throw runtime_error(std::string("Got interrupted."));
	}

	PosReader<vector<feature> >     pr_fea(cr.mFeaDesc);
	PosReader<vector<observation> > pr_obs(cr.mObsDesc);
	ifstream fea_stream("/tmp/erl/fea.txt");
	ifstream obs_stream("/tmp/erl/cla.txt");
	pr_fea.read(fea_stream);
	pr_obs.read(obs_stream);

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

	int maxiter = 100;
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
			float range = 0.1;
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
	int threshidx = (int)(gCfg().getFloat("code.view_sd_fact") * cr.mFeaDesc.size());
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
	foreach(observation& o, cr.mObsDesc){
		xstat += o.mPos[0];
		ystat += o.mPos[1];
	}
	float img_width  = gCfg().getFloat("code.img_width");
	float img_height = gCfg().getFloat("code.img_height");
	foreach(feature& f, cr.mFeaDesc){
		if(f.mIgnore) continue;
		f.mPos[0] = normalize_minmax(f.mPos[0], 0.0f, img_width,  xstat);
		f.mPos[1] = normalize_minmax(f.mPos[1], 0.0f, img_height, ystat);
		f.mColorFact = normalize_minmax(f.mColorFact,0.0f,255.0f,cstat);
		f.mSize = normalize_minmax(f.mSize,8.0f,20.0f,sstat);
	}

	fs::path pointsxml("/tmp/erl/points.xml");
	fs::ofstream xml(pointsxml);
	xml << "<?xml version='1.0' encoding='utf-8' ?>"<<endl
		<< "<data>"<<endl;
	string color;
	foreach(feature& f, cr.mFeaDesc){
		if(f.mIgnore) continue;
		string color;
		if(cr.mKlasses.size()<4){
			switch(f.mBestKlass){
				case 0: color="0x%02X0101"; break;
				case 1: color="0x01%02X01"; break;
				case 2: color="0x0101%02X"; break;
				default: throw runtime_error(string("too few colors available!"));
			}
		}else{
			color="0x%02X0101";
		}
		color = (boost::format(color)%(int)f.mColorFact).str();
		xml << boost::format(
				"<node id='%s' color='%s' size='%1.0f' x='%03.0f' y='%03.0f' objtype='%d'>")
			%          f.mId   % color   % f.mSize   % f.mPos[0] % f.mPos[1] % 1 
			<< endl;
		writeFeaChildren(xml,cr,f.mRunningNumber);
		xml << "</node>"<<endl;
	}
	vector<string> obscolor;
	if(cr.mKlasses.size() <= 4)
		obscolor += "0x0000FF", "0xFFAAAA", "0x0000AA", "0xAA0808";
	else if(cr.mKlasses.size() <= 5)
		obscolor += "0x000000", "0x00003E", "0x00007C", "0x0000B2", "0x0000FF";
	else
		throw runtime_error(string("not enough class colors available"));
	foreach(observation& o, cr.mObsDesc){
		o.mPos[0] = normalize_minmax(o.mPos[0], 0.0f, img_width,  xstat);
		o.mPos[1] = normalize_minmax(o.mPos[1], 0.0f, img_height, ystat);
		string color = obscolor[o.mKlass];
		xml << boost::format(
				"<node id='obs-%d' color='%s' size='%1.0f' x='%03.0f' y='%03.0f' objtype='%d'>")
			%      o.mRunningNumber     % color       % 4    % o.mPos[0]  % o.mPos[1] % 0 
			<< endl;
		writeObsChildren(xml,cr,o.mRunningNumber);
		xml << "</node>"<<endl;
	}

	xml << "</data>"<<endl;
	xml.close();

	system("cd /home/sarx/prog/uni/ma/molemb1 ; perl pointsxml2svg.pl /tmp/erl/points.xml");
	system("convert /tmp/erl/points-static.svg /tmp/erl/points.png");
	system("/usr/bin/play -q /usr/lib/xcdroast/sound/test.wav 2>&1>/dev/null");
}

void CODE_data_gen::writeObsChildren(ostream& os, CoocReader&cr, int obsnr)
{
	CoocReader::matrix_itype& ofmat = *cr.getObsFeatMat();
	for(unsigned int f=0;f<ofmat.size2();f++){
		if(ofmat(obsnr,f) == 0)
			continue;
		if(cr.mFeaDesc[f].mIgnore)
			continue;
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
