#ifndef RESPONSE
#define RESPONSE
#include <helper.h>
#include <vector>
#include <string>
#include <map>

using namespace std;

class ttbar;
class TDirectory;

class TTBarResponse
{
	protected:
		string prefix_;
		ttbar* an_;
		TDirectory* dir;
		TH1DCollection plot1d;	
		TH2DCollection plot2d;	
	  map<string, double> defaults_;
  	map<string, pair<double, double> > values_; //gen, reco
  	double weight_=1;

	public:
		TTBarResponse(string prefix, ttbar* an);
		~TTBarResponse();
  	void AddMatrix(string name, double default_val, const vector<double>& Mbins, const vector<double>& Tbins, string label);
	  void FillTruth(string name, double val, double weight);
		void FillReco(string name, double val, double weight);
		void FillTruthReco(string name, double tval, double rval, double weight);
		void Flush();
};

#endif
