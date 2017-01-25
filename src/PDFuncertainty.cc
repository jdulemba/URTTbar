#include "Analyses/URTTbar/interface/PDFuncertainty.h"
#include "TDirectory.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"

using namespace std;
using namespace LHAPDF;        
PDFuncertainty::PDFuncertainty(const string setorigname, int memorig, const vector<string>& setnames) : 
  oldx1_(0.),
  use_evt_weights_(false),
  nweights_(-1) {
	setorig_ = new PDFSet(setorigname);
	pdforig_ = setorig_->mkPDF(memorig);
	for(const string& setname : setnames)
	{
		sets_.push_back(new PDFSet(setname));
		pdfs_.push_back(sets_.back()->mkPDFs());
		weights_.push_back(vector<double>(pdfs_.back().size()));
	}
}

PDFuncertainty::PDFuncertainty(int nweights):
  oldx1_(0.),
  use_evt_weights_(true),
  nweights_(nweights),
  setorig_(),
  pdforig_(),
  sets_(),
  pdfs_(),
  weights_() {}

PDFuncertainty::~PDFuncertainty() 
{

}

void PDFuncertainty::SetupWeights(URStreamer& streamer) {
	const Geninfo& info = streamer.genInfo();
	double x1 = info.x1();	
	if(oldx1_ == x1) return;
	oldx1_ = x1;
	double x2 = info.x2();
	double Q = info.renScale();
	int id1 = info.pdfid1();	
	int id2 = info.pdfid2();	
	for(size_t s = 0 ; s < sets_.size() ; ++s) {
		for(size_t p = 0 ; p < pdfs_[s].size() ; ++p) {
			weights_[s][p] = pdfs_[s][p]->xfxQ(id1,x1,Q)/pdforig_->xfxQ(id1,x1,Q) * pdfs_[s][p]->xfxQ(id2,x2,Q)/pdforig_->xfxQ(id2,x2,Q);
		}
	}
}

void PDFuncertainty::fill_replicas(string dirname, string name, double val, double weight, URStreamer& streamer) {
  if(!use_evt_weights_) {
    SetupWeights(streamer);
		size_t counter=0;
    for(size_t s = 0 ; s < sets_.size() ; ++s) {
      for(size_t p = 0 ; p < pdfs_[s].size() ; ++p) {
        hists_[dirname][name][counter].fill(val, weight*weights_[s][p]);
				counter++;
      }
    }
    return;
  }

	const vector<Mcweight>& ws =  streamer.MCWeights();
	if(hists_.find(dirname) == hists_.end()) {
    Logger::log().fatal() << "You are asking to store PDFs replicas for the directory " << dirname << " which does not exist" << endl; 
    throw 42;		
	}
	if(hists_[dirname].find(name) == hists_[dirname].end()) {
    Logger::log().fatal() << "You are asking to store PDFs replicas for the plot " << name << " which does not exist" << endl; 
    throw 42;		
	}
	if(hists_[dirname][name].size() != ws.size()) {
    Logger::log().fatal() << "I got " << ws.size() << " pdf shifts, which is not what I expected! (" << hists_[dirname][name].size() << ")" << endl; 
    throw 42;
	}

	for(size_t h = 0 ; h < ws.size() ; ++h) {
		hists_[dirname][name][h].fill(val, weight*ws[h].weights()/ws[0].weights());
	}
}


void PDFuncertainty::fill_replicas2D(string dirname, string name, double xval, double yval, double weight, URStreamer& streamer) {
  if(!use_evt_weights_) {
    SetupWeights(streamer);
		size_t counter=0;
    for(size_t s = 0 ; s < sets_.size() ; ++s) {
      for(size_t p = 0 ; p < pdfs_[s].size() ; ++p) {
        hists_[dirname][name][counter].fill(xval, yval, weight*weights_[s][p]);
				counter++;
      }
    }
    return;
  }

	const vector<Mcweight>& ws =  streamer.MCWeights();
	if(hists_.find(dirname) == hists_.end()) {
    Logger::log().fatal() << "You are asking to store PDFs replicas for the directory " << dirname << " which does not exist" << endl; 
    throw 42;		
	}
	if(hists_[dirname].find(name) == hists_[dirname].end()) {
    Logger::log().fatal() << "You are asking to store PDFs replicas for the plot " << name << " which does not exist" << endl; 
    throw 42;		
	}
	if(hists_[dirname][name].size() != ws.size()) {
    Logger::log().fatal() << "I got " << ws.size() << " pdf shifts, which is not what I expected! (" << hists_[dirname][name].size() << ")" << endl; 
    throw 42;
	}

	for(size_t h = 0 ; h < ws.size() ; ++h) {
		hists_[dirname][name][h].fill(xval, yval, weight*ws[h].weights()/ws[0].weights());
	}
}
