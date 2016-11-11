#ifndef HPDFUNCERTAINTY
#define HPDFUNCERTAINTY
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "URAnalysis/AnalysisFW/interface/RObject.h"

#include <LHAPDF/LHAPDF.h>
#include <sstream>

class URStreamer;
class TDirectory;

class PDFuncertainty
{
private:
  LHAPDF::PDFSet* setorig_;
  LHAPDF::PDF* pdforig_;
  std::vector<LHAPDF::PDFSet*> sets_;
  std::vector< std::vector<LHAPDF::PDF*> > pdfs_;
  std::map<std::string, std::map<std::string, std::vector<RObject> > > hists_;
  std::vector< std::vector<double> > weights_;
  double oldx1_;
  bool use_evt_weights_;
  int nweights_;
  void SetupWeights(URStreamer& streamer);

public:
  PDFuncertainty(const std::string setorigname, int memorig, const std::vector<std::string>& setnames); //custom PDF settings
  PDFuncertainty(int nweights); //uses weights stored in MC
  ~PDFuncertainty();

  template <class H, typename ... Args>
	void book_replicas(std::string dirname, std::string name, Args ... args) {
    size_t nsets = (use_evt_weights_) ? 1 : sets_.size();
		size_t nhistos = (use_evt_weights_) ? nweights_ : 0;
		if(!use_evt_weights_) {
			for(size_t s = 0 ; s < nsets ; ++s) {
				nhistos += pdfs_[s].size();
			}
		}
    hists_[dirname][name];
		hists_[dirname][name].reserve(nhistos);
    for(size_t s = 0 ; s < nsets ; ++s) {
      std::string setname = (use_evt_weights_) ? "mcws" : sets_[s]->name();
      size_t npdfs = (use_evt_weights_) ? nweights_ : pdfs_[s].size();
      for(size_t p = 0 ; p < npdfs ; ++p) {
        std::stringstream hname;
        hname  <<  name << "_" << setname << "_pdf_" << p;
        hists_[dirname][name].push_back(
					RObject::book<H>(hname.str().c_str(), hname.str().c_str(), args ...)
          );
      }
    }
  }

	void book_replicas(std::string dirname, std::string name, RObject& stamp) {
    size_t nsets = (use_evt_weights_) ? 1 : sets_.size();
		size_t nhistos = (use_evt_weights_) ? nweights_ : 0;
		if(!use_evt_weights_) {
			for(size_t s = 0 ; s < nsets ; ++s) {
				nhistos += pdfs_[s].size();
			}
		}
    hists_[dirname][name];
		hists_[dirname][name].reserve(nhistos);
    for(size_t s = 0 ; s < nsets ; ++s) {
      std::string setname = (use_evt_weights_) ? "mcws" : sets_[s]->name();
      size_t npdfs = (use_evt_weights_) ? nweights_ : pdfs_[s].size();
      for(size_t p = 0 ; p < npdfs ; ++p) {
        std::stringstream hname;
        hname  <<  name << "_" << setname << "_pdf_" << p;
        hists_[dirname][name].push_back(
					stamp.clone(hname.str())
          );
      }
    }
  }

  void fill_replicas(std::string dirname, std::string name, double val, double weight, URStreamer& streamer);
  void fill_replicas2D(std::string dirname, std::string name, double xval, double yval, double weight, URStreamer& streamer);
};

#endif
