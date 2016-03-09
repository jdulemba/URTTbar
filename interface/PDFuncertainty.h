#ifndef HPDFUNCERTAINTY
#define HPDFUNCERTAINTY
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include <TH1D.h>
#include <TH2D.h>

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
  std::map<std::string, std::map<std::string,std::vector< std::vector<TH1F*> > > > hist1d_;
  std::vector< std::vector<double> > weights_;
  double oldx1_;
  bool use_evt_weights_;
  int nweights_;
  void SetupWeights(URStreamer& streamer);

  std::map<std::string, TDirectory*> histdir_;
  std::map<std::string, std::vector<TH1D*> > Whist1d_;


public:
  PDFuncertainty(const std::string setorigname, int memorig, const std::vector<std::string>& setnames); //custom PDF settings
  PDFuncertainty(int nweights); //uses weights stored in MC
  ~PDFuncertainty();

  template <typename ... Args>
	void book_replicas(std::string dirname, std::string name, Args ... args) {
    size_t nsets = (use_evt_weights_) ? 1 : sets_.size();
    hist1d_[dirname][name].resize(nsets);
    for(size_t s = 0 ; s < nsets ; ++s) {
      std::string setname = (use_evt_weights_) ? "mcws" : sets_[s]->name();
      size_t npdfs = (use_evt_weights_) ? nweights_ : pdfs_[s].size();
      for(size_t p = 0 ; p < npdfs ; ++p) {
        std::stringstream hname;
        hname  <<  name << "_" << setname << "_pdf_" << p;
        hist1d_[dirname][name][s].push_back(
          new TH1F(hname.str().c_str(), hname.str().c_str(), args ...)
          );
      }
    }
  }

  void fill_replicas(std::string dirname, std::string name, double val, double weight, URStreamer& streamer);
};

#endif
