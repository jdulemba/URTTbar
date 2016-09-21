#ifndef HTTBSOLVER
#define HTTBSOLVER
#include <TMatrixD.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <iostream>
#include <memory>
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"


using namespace std;
using namespace TMath;

class IDMet;
class TDirectory;
class Permutation;

class TTBarSolver {
private:
	std::shared_ptr<TH2D> WTmass_right_;
	std::shared_ptr<TH1D> N_right_; 

	// std::shared_ptr<TH1D> lep_b_ratio_right_; 
	// std::shared_ptr<TH1D> wj2_b_ratio_right_; 
	// std::shared_ptr<TH1D> wj1_btag_right_; 
	// std::shared_ptr<TH1D> wj2_btag_right_; 
	// std::shared_ptr<TH1D> wj1_qgtag_right_; 
	// std::shared_ptr<TH1D> wj2_qgtag_right_; 
	// std::shared_ptr<TH1D> bj2_btag_right_; 

	bool USEBTAG_;
	bool USENS_;
	bool USEMASS_;

	const double mtop_ = 173.;
	const double mw_ = 80.;	
public:
	double Test(double* par);
	TTBarSolver();
	TTBarSolver(bool dummy);
	~TTBarSolver();

  template <class T>
  std::shared_ptr<T> preproccess_histo(TDirectory* dir, std::string path, std::string newname) {
    T* original = (T*) dir->Get( path.c_str() );
    if(!original) {
      Logger::log().fatal() << "Could not get " << path << " from the file " << dir->GetName() << "!" << std::endl;
      throw 42;
    }
    std::shared_ptr<T> ptr((T*) original->Clone(newname.c_str()));
		ptr->Scale(1./ptr->Integral("width"));
    return ptr;
  }

	void Init(string filename, bool usebtag = true, bool usens = true, bool usemass = true);//provide root file with probability distribution, switches if btag and neutrino solver information should be used for final discriminant Res()
	void Init(TDirectory* dir=0, bool usebtag=false, bool usens=false, bool usemass=false, string wtname="mWhad_vs_mtophad", string bname="btag_value", string nuname="nusolver_chi2");

	void Solve(Permutation &hyp, bool lazy=false);
};

#endif
