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

void myfuncln(Int_t& npar, Double_t* deriv, Double_t& val, Double_t* par, Int_t flag);

class TTBarSolver
{
	private:
		TMinuit minuit_;
		std::shared_ptr<TH2D> WTmass_right_;
		std::shared_ptr<TH1D> BTag_right_; 
		std::shared_ptr<TH1D> N_right_; 
		std::shared_ptr<TH2D> WTmass_wrong_;
		std::shared_ptr<TH1D> BTag_wrong_; 
		std::shared_ptr<TH1D> N_wrong_; 

		TLorentzVector* bhad_=0;
		TLorentzVector* j1had_=0;
		TLorentzVector* j2had_=0;
		TLorentzVector* blep_=0;
		TLorentzVector* llep_=0;
		TLorentzVector* met_=0;

		TLorentzVector bhadT_;
		TLorentzVector j1hadT_;
		TLorentzVector j2hadT_;
		TLorentzVector blepT_;
		TLorentzVector llepT_;
		TLorentzVector metT_;
		double ubhad_;
		double uj1had_;
		double uj2had_;
		double ublep_;
		double ullep_;
		double umetx_;
		double umety_;
		double rhomet_;

		double nschi_;

		double btagtest_;
		double nstest_;
		double masstest_;
		double res_;

		bool USEBTAG_;
		bool USENS_;
		bool USEMASS_;
	public:
		double Test(double* par);
		static TTBarSolver* TTBS; 
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

		void Solve(Jet* bhad, Jet* j1had, Jet* j2had, Jet* blep, TLorentzVector* llep, IDMet* met);

		//extrem unlikely hypothesis will return a value >= 1E10
		double Res() const {return res_;}//final discriminant
		double NSRes() const {return nstest_;}//-log(l) of neutriosolver 
		double BTagRes() const {return btagtest_;} //-log(l) of btagging
		double MassRes() const {return masstest_;} //-log(l) of 2d-mass test

		double NSChi2() const {return nschi_;}//chi2 of neutrinosolver

		//improved objects: currently only usefull for neutrino
		TLorentzVector BHad() const {return bhadT_;}
		TLorentzVector Wja() const {return j1hadT_;}
		TLorentzVector Wjb() const {return j2hadT_;}
		TLorentzVector BLep() const {return blepT_;}
		TLorentzVector L() const {return llepT_;}
		TLorentzVector Nu() const {return metT_;}
};

#endif
