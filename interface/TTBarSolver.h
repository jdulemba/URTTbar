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
#include "JetMETCorrections/Modules/interface/JetResolution.h"

using namespace std;
using namespace TMath;

class IDMet;
class TDirectory;
class Permutation;

class TTBarSolver {
private:
	std::shared_ptr<TH2D> WTmass_right_;
	std::shared_ptr<TH1D> N_right_; 

    std::shared_ptr<TH2D> Rel_Delt_WTmass_right_; // replaces WTmass_right_ in mass discriminant


        // Joseph added for discriminants
    std::shared_ptr<TH2D> Max_Mjet_3J_right_;
    std::shared_ptr<TH2D> Max_Mjet_4J_right_;
	std::shared_ptr<TH1F> N_3J_right_; 
        //


	std::shared_ptr<TH1D> lep_b_ratio_right_; 
	std::shared_ptr<TH1D> wj2_b_ratio_right_; 
	std::shared_ptr<TH1D> wj1_btag_right_; 
	std::shared_ptr<TH1D> wj2_btag_right_; 
	std::shared_ptr<TH1D> wj1_qgtag_right_; 
	std::shared_ptr<TH1D> wj2_qgtag_right_; 
	// std::shared_ptr<TH1D> bj2_btag_right_; 

	bool USEBTAG_    	 = false;
	bool USENS_      	 = false;
	bool USEMASS_    	 = false;
	bool useptratios_	 = false;
	bool usewjetqgtag_   = false;

        // Joseph added for perm discriminant
    bool USEPERM_        = false;
        //

        // find jet resolution and scale factors
//    JME::JetResolution resolution_;
//    JME::JetResolutionScaleFactor jer_sf_;

	const double mtop_ = 173.;
	const double mw_ = 80.;	
public:
	double Test(double* par);
	TTBarSolver(bool active=true);
	~TTBarSolver();

  template <class T>
  std::shared_ptr<T> preproccess_histo(TDirectory* dir, std::string path, std::string newname) {
    T* original = (T*) dir->Get( path.c_str() );
    if(!original) {
      Logger::log().fatal() << "Could not get " << path << " from the file " << dir->GetName() << "!" << std::endl;
      throw 42;
    }
    std::shared_ptr<T> ptr((T*) original->Clone(newname.c_str()));
		ptr->Scale(1./ptr->Integral(), "width");
    return ptr;
  }

  template <class T>
  std::shared_ptr<T> preproccess_histo(TDirectory* dir, std::string path_right, std::string path_wrong, std::string newname) {
    T* right = (T*) dir->Get( path_right.c_str() );
    T* wrong = (T*) dir->Get( path_wrong.c_str() );		
    if(!right || !wrong) {
      Logger::log().fatal() << "Could not get " << path_right << " or " << path_wrong << " from the file " << dir->GetName() << "!" << std::endl;
      throw 42;
    }
		right->Scale(1./right->Integral(), "width");
		wrong->Scale(1./wrong->Integral(), "width");
		right->Divide(wrong);
    std::shared_ptr<T> ptr((T*) right->Clone(newname.c_str()));
    return ptr;
  }

	void Solve(Permutation &hyp, bool lazy=false);


        // Joseph added for perm discriminant
	void Solve_3J(Permutation &hyp, bool lazy=true);
	void Solve_4J(Permutation &hyp, bool lazy=true);
        //
};

#endif
