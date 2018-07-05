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


        // Joseph added 3J distributions (merged and lost) for discriminants
        // merged events
        std::shared_ptr<TH2D> Max_Mjet_3J_merged_right_;// mbpjet_vs_maxmjet_merged_right
        std::shared_ptr<TH1F> NSchi2_3J_merged_right_;// nusolver_chi2_merged_right
        std::shared_ptr<TH1F> NSdist_3J_merged_right_;// nusolver_dist_merged_right

        // lost events
        std::shared_ptr<TH1F> Mbpjet_3J_lost_right_;// mbpjet_lost_right
        std::shared_ptr<TH1F> NSchi2_3J_lost_right_;// nusolver_chi2_lost_right
        std::shared_ptr<TH1F> NSdist_3J_lost_right_;// nusolver_dist_lost_right
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

        // Joseph added for 3 jet events (merged and lost)
        bool USE3JMERGED_        = false;
        bool USE3JLOST_        = false;
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

        void Solve(Permutation &hyp, bool lazy=true);

        // Joseph added for 3 jet events (merged and lost) 
        void Solve_3J_Merged(Permutation &hyp, bool lazy=true);
        void Solve_3J_Lost(Permutation &hyp, bool lazy=true);
};

#endif
