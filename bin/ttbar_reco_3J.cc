#include "TROOT.h"
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "Analyses/URTTbar/interface/DR_TTGenMatcher.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/GenObject.h"
#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include <algorithm>
#include <list>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include <boost/algorithm/string/predicate.hpp>
#include "TMath.h"
#include "TRandom3.h"
#include <set>
#include "Analyses/URTTbar/interface/IDMet.h"
#include "TUUID.h"
#include "Analyses/URTTbar/interface/systematics.h"
//#include "Analyses/URTTbar/interface/TTBarPlots.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"

#include "TTree.h"
//#include <string>

using namespace TMath;
using namespace std;

class ttbar_reco_3J : public AnalyzerBase
{

    public:
        //enum TTNaming_Merged {RIGHT, MERGE_SWAP, MERGE, WRONG};

    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        //	    double njets_ = 0;

        double disc_cut_ = 2.;

        //Initialize ttree vars
        //merged treat merged
        double M_TM_THad_Pt = -1e10;
        double M_TM_Massdisc = 1e10;
        double M_TM_NSdisc = -1e10;
        double M_TM_Totaldisc = -1e10;
        //double M_TM_Totaldisc_vs_THadPt[2];
        //merged treat lost
        double M_TL_THad_Pt = -1e10;
        double M_TL_Massdisc = -1e10;
        double M_TL_NSdisc = -1e10;
        double M_TL_Totaldisc = -1e10;
        //lost treat merged
        double L_TM_THad_Pt = -1e10;
        double L_TM_Massdisc = 1e10;
        double L_TM_NSdisc = -1e10;
        double L_TM_Totaldisc = -1e10;
        //double L_TM_Totaldisc_vs_THadPt[2];
        //lost treat lost
        double L_TL_THad_Pt = -1e10;
        double L_TL_Massdisc = -1e10;
        double L_TL_NSdisc = -1e10;
        double L_TL_Totaldisc = -1e10;

        // bins for hists
        int mass_bins_ = 100;
        int pt_bins_ = 50;
        int eta_bins_ = 50;
        int costh_bins_ = 50;

        double pt_min_ = 0.;
        double pt_max_ = 1000.;
        double eta_min_ = -2.5;
        double eta_max_ = 2.5;
        double costh_min_ = -1.;
        double costh_max_ = 1.;

        // mass disc
        int massdisc_bins_ = 120;
        double massdisc_min_ = -10.;
        double massdisc_max_ = 20.;

        // ns disc
        int nsdisc_bins_ = 80;
        double nsdisc_min_ = -5.;
        double nsdisc_max_ = 15.;

        // total disc
        int combdisc_bins_ = 160;
        double combdisc_min_ = -10.;
        double combdisc_max_ = 30.;

        //histograms
        map< string, map < string, map< string, map< string, map< string, RObject > > > > > histos_merged_;
        map< string, map < string, map< string, map< string, map< string, map< string, RObject > > > > > > histos_disc_;
        //unordered_map<string, map< string, RObject> > histos_;
        //unordered_map<string, map< string, RObject> > histos_;
        map<string, map< string, RObject> > histos_;

        CutFlowTracker tracker_; //tracks how many events pass cuts

        //switches
        bool isData_, isTTbar_;

        //selectors and helpers
        TTObjectSelector object_selector_; //selects ttbar objects
        TTGenParticleSelector genp_selector_; //selects generator level objects
        TTBarSolver solver_; //solves ttbar events
        TTGenMatcher matcher_; //matches particles on generator level
        TTPermutator permutator_;

        DR_TTGenMatcher dr_matcher_;

        float evt_weight_;
        MCWeightProducer mc_weights_;
        IDJet::BTag cut_tight_b_ = IDJet::BTag::CSVTIGHT;
        IDJet::BTag cut_medium_b_ = IDJet::BTag::CSVMEDIUM;
        IDJet::BTag cut_loose_b_ = IDJet::BTag::CSVLOOSE;


        // TTrees for multidim analysis
        TTree *Merged_TreatMerged = 0;
        TTree *Merged_TreatLost = 0;
        TTree *Lost_TreatMerged = 0;
        TTree *Lost_TreatLost = 0;

    public:
        ttbar_reco_3J(const std::string output_filename):
            AnalyzerBase("ttbar_reco_3J", output_filename),
            tracker_(),
            object_selector_(),
            genp_selector_(TTGenParticleSelector::SelMode::LHE), //for ttJets files
            //genp_selector_(), // for other files
            solver_(),
            evt_weight_(1.),
            mc_weights_()
    {
        //set tracker
        tracker_.use_weight(&evt_weight_);
        //find out which sample we're running on

        opts::variables_map &values = URParser::instance().values();
        string output_file = values["output"].as<std::string>();
        string sample = systematics::get_sample(output_file);
        bool isSignal = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
        isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal;

        isData_  = boost::starts_with(sample, "data");
        Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << ", isSignal: " << isSignal << endl;
        if(isData_) {
            if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
            else object_selector_.lepton_type(1);
        }
        if(!isData_) mc_weights_.init(sample);
    };

        TDirectory* getDir(string path){
            TDirectory* out = (TDirectory*) outFile_.Get(path.c_str());
            if(out) return out;
            outFile_.mkdir(path.c_str());
            return (TDirectory*) outFile_.Get(path.c_str());
        }

        template <class H, typename ... Args>
            void book(string folder, string name, Args ... args)
            {
                getDir(folder)->cd();
                histos_[folder][name] = RObject::book<H>(name.c_str(), args ...);
            }

        //This method is called once per job at the beginning of the analysis
        //book here your histograms/tree and run every initialization needed

        /// book and fill plots at gen-level
        void book_gen_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", 100, 200., 2000.);// nbins, mass_min, mass_max
            book<TH1F>(folder+"/Mass", "THad", "", mass_bins_, 150., 200.);
            book<TH1F>(folder+"/Mass", "TLep", "", mass_bins_, 150., 200.);
            book<TH1F>(folder+"/Costh", "THad", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(folder+"/Costh", "TLep", "", costh_bins_, costh_min_, costh_max_);
            book<TH2D>(folder+"/Costh", "Mmerged_vs_Costh", "", costh_bins_, costh_min_, costh_max_, mass_bins_, 0., 300.);
            book<TH2D>(folder+"/Costh", "Ptmerged_vs_Costh", "", costh_bins_, costh_min_, costh_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(folder+"/Costh", "Etamerged_vs_Costh", "", costh_bins_, costh_min_, costh_max_, eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Pt", "THad_PT_P", "", pt_bins_, -0.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad_PZ_P", "", pt_bins_, -1.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TLep", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TTbar", "",pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Eta", "THad", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Eta", "TLep", "", eta_bins_, eta_min_, eta_max_);          
            //book<TH1F>(folder+"/Eta", "TTbar", "", 100, -2.5, 2.5);

            book<TH1F>(folder+"/BHad", "Mass", "", mass_bins_, 0., 50.);
            book<TH1F>(folder+"/BHad", "Pt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/BHad", "Eta", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/BLep", "Mass", "", mass_bins_, 0., 50.);
            book<TH1F>(folder+"/BLep", "Pt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/BLep", "Eta", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/WJa", "Mass", "", mass_bins_, 0., 50.);
            book<TH1F>(folder+"/WJa", "Pt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/WJa", "Eta", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/WJb", "Mass", "", mass_bins_, 0., 50.);
            book<TH1F>(folder+"/WJb", "Pt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/WJb", "Eta", "", eta_bins_, eta_min_, eta_max_);

        }
        void fill_gen_plots( string folder, GenTTBar &ttbar ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            auto bhad_dir = histos_.find(folder+"/BHad");
            auto blep_dir = histos_.find(folder+"/BLep");
            auto wja_dir = histos_.find(folder+"/WJa");
            auto wjb_dir = histos_.find(folder+"/WJb");

            mass_dir->second["TTbar"].fill(ttbar.M());
            pt_dir->second["TTbar"].fill(ttbar.Pt());
            //eta_dir->second["TTbar"].fill(ttbar.Eta());

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            if( ttbar.had_top() ){
                double merged_parton_mass = -10.;
                double merged_parton_pt = -10.;
                double merged_parton_eta = -10.;
                if( ttbar.only_merged_bhadwja(0.4) ){
                    merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->first).M();
                    merged_parton_pt = (*ttbar.had_b() + *ttbar.had_W()->first).Pt();
                    merged_parton_eta = (*ttbar.had_b() + *ttbar.had_W()->first).Eta();
                }
                if( ttbar.only_merged_bhadwjb(0.4) ){
                    merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->second).M();
                    merged_parton_pt = (*ttbar.had_b() + *ttbar.had_W()->second).Pt();
                    merged_parton_eta = (*ttbar.had_b() + *ttbar.had_W()->second).Eta();
                }
                if( ttbar.only_merged_wjawjb(0.4) ){
                    merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->first).M();
                    merged_parton_pt = (*ttbar.had_b() + *ttbar.had_W()->first).Pt();
                    merged_parton_eta = (*ttbar.had_b() + *ttbar.had_W()->first).Eta();
                }

                mass_dir->second["THad"].fill(ttbar.had_top()->M());
                pt_dir->second["THad"].fill(ttbar.had_top()->Pt());
                pt_dir->second["THad_PT_P"].fill(ttbar.had_top()->Pt()/ttbar.had_top()->P());
                pt_dir->second["THad_PZ_P"].fill(ttbar.had_top()->Pz()/ttbar.had_top()->P());
                eta_dir->second["THad"].fill(ttbar.had_top()->Eta());
                costh_dir->second["THad"].fill(gen_cosths.first);
                costh_dir->second["Mmerged_vs_Costh"].fill(gen_cosths.first, merged_parton_mass);
                costh_dir->second["Ptmerged_vs_Costh"].fill(gen_cosths.first, merged_parton_pt);
                costh_dir->second["Etamerged_vs_Costh"].fill(gen_cosths.first, merged_parton_eta);

                bhad_dir->second["Mass"].fill(ttbar.had_b()->M());
                bhad_dir->second["Pt"].fill(ttbar.had_b()->Pt());
                bhad_dir->second["Eta"].fill(ttbar.had_b()->Eta());

                wja_dir->second["Mass"].fill(ttbar.had_W()->first->M());
                wja_dir->second["Pt"].fill(ttbar.had_W()->first->Pt());
                wja_dir->second["Eta"].fill(ttbar.had_W()->first->Eta());

                wjb_dir->second["Mass"].fill(ttbar.had_W()->second->M());
                wjb_dir->second["Pt"].fill(ttbar.had_W()->second->Pt());
                wjb_dir->second["Eta"].fill(ttbar.had_W()->second->Eta());

            }
            if( ttbar.lep_top() ){
                //cout << "tlep m: " << ttbar.lep_top()->M() << endl;
                mass_dir->second["TLep"].fill(ttbar.lep_top()->M());
                pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt());
                eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta());
                costh_dir->second["TLep"].fill(gen_cosths.second);

                blep_dir->second["Mass"].fill(ttbar.lep_b()->M());
                blep_dir->second["Pt"].fill(ttbar.lep_b()->Pt());
                blep_dir->second["Eta"].fill(ttbar.lep_b()->Eta());
            }
        }
        /// book and fill plots at gen-level

        /// book and fill plots at reco-level
        void book_reco_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", mass_bins_, 200., 2000.);
            book<TH1F>(folder+"/Mass", "THad", "", mass_bins_, 50., 250.);
            book<TH1F>(folder+"/Costh", "THad", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(folder+"/Costh", "TLep", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(folder+"/Pt", "THad_PT_P", "", pt_bins_, -0.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad_PZ_P", "", pt_bins_, -1.1, 1.1);
            book<TH1F>(folder+"/Pt", "THad", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TLep", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Pt", "TTbar", "",pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder+"/Eta", "THad", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Eta", "TLep", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(folder+"/Eta", "TTbar", "",eta_bins_, eta_min_, eta_max_);
        }
        void fill_reco_plots( string folder, Permutation &perm ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(perm.LVect().M());
            pt_dir->second["TTbar"].fill(perm.LVect().Pt());
            eta_dir->second["TTbar"].fill(perm.LVect().Eta());

            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(perm.THad().M());
            pt_dir->second["THad"].fill(perm.THad().Pt());
            pt_dir->second["THad_PT_P"].fill(perm.THad().Pt()/perm.THad().P());
            pt_dir->second["THad_PZ_P"].fill(perm.THad().Pz()/perm.THad().P());
            eta_dir->second["THad"].fill(perm.THad().Eta());
            costh_dir->second["THad"].fill(reco_cosths.first);
            pt_dir->second["TLep"].fill(perm.TLep().Pt());
            eta_dir->second["TLep"].fill(perm.TLep().Eta());
            costh_dir->second["TLep"].fill(reco_cosths.second);
        }
        /// book and fill plots at reco-level

        /// book and fill resolution plots
        void book_reso_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder+"/Mass", "THad", "", mass_bins_, -1000., 500.);
            book<TH1F>(folder+"/Costh", "THad", "", costh_bins_, -2*costh_max_, 2*costh_max_);
            book<TH1F>(folder+"/Costh", "TLep", "", costh_bins_, -2*costh_max_, 2*costh_max_);
            book<TH1F>(folder+"/Pt", "THad_PT_P", "", pt_bins_, -0.2, 1.1);
            book<TH1F>(folder+"/Pt", "THad_PZ_P", "", pt_bins_, -0.5, 1.1);
            book<TH1F>(folder+"/Pt", "THad", "", pt_bins_, -1*pt_max_, pt_max_);
            book<TH1F>(folder+"/Pt", "TLep", "", pt_bins_, -1*pt_max_, pt_max_);
            book<TH1F>(folder+"/Pt", "TTbar", "",pt_bins_, -1*pt_max_, pt_max_);
            book<TH1F>(folder+"/Eta", "THad", "", eta_bins_, -2*eta_max_, 2*eta_max_);
            book<TH1F>(folder+"/Eta", "TLep", "", eta_bins_, -2*eta_max_, 2*eta_max_);
            //book<TH1F>(folder+"/Eta", "TTbar", "",eta_bins_, -2*eta_max_, 2*eta_max_);
        }
        void fill_reso_plots( string folder, Permutation &perm, GenTTBar &ttbar ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(ttbar.M() - perm.LVect().M());
            pt_dir->second["TTbar"].fill(ttbar.Pt() - perm.LVect().Pt());
            //eta_dir->second["TTbar"].fill(ttbar.Eta() - perm.LVect().Eta());

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(ttbar.had_top()->M() - perm.THad().M());
            pt_dir->second["THad"].fill(ttbar.had_top()->Pt() - perm.THad().Pt());
            pt_dir->second["THad_PT_P"].fill(ttbar.had_top()->Pt()/ttbar.had_top()->P() - perm.THad().Pt()/perm.THad().P());
            pt_dir->second["THad_PZ_P"].fill(ttbar.had_top()->Pz()/ttbar.had_top()->P() - perm.THad().Pz()/perm.THad().P());
            eta_dir->second["THad"].fill(ttbar.had_top()->Eta() - perm.THad().Eta());
            costh_dir->second["THad"].fill(gen_cosths.first - reco_cosths.first);
            pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt() - perm.TLep().Pt());
            eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta() - perm.TLep().Eta());
            costh_dir->second["TLep"].fill(gen_cosths.second - reco_cosths.second);
        }
        /// book and fill resolution plots

        /// book and fill kin plots for jets
        void book_jet_kin_plots( string folder ){
            book<TH1F>(folder, "Mass", "", mass_bins_, 0., 50.);
            book<TH1F>(folder, "Pt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(folder, "Eta", "", eta_bins_, eta_min_, eta_max_);
        }
        void fill_jet_kin_plots( string folder, IDJet* jet ){
            auto dir = histos_.find(folder);

            dir->second["Mass"].fill(jet->M());
            dir->second["Pt"].fill(jet->Pt());
            dir->second["Eta"].fill(jet->Eta());
        }
        /// book and fill kin plots for jets

        /// book and fill plots for permutation discriminants
        void book_disc_plots( string folder ){
            book<TH1F>(folder, "Massdisc", "", massdisc_bins_, massdisc_min_, massdisc_max_);
            book<TH1F>(folder, "NSdisc", "", nsdisc_bins_, nsdisc_min_, nsdisc_max_);
            book<TH1F>(folder, "Totaldisc", "", combdisc_bins_, combdisc_min_, combdisc_max_);
        }
        void fill_disc_plots( string folder, Permutation &perm ){
            auto dir = histos_.find(folder);

            dir->second["Massdisc"].fill(perm.MassDiscr());
            dir->second["NSdisc"].fill(perm.NuDiscr());
            dir->second["Totaldisc"].fill(perm.Prob());
        }
        /// book and fill plots for permutation discriminants

        /// book and fill plots for merged and lost permutation discriminants with different numbers of solutions
        void book_disc_solution_type_plots( string folder ){
            book<TH1F>(folder, "Only_Merged_BestPerm_Solution", "", combdisc_bins_, combdisc_min_, combdisc_max_);
            book<TH1F>(folder, "Only_Lost_BestPerm_Solution", "", combdisc_bins_, combdisc_min_, combdisc_max_);
            book<TH1F>(folder, "Neither_BestPerm_Solution", "", 1, 0.5, 1.5);
            book<TH2D>(folder, "Lost_BestPerm_Solution_vs_Merged_BestPerm_Solution", "", combdisc_bins_, combdisc_min_, combdisc_max_, combdisc_bins_, combdisc_min_, combdisc_max_);
        }
        void fill_disc_solution_type_plots( string folder, Permutation &merged_perm, Permutation &lost_perm ){
            auto dir = histos_.find(folder);

            if( !merged_perm.IsEmpty() && lost_perm.IsEmpty() ){ // only TreatMerged has solution
                dir->second["Only_Merged_BestPerm_Solution"].fill(merged_perm.Prob());
            }
            else if( merged_perm.IsEmpty() && !lost_perm.IsEmpty() ){ // only TreatLost has solution
                dir->second["Only_Lost_BestPerm_Solution"].fill(lost_perm.Prob());
            }
            else if( !merged_perm.IsEmpty() && !lost_perm.IsEmpty() ){ // both solutions exist
                dir->second["Lost_BestPerm_Solution_vs_Merged_BestPerm_Solution"].fill( merged_perm.Prob(), lost_perm.Prob() );
            }
            else{
                dir->second["Neither_BestPerm_Solution"].fill( 1. );
            }
        }
        /// book and fill plots for merged and lost permutation discriminants with different numbers of solutions


        /// book and fill plots for btagging
        void book_perm_btag( string folder ){
            book<TH1F>(folder, "BHad_BTag", "", 20, 0.0, 1.0);
            book<TH1F>(folder, "BLep_BTag", "", 20, 0.0, 1.0);
            book<TH1F>(folder, "WJa_BTag", "", 20, 0.0, 1.0);
            book<TH1F>(folder, "WJb_BTag", "", 20, 0.0, 1.0);
            book<TH1F>(folder, "num_BTag", "", 4, 0.0, 4.0);
        }
        void fill_perm_btag( string folder, Permutation &perm ){
            auto dir = histos_.find(folder);

            dir->second["BHad_BTag"].fill( perm.BHad()->csvIncl() );
            dir->second["BLep_BTag"].fill( perm.BLep()->csvIncl() );
            if( perm.WJa() ) dir->second["WJa_BTag"].fill( perm.WJa()->csvIncl() );
            if( perm.WJb() ) dir->second["WJb_BTag"].fill( perm.WJb()->csvIncl() );

            int n_btag = num_btag(perm);

            dir->second["num_BTag"].fill( n_btag );
        }
        /// book and fill plots for btagging

        /// book and fill plots for hadronic top mass corrections for lost-jet events
        void book_alpha_correction_plots( string folder ){
            book<TH2D>(folder, "Alpha_THad_P", "", 32, 0.9, 2.5, 500, 0., 50.);
            book<TH2D>(folder, "Alpha_THad_E", "", 32, 0.9, 2.5, 500, 0., 50.);

            // gen vs reco plots for bins of 173.1/reco M(thad)
                // Energy
            book<TH2D>(folder, "Gen_vs_Reco_THadE_0.9to1.1", "", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
            book<TH2D>(folder, "Gen_vs_Reco_THadE_1.1to1.3", "", 75, 0., 1500., 100, 0., 2000.); // 1.1 < 173.1/reco M(thad) < 1.3
            book<TH2D>(folder, "Gen_vs_Reco_THadE_1.3to1.5", "", 75, 0., 1500., 100, 0., 2000.); // 1.3 < 173.1/reco M(thad) < 1.5
            book<TH2D>(folder, "Gen_vs_Reco_THadE_1.5to1.7", "", 75, 0., 1500., 100, 0., 2000.); // 1.5 < 173.1/reco M(thad) < 1.7
            book<TH2D>(folder, "Gen_vs_Reco_THadE_1.7to1.9", "", 75, 0., 1500., 100, 0., 2000.); // 1.7 < 173.1/reco M(thad) < 1.9
            book<TH2D>(folder, "Gen_vs_Reco_THadE_1.9to2.1", "", 75, 0., 1500., 100, 0., 2000.); // 1.9 < 173.1/reco M(thad) < 2.1
            book<TH2D>(folder, "Gen_vs_Reco_THadE_2.1to2.3", "", 75, 0., 1500., 100, 0., 2000.); // 2.1 < 173.1/reco M(thad) < 2.3
            book<TH2D>(folder, "Gen_vs_Reco_THadE_2.3to2.5", "", 75, 0., 1500., 100, 0., 2000.); // 2.3 < 173.1/reco M(thad) < 2.5
               // Momentum
            book<TH2D>(folder, "Gen_vs_Reco_THadP_0.9to1.1", "", 75, 0., 1500., 100, 0., 2000.); // 0.9 < 173.1/reco M(thad) < 1.1
            book<TH2D>(folder, "Gen_vs_Reco_THadP_1.1to1.3", "", 75, 0., 1500., 100, 0., 2000.); // 1.1 < 173.1/reco M(thad) < 1.3
            book<TH2D>(folder, "Gen_vs_Reco_THadP_1.3to1.5", "", 75, 0., 1500., 100, 0., 2000.); // 1.3 < 173.1/reco M(thad) < 1.5
            book<TH2D>(folder, "Gen_vs_Reco_THadP_1.5to1.7", "", 75, 0., 1500., 100, 0., 2000.); // 1.5 < 173.1/reco M(thad) < 1.7
            book<TH2D>(folder, "Gen_vs_Reco_THadP_1.7to1.9", "", 75, 0., 1500., 100, 0., 2000.); // 1.7 < 173.1/reco M(thad) < 1.9
            book<TH2D>(folder, "Gen_vs_Reco_THadP_1.9to2.1", "", 75, 0., 1500., 100, 0., 2000.); // 1.9 < 173.1/reco M(thad) < 2.1
            book<TH2D>(folder, "Gen_vs_Reco_THadP_2.1to2.3", "", 75, 0., 1500., 100, 0., 2000.); // 2.1 < 173.1/reco M(thad) < 2.3
            book<TH2D>(folder, "Gen_vs_Reco_THadP_2.3to2.5", "", 75, 0., 1500., 100, 0., 2000.); // 2.3 < 173.1/reco M(thad) < 2.5


            // reco corrected and uncorrected plots
            book<TH1F>(folder, "Reco_MTHad", "", 50, 50., 250.);
            book<TH1F>(folder, "Reco_MTHad_Corrected", "", 100, 50., 450.);
            book<TH1F>(folder, "Reco_MTTbar", "", mass_bins_, 200., 2000.);
            book<TH1F>(folder, "Reco_MTTbar_Corrected", "", mass_bins_, 200., 2000.);
            book<TH2D>(folder, "Reco_vs_Gen_MTTbar", ";Gen M(ttbar); Reco M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);
            book<TH2D>(folder, "Reco_vs_Gen_Corrected_MTTbar", ";Gen M(ttbar); Reco Corrected M(ttbar)", mass_bins_, 200., 2000., mass_bins_, 200., 2000.);

            book<TH1F>(folder, "Reco_THad_Labframe_Costh", "", costh_bins_, -1., 1.);
            book<TH1F>(folder, "Reco_THad_Labframe_Costh_Corrected", "", costh_bins_, -1., 1.);
            book<TH1F>(folder, "Reco_THad_Costh", "", costh_bins_, -1., 1.);
            book<TH1F>(folder, "Reco_THad_Costh_Corrected", "", costh_bins_, -1., 1.);
            book<TH2D>(folder, "Reco_vs Gen_Corrected_THad_Costh", ";Gen Cos(theta)(t_{h}); Reco Corrected Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);
            book<TH2D>(folder, "Reco_vs Gen_THad_Costh", ";Gen Cos(theta)(t_{h}); Reco Cos(theta)(t_{h})", costh_bins_, -1., 1., costh_bins_, -1., 1.);

            // reso corrected and uncorrected plots
            book<TH1F>(folder, "Reso_MTHad", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTHad_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_2Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_vs_Gen_MTTbar_2Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_3Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_vs_Gen_MTTbar_3Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_4Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_vs_Gen_MTTbar_4Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_Corrected", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_Corrected_2Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_vs_Gen_MTTbar_2Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_Corrected_3Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_vs_Gen_MTTbar_3Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Reso_MTTbar_Corrected_4Partons", "", mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_vs_Gen_MTTbar_4Partons", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH1F>(folder, "Frac_Reso_MTHad", "", mass_bins_, -2., 2.);
            book<TH1F>(folder, "Frac_Reso_MTHad_Corrected", "", mass_bins_, -2., 2.);
            book<TH1F>(folder, "Frac_Reso_MTTbar", "", mass_bins_, -2., 2.);
            book<TH1F>(folder, "Frac_Reso_MTTbar_Corrected", "", mass_bins_, -2., 2.);

            book<TH1F>(folder, "Reso_THad_Costh", "", costh_bins_, -2., 2.);
            book<TH1F>(folder, "Reso_THad_Costh_Corrected", "", costh_bins_, -2., 2.);

            
            book<TH2D>(folder, "Reso_MTTbar_Corrected_LCut_THadPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_GCut_THadPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_LCut_THadMass_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_GCut_THadMass_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_LCut_TLepPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_GCut_TLepPt_vs_Gen_MTTbar", "", 80, 0., 2000., 50, 0., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_LCut_NuPz_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);
            book<TH2D>(folder, "Reso_MTTbar_Corrected_GCut_NuPz_vs_Gen_MTTbar", "", 80, 0., 2000., mass_bins_, -1000., 1000.);

        }
        void fill_alpha_correction_plots( string folder, GenTTBar &ttbar, Permutation &perm ){
            auto dir = histos_.find(folder);

            if( perm.THad().M() > 180.0 ) return;
            //dir->second["Reco_vs_Gen_MTTbar"].fill( ttbar.M(), perm.LVect().M() );

            dir->second["Alpha_THad_P"].fill( 173.1/perm.THad().M(), ttbar.had_top()->P()/perm.THad().P() );//mthad taken from pdg
            dir->second["Alpha_THad_E"].fill( 173.1/perm.THad().M(), ttbar.had_top()->E()/perm.THad().E() );

            if( 0.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.1 ){
                dir->second["Gen_vs_Reco_THadE_0.9to1.1"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_0.9to1.1"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.3 ){
                dir->second["Gen_vs_Reco_THadE_1.1to1.3"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_1.1to1.3"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.5 ){
                dir->second["Gen_vs_Reco_THadE_1.3to1.5"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_1.3to1.5"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.5 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.7 ){
                dir->second["Gen_vs_Reco_THadE_1.5to1.7"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_1.5to1.7"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.7 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 1.9 ){
                dir->second["Gen_vs_Reco_THadE_1.7to1.9"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_1.7to1.9"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 1.9 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.1 ){
                dir->second["Gen_vs_Reco_THadE_1.9to2.1"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_1.9to2.1"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 2.1 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.3 ){
                dir->second["Gen_vs_Reco_THadE_2.1to2.3"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_2.1to2.3"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }
            if( 2.3 <= 173.1/perm.THad().M() && 173.1/perm.THad().M() < 2.5 ){
                dir->second["Gen_vs_Reco_THadE_2.3to2.5"].fill( perm.THad().E(), ttbar.had_top()->E() );
                dir->second["Gen_vs_Reco_THadP_2.3to2.5"].fill( perm.THad().P(), ttbar.had_top()->P() );
            }


                //alphas found by fitting Alpha_THad_P/E hists
                    //values taken from 1degree vals ttp://home.fnal.gov/~jdulemba/Plots/ttbar_reco_3J/2018/Compare_Lost_Merged_Jets/Full/ttJetsM0/3J_Event_Plots/Final_Reco/Clear_and_MassCut_Classes/Class_Lost/Alpha_Correction/fit_parameters.json
            double alpha_E = 0.6356*( 173.1/perm.THad().M() ) + 0.2862; // only alpha_E used because it's more consistent over mtt spectrum 
            //double alpha_P = 0.3853*( 173.1/perm.THad().M() ) + 0.5502;

            //TLorentzVector Alpha_THad(alpha_P*perm.THad().Px(), alpha_P*perm.THad().Py(), alpha_P*perm.THad().Pz(), alpha_E*perm.THad().E());
            TLorentzVector Alpha_THad(alpha_E*perm.THad().Px(), alpha_E*perm.THad().Py(), alpha_E*perm.THad().Pz(), alpha_E*perm.THad().E());


            // costheta variables
            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
                // perm thad
            hyp::TTbar reco_ttang(perm); // doesn't work because not all wjets defined
            auto reco_ttcm = reco_ttang.to_CM();
            double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
            double reco_thad_labframe_cth = reco_ttang.unit3D().Dot(perm.THad().Vect().Unit());
            //cout << "reco LF cth: " << reco_thad_labframe_cth << endl;

                // corrected thad
            TLorentzVector reco_corr_ttang = Alpha_THad+perm.TLep();
            TLorentzVector reco_ttcm_corr_thad = Alpha_THad;
            reco_ttcm_corr_thad.Boost(-1*reco_corr_ttang.BoostVector());
            double reco_corr_thad_cth = reco_corr_ttang.Vect().Unit().Dot(reco_ttcm_corr_thad.Vect().Unit());    
            double reco_corr_labframe_thad_cth = reco_corr_ttang.Vect().Unit().Dot(Alpha_THad.Vect().Unit());
            //cout << "LF corr cth: " << reco_corr_labframe_cth << endl;

                // reco plots
            dir->second["Reco_MTHad"].fill( perm.THad().M() );
            dir->second["Reco_MTHad_Corrected"].fill( Alpha_THad.M() );

            dir->second["Reco_MTTbar"].fill( perm.LVect().M() );
            dir->second["Reco_MTTbar_Corrected"].fill( (Alpha_THad + perm.TLep()).M() );

            dir->second["Reco_vs_Gen_MTTbar"].fill( ttbar.M(), perm.LVect().M() );
            dir->second["Reco_vs_Gen_Corrected_MTTbar"].fill( ttbar.M(), (Alpha_THad + perm.TLep()).M() );

            dir->second["Reco_THad_Labframe_Costh"].fill( reco_thad_labframe_cth );
            dir->second["Reco_THad_Labframe_Costh_Corrected"].fill( reco_corr_labframe_thad_cth );
            dir->second["Reco_THad_Costh"].fill( reco_thad_cth );
            dir->second["Reco_THad_Costh_Corrected"].fill( reco_corr_thad_cth );
            dir->second["Reco_vs Gen_Corrected_THad_Costh"].fill( gen_cosths.first, reco_corr_thad_cth );
            dir->second["Reco_vs Gen_THad_Costh"].fill(gen_cosths.first, reco_thad_cth);
                // reso plots
            dir->second["Reso_MTHad"].fill( ttbar.had_top()->M() - perm.THad().M() );
            dir->second["Reso_MTHad_Corrected"].fill( ttbar.had_top()->M() - Alpha_THad.M() );

            dir->second["Reso_MTTbar"].fill( ttbar.M() - perm.LVect().M() );
            dir->second["Reso_MTTbar_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

            dir->second["Reso_MTTbar_Corrected"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );

            if( ttbar.two_partons_in_acceptance(20., 2.4) ){
                dir->second["Reso_MTTbar_2Partons"].fill( ttbar.M() - perm.LVect().M() );
                dir->second["Reso_MTTbar_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

                dir->second["Reso_MTTbar_Corrected_2Partons"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
                dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar_2Partons"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            }
            if( ttbar.three_partons_in_acceptance(20., 2.4) ){
                dir->second["Reso_MTTbar_3Partons"].fill( ttbar.M() - perm.LVect().M() );
                dir->second["Reso_MTTbar_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

                dir->second["Reso_MTTbar_Corrected_3Partons"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
                dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar_3Partons"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            }
            if( ttbar.all_partons_in_acceptance(20., 2.4) ){
                dir->second["Reso_MTTbar_4Partons"].fill( ttbar.M() - perm.LVect().M() );
                dir->second["Reso_MTTbar_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.M() - perm.LVect().M() );

                dir->second["Reso_MTTbar_Corrected_4Partons"].fill( ttbar.M() - (Alpha_THad + perm.TLep()).M() );
                dir->second["Reso_MTTbar_Corrected_vs_Gen_MTTbar_4Partons"].fill( ttbar.M(), ttbar.M() - (Alpha_THad + perm.TLep()).M() );
            }

            dir->second["Frac_Reso_MTHad"].fill( (ttbar.had_top()->M() - perm.THad().M())/ttbar.had_top()->M() );
            dir->second["Frac_Reso_MTHad_Corrected"].fill( (ttbar.had_top()->M() - Alpha_THad.M())/ttbar.had_top()->M() );
            dir->second["Frac_Reso_MTTbar"].fill( (ttbar.M() - perm.LVect().M())/ttbar.M() );
            dir->second["Frac_Reso_MTTbar_Corrected"].fill( (ttbar.M() - (Alpha_THad + perm.TLep()).M())/ttbar.M() );

            dir->second["Reso_THad_Costh"].fill( gen_cosths.first - reco_thad_cth );
            dir->second["Reso_THad_Costh_Corrected"].fill( gen_cosths.first - reco_corr_thad_cth );


            if( ttbar.M() - (Alpha_THad + perm.TLep()).M() > ttbar.M()/3. ){
                dir->second["Reso_MTTbar_Corrected_GCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt() );
                dir->second["Reso_MTTbar_Corrected_GCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M() );
                dir->second["Reso_MTTbar_Corrected_GCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt() );
                dir->second["Reso_MTTbar_Corrected_GCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz() );
            }
            else{
                dir->second["Reso_MTTbar_Corrected_LCut_THadPt_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.Pt() );
                dir->second["Reso_MTTbar_Corrected_LCut_THadMass_vs_Gen_MTTbar"].fill( ttbar.M(), Alpha_THad.M() );
                dir->second["Reso_MTTbar_Corrected_LCut_TLepPt_vs_Gen_MTTbar"].fill( ttbar.M(), perm.TLep().Pt() );
                dir->second["Reso_MTTbar_Corrected_LCut_NuPz_vs_Gen_MTTbar"].fill( ttbar.M(), perm.Nu().Pz() );
            }

        }
        /// book and fill plots for hadronic top mass corrections for lost-jet events


        virtual void begin()
        {
            Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
            outFile_.cd();


            opts::variables_map &values = URParser::instance().values();
            string output_file = values["output"].as<std::string>();
            string sample = systematics::get_sample(output_file);
            Logger::log().debug() << "		" << sample << endl;



            // TTrees for Lost and Merged events
            Merged_TreatMerged = new TTree("Merged_TreatMerged", "Tree with merged events using merged discr.");
            Merged_TreatMerged->Branch("M_TM_THad_Pt", &M_TM_THad_Pt, "M_TM_THad_Pt/D");
            Merged_TreatMerged->Branch("M_TM_Massdisc", &M_TM_Massdisc, "M_TM_Massdisc/D");
            Merged_TreatMerged->Branch("M_TM_NSdisc", &M_TM_NSdisc, "M_TM_NSdisc/D");
            Merged_TreatMerged->Branch("M_TM_Totaldisc", &M_TM_Totaldisc, "M_TM_Totaldisc/D");
            //Merged_TreatMerged->Branch("Totaldisc_vs_THadPt", &M_TM_Totaldisc_vs_THadPt, "M_TM_Totaldisc_vs_THadPt[2]/D");

            Merged_TreatLost = new TTree("Merged_TreatLost", "Tree with merged events using merged discr.");
            Merged_TreatLost->Branch("M_TL_THad_Pt", &M_TL_THad_Pt, "M_TL_THad_Pt/D");
            Merged_TreatLost->Branch("M_TL_Massdisc", &M_TL_Massdisc, "M_TL_Massdisc/D");
            Merged_TreatLost->Branch("M_TL_NSdisc", &M_TL_NSdisc, "M_TL_NSdisc/D");
            Merged_TreatLost->Branch("M_TL_Totaldisc", &M_TL_Totaldisc, "M_TL_Totaldisc/D");

            Lost_TreatMerged = new TTree("Lost_TreatMerged", "Tree with lost events using merged discr.");
            Lost_TreatMerged->Branch("L_TM_THad_Pt", &L_TM_THad_Pt, "L_TM_THad_Pt/D");
            Lost_TreatMerged->Branch("L_TM_Massdisc", &L_TM_Massdisc, "L_TM_Massdisc/D");
            Lost_TreatMerged->Branch("L_TM_NSdisc", &L_TM_NSdisc, "L_TM_NSdisc/D");
            Lost_TreatMerged->Branch("L_TM_Totaldisc", &L_TM_Totaldisc, "L_TM_Totaldisc/D");
            //Lost_TreatMerged->Branch("Totaldisc_vs_THadPt", &L_TM_Totaldisc_vs_THadPt, "_L_TM_Totaldisc_vs_THadPt[2]/D");

            Lost_TreatLost = new TTree("Lost_TreatLost", "Tree with lost events using merged discr.");
            Lost_TreatLost->Branch("L_TL_THad_Pt", &L_TL_THad_Pt, "L_TL_THad_Pt/D");
            Lost_TreatLost->Branch("L_TL_Massdisc", &L_TL_Massdisc, "L_TL_Massdisc/D");
            Lost_TreatLost->Branch("L_TL_NSdisc", &L_TL_NSdisc, "L_TL_NSdisc/D");
            Lost_TreatLost->Branch("L_TL_Totaldisc", &L_TL_Totaldisc, "L_TL_Totaldisc/D");


            //int nbins;
            //double mass_max, mass_min, nbins;


            //if( sample == "ttJetsM0" ){
            //    nbins = 40;
            //    ////			m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
            //    mass_min = 0.;
            //    mass_max = 2000.;
            //}


            //else if( sample == "ttJetsM700" ){
            //    nbins = 30;
            //    ////			m_bins[nbins+1] = {700, 800, 900, 1000};
            //    mass_min = 700.;
            //    mass_max = 1000.;
            //}
            //else if( sample == "ttJetsM1000" ){
            //    nbins = 10;
            //    ////			m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
            //    mass_min = 1000.;
            //    mass_max = 2000.;
            //}
            //else{
            //    nbins = 80;
            //    mass_min = 200.;
            //    mass_max = 1000.;
            //}

            string merged_expected = "Expected_Plots/MERGED";
            book<TH1F>(merged_expected, "Expected_Event_Categories_3J", "", 5, 0.5, 5.5);

            //book<TH1F>("Kin_Checks", "Lead_Jet_Pt", "", 20, 0., 1000.);
            //book<TH1F>("Kin_Checks", "Jet_Pt", "", 30, 0., 900.);

            vector<string> gen_folders;
            if( isTTbar_ ) gen_folders = {"FULLHAD", "SEMILEP", "FULLLEP"};
            else gen_folders = {""};

            for( auto event_type : gen_folders ){
                book_gen_plots( "Gen_Plots/"+event_type );
            }

            vector<string> event_types_lower = {"Merged", "Lost"};
            vector<string> event_treatments = {"TreatMerged", "TreatLost"};
            vector<string> disc_cats = {"LCut", "GCut"};
            vector<string> objects = {"BHad", "BLep", "WJa", "WJb"};

                // plots for each jet in event
            vector<string> jets = {"j1", "j2", "j3"};

                // cut flow categories
            vector<string> cut_flows = {"Before_Matching", "Bs_and_WJet", "3_Unique", "Merged_or_Lost", "Merged_Bs", "After_Matching", "After_Matching_merged_bp"};

                // classify events based on btagging
            vector<string> btag_cats = {"BTag_LCut", "BTag_GCut"};

                // classify events based on gen-level
            vector<string> gen_acceptances = {"In_Acceptance", "Not_In_Acceptance"};
            vector<string> partial_hadtop_merges = {"THad_Partial_Merge", "THad_Not_Partial_Merge"};

            //for( auto gen_acceptance : gen_acceptances ){
            //    for( auto partial_hadtop_merge : partial_hadtop_merges ){
            //        book_gen_plots( "Gen_Cat_Plots/"+gen_acceptance+"/"+partial_hadtop_merge+"/Gen" );
            //        book_gen_plots( "Gen_Cat_Plots/"+gen_acceptance+"/"+partial_hadtop_merge+"/TreatMerged/Gen" );
            //        book_reco_plots("Gen_Cat_Plots/"+gen_acceptance+"/"+partial_hadtop_merge+"/TreatMerged/Reconstruction" );
            //        book_reso_plots("Gen_Cat_Plots/"+gen_acceptance+"/"+partial_hadtop_merge+"/TreatMerged/Resolution" );

            //        for( auto cut_flow : cut_flows ){
            //            for( auto jet : jets ){
            //                book_jet_kin_plots("Gen_Cat_Plots/"+gen_acceptance+"/"+partial_hadtop_merge+"/Jet_Kin/"+cut_flow+"/Costh_Neg/"+jet );
            //                book_jet_kin_plots("Gen_Cat_Plots/"+gen_acceptance+"/"+partial_hadtop_merge+"/Jet_Kin/"+cut_flow+"/Costh_Pos/"+jet );
            //            }
            //        }
            //    }
            //}



            for( auto evt_type : event_types_lower ){ // Merged/Lost
                book_disc_solution_type_plots( "Matched_Perm_Plots/"+evt_type );
                //for( auto gen_acceptance : gen_acceptances ){
                //    for( auto partial_hadtop_merge : partial_hadtop_merges ){
                //        for( auto evt_treat : event_treatments ){ // Treat event as Merged or Lost
                //            book_perm_btag("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/BTag" );

                //            book_disc_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge );
                //            book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/Gen" );
                //            book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/Reconstruction" );
                //            book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/Resolution" );

                //            for( auto disc_cat : disc_cats ){
                //                book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/Gen/"+disc_cat );
                //                book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/Reconstruction/"+disc_cat );
                //                book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen_Reqs/"+gen_acceptance+"/"+partial_hadtop_merge+"/Resolution/"+disc_cat );
                //            }
                //        }
                //    }
                //}

                for( auto evt_treat : event_treatments ){ // Treat event as Merged or Lost
                    //book_perm_btag( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/BTag"  );

                    book_disc_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Discr");
                    book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen" );
                    book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Reconstruction" );
                    book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Resolution" );

                    for( auto disc_cat : disc_cats ){
                        book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen/"+disc_cat );
                        book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Reconstruction/"+disc_cat );
                        book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Resolution/"+disc_cat );
                    }

                    //for( auto object : objects ){
                    //    book_jet_kin_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+object );
                    //}

                    //for( auto btag_cat : btag_cats ){
                    //    book_perm_btag( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/BTag"  );

                    //    book_disc_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"" );
                    //    book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/Gen" );
                    //    book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/Reconstruction" );
                    //    book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/Resolution" );

                    //    for( auto disc_cat : disc_cats ){
                    //        book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/Gen/"+disc_cat );
                    //        book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/Reconstruction/"+disc_cat );
                    //        book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/"+btag_cat+"/Resolution/"+disc_cat );
                    //    }
                    //}
                }
            }


        //// plots for 3-jet events using event jets as perm
            vector<string> merged_evt_type_categories = {"RIGHT", "MERGED_SWAP", "MERGED", "WRONG"};
            vector<string> merged_best_perm_solutions = {"Clear_Classification/Class_Merged", "Final_Classification/Class_Merged",
                                                         "Only_Merged", "Merged_BP", "Both_BP/Merged_BP/Class_Merged",
                                                         "Both_BP/Merged_BP/Class_Lost", "Both_BP/Class_Merged",
                                                         "Both_BP/Opposite_Class/Merged_BP/Class_Merged",
                                                         "Both_BP/Opposite_Class/Merged_BP/Class_Lost"
                                                        };

            //// reco/resolution plots for events that have clear classifications
            //book_reco_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Reconstruction" );
            //book_reso_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Resolution" );

            //// final reco/resolution for all event classifications with cuts (based on thad mass)
            //book_reco_plots("3J_Event_Plots/Final_Classification/Class_Lost/Reconstruction" );
            //book_reso_plots("3J_Event_Plots/Final_Classification/Class_Lost/Resolution" );


                // MERGED 3-jet events
            for( auto merged_bp_solution : merged_best_perm_solutions ){

                    // plots for merged_bp from 3-jet events
                book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr" );
                book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen" );
                book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction" );
                book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution" );

                for( auto merged_evt_type_cat : merged_evt_type_categories ){
                        // plots for merged_bp from 3-jet events
                    book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr/"+merged_evt_type_cat );
                    book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen/"+merged_evt_type_cat );
                    book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction/"+merged_evt_type_cat );
                    book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution/"+merged_evt_type_cat );

                    for( auto disc_cat : disc_cats ){
                            // plots for merged_bp from 3-jet events
                        book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr/"+merged_evt_type_cat+"/"+disc_cat );
                        book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen/"+merged_evt_type_cat+"/"+disc_cat );
                        book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction/"+merged_evt_type_cat+"/"+disc_cat );
                        book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution/"+merged_evt_type_cat+"/"+disc_cat );

                        if( merged_bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Lost" || merged_bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Merged" ){
                            book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr/THad_Mass/"+disc_cat+"/"+merged_evt_type_cat );
                            book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen/THad_Mass/"+disc_cat+"/"+merged_evt_type_cat );
                            book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction/THad_Mass/"+disc_cat+"/"+merged_evt_type_cat );
                            book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution/THad_Mass/"+disc_cat+"/"+merged_evt_type_cat );
                            //cout << "3J_Event_Plots/" << merged_bp_solution << "/Resolution/THad_Mass/"+disc_cat+"/"+merged_evt_type_cat
                            //cout << "3J_Event_Plots/" << merged_bp_solution << "/Resolution/THad_Mass/" << disc_cat << "/" << merged_evt_type_cat << endl;
                        }

                    }
                }
                for( auto disc_cat : disc_cats ){
                        // plots for merged_bp from 3-jet events
                    book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen/"+disc_cat );
                    book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction/"+disc_cat );
                    book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution/"+disc_cat );

                    if( merged_bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Lost" || merged_bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Merged" ){
                        book_disc_plots("3J_Event_Plots/"+merged_bp_solution+"/Discr/THad_Mass/"+disc_cat );
                        book_gen_plots( "3J_Event_Plots/"+merged_bp_solution+"/Gen/THad_Mass/"+disc_cat );
                        book_reco_plots("3J_Event_Plots/"+merged_bp_solution+"/Reconstruction/THad_Mass/"+disc_cat );
                        book_reso_plots("3J_Event_Plots/"+merged_bp_solution+"/Resolution/THad_Mass/"+disc_cat );
                        //cout << "3J_Event_Plots/" << merged_bp_solution << "/Resolution/THad_Mass/"+disc_cat
                        //cout << "3J_Event_Plots/" << merged_bp_solution << "/Resolution/THad_Mass/" << disc_cat << endl;
                    }
                }
            }


                // LOST 3-jet events
            vector<string> lost_evt_type_categories = {"RIGHT", "LOST_SWAP", "LOST", "WRONG"};
            vector<string> lost_best_perm_solutions = { "Clear_Classification/Class_Lost", "Final_Classification/Class_Lost",
                                                        "Only_Lost", "Lost_BP", "Both_BP/Lost_BP/Class_Merged",
                                                        "Both_BP/Lost_BP/Class_Lost", "Both_BP/Class_Lost",
                                                        "Both_BP/Opposite_Class/Lost_BP/Class_Merged",
                                                        "Both_BP/Opposite_Class/Lost_BP/Class_Lost"
                                                      };

            for( auto lost_bp_solution : lost_best_perm_solutions ){

                    // plots for lost_bp from 3-jet events
                book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr" );
                book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen" );
                book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction" );
                book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution" );
                book_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Alpha_Correction" );

                for( auto lost_evt_type_cat : lost_evt_type_categories ){
                        // plots for lost_bp from 3-jet events
                    book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr/"+lost_evt_type_cat );
                    book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/"+lost_evt_type_cat );
                    book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/"+lost_evt_type_cat );
                    book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/"+lost_evt_type_cat );
                    book_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Alpha_Correction/"+lost_evt_type_cat );

                    for( auto disc_cat : disc_cats ){
                    //        // plots for lost_bp from 3-jet events
                    //    book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr/"+lost_evt_type_cat+"/"+disc_cat );
                    //    book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/"+lost_evt_type_cat+"/"+disc_cat );
                    //    book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/"+lost_evt_type_cat+"/"+disc_cat );
                    //    book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/"+lost_evt_type_cat+"/"+disc_cat );

                        if( lost_bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Merged" || lost_bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Lost" ){
                            book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr/THad_Mass/"+disc_cat+"/"+lost_evt_type_cat );
                            book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/THad_Mass/"+disc_cat+"/"+lost_evt_type_cat );
                            book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/THad_Mass/"+disc_cat+"/"+lost_evt_type_cat );
                            book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/THad_Mass/"+disc_cat+"/"+lost_evt_type_cat );
                            book_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Alpha_Correction/THad_Mass/"+disc_cat+"/"+lost_evt_type_cat );
                            //cout << "3J_Event_Plots/" << lost_bp_solution << "/Resolution/THad_Mass/"+disc_cat+"/"+lost_evt_type_cat
                        }

                    }
                }
                for( auto disc_cat : disc_cats ){
                //        // plots for lost_bp from 3-jet events
                //    book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/"+disc_cat );
                //    book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/"+disc_cat );
                //    book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/"+disc_cat );

                    if( lost_bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Merged" || lost_bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Lost" ){
                        book_disc_plots("3J_Event_Plots/"+lost_bp_solution+"/Discr/THad_Mass/"+disc_cat );
                        book_gen_plots( "3J_Event_Plots/"+lost_bp_solution+"/Gen/THad_Mass/"+disc_cat );
                        book_reco_plots("3J_Event_Plots/"+lost_bp_solution+"/Reconstruction/THad_Mass/"+disc_cat );
                        book_reso_plots("3J_Event_Plots/"+lost_bp_solution+"/Resolution/THad_Mass/"+disc_cat );
                        book_alpha_correction_plots("3J_Event_Plots/"+lost_bp_solution+"/Alpha_Correction/THad_Mass/"+disc_cat );
                        //cout << "3J_Event_Plots/" << lost_bp_solution << "/Resolution/THad_Mass/"+disc_cat
                        //cout << "3J_Event_Plots/" << lost_bp_solution << "/Resolution/THad_Mass/" << disc_cat << endl;
                    }
                }
            }



            Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }


        // find costheta values for thad and tlep
        //reco
        std::pair< double, double > reco_costh_tops( Permutation &perm ){

            hyp::TTbar reco_ttang(perm); // doesn't work because not all wjets defined
            auto reco_ttcm = reco_ttang.to_CM();
            double reco_thad_cth = reco_ttang.unit3D().Dot(reco_ttcm.thad().unit3D());
            double reco_tlep_cth = reco_ttang.unit3D().Dot(reco_ttcm.tlep().unit3D());

            //cout << "reco_ttcm thad: " << reco_ttcm.thad() << endl;
            //cout << "reco_ttcm tlep: " << reco_ttcm.tlep() << endl;
            //cout << "reco_ttcm: " << reco_ttcm << endl;
            //cout << "reco_thad_cth: " << reco_thad_cth << endl;
            //cout << "reco_tlep_cth: " << reco_tlep_cth << endl;

            return std::make_pair(reco_thad_cth, reco_tlep_cth);
        }

        //gen
        std::pair< double, double > gen_costh_tops( GenTTBar &ttbar ){

            hyp::TTbar gen_ttang(ttbar);
            auto gen_ttcm = gen_ttang.to_CM();
            double gen_thad_cth = gen_ttang.unit3D().Dot(gen_ttcm.thad().unit3D());
            double gen_tlep_cth = gen_ttang.unit3D().Dot(gen_ttcm.tlep().unit3D());

            //cout << "gen ttcm had: " << gen_ttcm.thad() << endl;
            //cout << "gen ttcm lep: " << gen_ttcm.tlep() << endl;
            //cout << "gen thad: " << *ttbar.had_top() << endl;
            //cout << "gen tlep: " << ttbar.lep_top() << endl;
            return std::make_pair(gen_thad_cth, gen_tlep_cth);
        }

        // find best perm assuming it's merged
        Permutation merged_best_perm( Permutation &perm ){

            Permutation merged_bp;
            double merged_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Merged(test_perm);
                if( test_perm.Prob() < merged_lowest_Totaldisc_3J ){
                    merged_lowest_Totaldisc_3J = test_perm.Prob();
                    merged_bp = test_perm;
                }
            }
            return merged_bp;
        }

        // find best perm assuming it's lost 
        Permutation lost_best_perm( Permutation &perm ){

            Permutation lost_bp;
            double lost_lowest_Totaldisc_3J = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Lost(test_perm);
                if( test_perm.Prob() < lost_lowest_Totaldisc_3J ){
                    lost_lowest_Totaldisc_3J = test_perm.Prob();
                    lost_bp = test_perm;
                }
            }
            return lost_bp;
        }

        // find if a parton has a jet matched to it
        IDJet* genp_match( const GenObject* gen, vector<IDJet*> &jets ){

            IDJet* match = 0;
            double dr = 1e100;
            for( auto jet : jets ){
                if( jet->DeltaR(*gen) > 0.4 ) continue;
                if( jet->DeltaR(*gen) < dr ){
                    dr = jet->DeltaR(*gen);
                    match = jet;
                }
            }
            //cout << "gen match: " << match << endl;
            return match;
        }


        // find number of objects that pass btagging
        int num_btag( Permutation &perm){
            int num_btag = 0;
            if( perm.BHad()->BTagId(IDJet::BTag::CSVMEDIUM) ) num_btag++; // if bhad passes btag wp
            if( perm.BLep() != perm.BHad() && perm.BLep()->BTagId(IDJet::BTag::CSVMEDIUM) ) num_btag++; // if blep passes btag wp and isn't same object as bhad
            if( perm.WJa() && perm.WJa() != perm.BHad() && perm.WJa() != perm.BLep() && perm.WJa()->BTagId(IDJet::BTag::CSVMEDIUM) ) num_btag++; // wja passes btag wp and isn't same as either b
            if( perm.WJb() && perm.WJb() != perm.BHad() && perm.WJb() != perm.BLep() && perm.WJb() != perm.WJa() && perm.WJb()->BTagId(IDJet::BTag::CSVMEDIUM) ) num_btag++;

            return num_btag;
        }

        // events with 3 jets
        std::pair< Permutation, Permutation >  process_all3J_evt(URStreamer &event)
            //std::pair< Permutation, Permutation >  process_all3J_evt(URStreamer &event, Permutation &matched_perm )
        {
            //auto merged_disc_tp_dir = histos_.find("Test_Perm_Disc_Plots/MERGED");
            //auto lost_disc_tp_dir = histos_.find("Test_Perm_Disc_Plots/LOST");
            //initialize permutation objects
            //jets
            IDJet* wj1 = 0;
            IDJet* wj2 = 0;
            IDJet* bj1 = 0;
            IDJet* bj2 = 0;
            Permutation empty_perm; // perm.WJa(), WJb(), BHad(), BLep()

            vector<IDJet*> jets_vector;
            for( auto jet : object_selector_.clean_jets() ){
                jets_vector.push_back(jet);
            }
            sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });

            tracker_.track("All Num BTag");

            if( !jets_vector[1]->BTagId(IDJet::BTag::CSVMEDIUM) ) return std::make_pair(empty_perm, empty_perm); // require at least 2 jets that pass btag

            tracker_.track("Num BTag >= 2");

            bj1 = jets_vector[0];
            bj2 = jets_vector[1];
            wj1 = jets_vector[2];

            // create perm from the event and find best perms from merged/lost TTBarSolver
            Permutation event_perm(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wjat, wjb, bhad, blep
            Permutation merged_bp = merged_best_perm( event_perm );
            Permutation lost_bp = lost_best_perm( event_perm );

            return std::make_pair(merged_bp, lost_bp);
        } // end of process_all3J_evt()


        void matched_perm_info( URStreamer &event ){
            permutator_.reset_3J();

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) return; // skip to next event if it's not semilep

            // get matched perm from event
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());

            if( mp.IsEmpty() || !(mp.BHad() && mp.BLep()) || !(mp.WJa() || mp.WJb()) ) return; //make sure matched_perm has BHad, BLep, and at least one WJet
            if( mp.unique_matches() < 3 ) return; // require perm to have at least 3 distinct objects
            if( !(mp.Merged_Event() || mp.Lost_Event()) ) return;// require event to be merged or lost
            if( mp.Merged_BHadBLep() ) return; // event must have separate bhad and blep

            /// event is merged or lost with two b's and a wjet after this point

            // initialize different event categories/classifications
            string event_perm_type;

            string disc_cut;
            string btag_cat;
            string acceptance;
            string partial_merge;
            //string event_perm_treatement;

            if( mp.Merged_Event() ) event_perm_type = "Merged";
            if( mp.Lost_Event() ) event_perm_type = "Lost";


                /// classify based on gen-level info
            if( ttbar.all_partons_in_acceptance(20, 2.4) ) acceptance = "In_Acceptance";
            else acceptance = "Not_In_Acceptance";
            if(ttbar.partial_hadronic_merged(0.4) ) partial_merge = "THad_Partial_Merge";
            else partial_merge = "THad_Not_Partial_Merge";
                        
            ///


            // reset vars for ttrees
            M_TM_THad_Pt = -1e10;
            M_TM_Massdisc = 1e10;
            M_TM_NSdisc = -1e10;
            M_TM_Totaldisc = -1e10;

            M_TL_THad_Pt = -1e10;
            M_TL_Massdisc = 1e10;
            M_TL_NSdisc = -1e10;
            M_TL_Totaldisc = -1e10;

            L_TM_THad_Pt = -1e10;
            L_TM_Massdisc = 1e10;
            L_TM_NSdisc = -1e10;
            L_TM_Totaldisc = -1e10;

            L_TL_THad_Pt = -1e10;
            L_TL_Massdisc = 1e10;
            L_TL_NSdisc = -1e10;
            L_TL_Totaldisc = -1e10;



            //// Find best permutations based on Merged and Lost solvers
            // merged perm
            Permutation merged_bp = merged_best_perm(mp);

            // lost perm
            Permutation lost_bp = lost_best_perm(mp);


            // fill plots
            if( !merged_bp.IsEmpty() && lost_bp.IsEmpty() ){ // only TreatMerged has solution
                //mp_dir[event_perm_type]["Merged_BestPerm_Solution"][""][""][""].fill( merged_bp.Prob() );
                if( mp.Merged_Event() ){
                    M_TM_THad_Pt = merged_bp.THad().Pt();
                    M_TM_Massdisc = merged_bp.MassDiscr();
                    M_TM_NSdisc = merged_bp.NuDiscr();
                    M_TM_Totaldisc = merged_bp.Prob();
                    Merged_TreatMerged->Fill();

                    //M_TL_THad_Pt = lost_bp.THad().Pt();
                    M_TL_Massdisc = lost_bp.MassDiscr();
                    M_TL_NSdisc = lost_bp.NuDiscr();
                    M_TL_Totaldisc = lost_bp.Prob();
                    Merged_TreatLost->Fill();
                }

                if( mp.Lost_Event() ){
                    L_TM_Massdisc = merged_bp.MassDiscr();
                    L_TM_NSdisc = merged_bp.NuDiscr();
                    L_TM_Totaldisc = merged_bp.Prob();
                    Lost_TreatMerged->Fill();

                    //L_TL_THad_Pt = lost_bp.THad().Pt();
                    L_TL_Massdisc = lost_bp.MassDiscr();
                    L_TL_NSdisc = lost_bp.NuDiscr();
                    L_TL_Totaldisc = lost_bp.Prob();
                    Lost_TreatLost->Fill();
                }
            }
            if( merged_bp.IsEmpty() && !lost_bp.IsEmpty() ){ // only TreatLost has solution
                if( mp.Merged_Event() ){
                    //M_TM_THad_Pt = merged_bp.THad().Pt();
                    M_TM_Massdisc = merged_bp.MassDiscr();
                    M_TM_NSdisc = merged_bp.NuDiscr();
                    M_TM_Totaldisc = merged_bp.Prob();
                    Merged_TreatMerged->Fill();

                    M_TL_THad_Pt = lost_bp.THad().Pt();
                    M_TL_Massdisc = lost_bp.MassDiscr();
                    M_TL_NSdisc = lost_bp.NuDiscr();
                    M_TL_Totaldisc = lost_bp.Prob();
                    Merged_TreatLost->Fill();
                }

                if( mp.Lost_Event() ){
                    //L_TM_THad_Pt = merged_bp.THad().Pt();
                    L_TM_Massdisc = merged_bp.MassDiscr();
                    L_TM_NSdisc = merged_bp.NuDiscr();
                    L_TM_Totaldisc = merged_bp.Prob();
                    Lost_TreatMerged->Fill();

                    L_TL_THad_Pt = lost_bp.THad().Pt();
                    L_TL_Massdisc = lost_bp.MassDiscr();
                    L_TL_NSdisc = lost_bp.NuDiscr();
                    L_TL_Totaldisc = lost_bp.Prob();
                    Lost_TreatLost->Fill();
                }
            }
            if( !merged_bp.IsEmpty() && !lost_bp.IsEmpty() ){ // both solutions exist
                if( mp.Merged_Event() ){
                    M_TM_THad_Pt = merged_bp.THad().Pt();
                    M_TM_Massdisc = merged_bp.MassDiscr();
                    M_TM_NSdisc = merged_bp.NuDiscr();
                    M_TM_Totaldisc = merged_bp.Prob();
                    Merged_TreatMerged->Fill();

                    M_TL_THad_Pt = lost_bp.THad().Pt();
                    M_TL_Massdisc = lost_bp.MassDiscr();
                    M_TL_NSdisc = lost_bp.NuDiscr();
                    M_TL_Totaldisc = lost_bp.Prob();
                    Merged_TreatLost->Fill();
                }

                if( mp.Lost_Event() ){
                    L_TM_THad_Pt = merged_bp.THad().Pt();
                    L_TM_Massdisc = merged_bp.MassDiscr();
                    L_TM_NSdisc = merged_bp.NuDiscr();
                    L_TM_Totaldisc = merged_bp.Prob();
                    Lost_TreatMerged->Fill();

                    L_TL_THad_Pt = lost_bp.THad().Pt();
                    L_TL_Massdisc = lost_bp.MassDiscr();
                    L_TL_NSdisc = lost_bp.NuDiscr();
                    L_TL_Totaldisc = lost_bp.Prob();
                    Lost_TreatLost->Fill();
                }
            }

            fill_disc_solution_type_plots( "Matched_Perm_Plots/"+event_perm_type, merged_bp, lost_bp );

            //TreatMerged
            if( !merged_bp.IsEmpty() ){

                // don't classify plots based on anything
                //fill_perm_btag( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/BTag", merged_bp);

                fill_disc_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Discr", merged_bp);
                fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen", ttbar);
                fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Reconstruction", merged_bp);
                fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Resolution", merged_bp, ttbar);

                    //fill plots based on Prob value
                if( merged_bp.Prob() < 2 ) disc_cut = "LCut";
                else disc_cut = "GCut";
                fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen/"+disc_cut, ttbar);
                fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Reconstruction/"+disc_cut, merged_bp);
                fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Resolution/"+disc_cut, merged_bp, ttbar);

                //    //fill kin plots for perm objects
                //fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/BHad", merged_bp.BHad());
                //fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/BLep", merged_bp.BLep());
                //if( merged_bp.WJa() ) fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/WJa", merged_bp.WJa());
                //if( merged_bp.WJb() ) fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/WJb", merged_bp.WJb());
                /////

                //// classify based on btagging
                //int n_btag = num_btag(merged_bp);
                //if( n_btag < 2 ) btag_cat = "BTag_LCut";
                //else btag_cat = "BTag_GCut";

                //fill_perm_btag( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/BTag", merged_bp);

                //fill_disc_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"", merged_bp);
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/Gen", ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/Reconstruction", merged_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/Resolution", merged_bp, ttbar);

                //    //fill plots based on Prob value
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/Gen/"+disc_cut, ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/Reconstruction/"+disc_cut, merged_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/"+btag_cat+"/Resolution/"+disc_cut, merged_bp, ttbar);
                /////
                
                ///// classify based on gen-level info
                //fill_perm_btag( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/BTag", merged_bp);

                //fill_disc_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"", merged_bp);
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Gen", ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Reconstruction", merged_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Resolution", merged_bp, ttbar);

                //    //fill plots based on Prob value
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Gen/"+disc_cut, ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Reconstruction/"+disc_cut, merged_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatMerged/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Resolution/"+disc_cut, merged_bp, ttbar);
                /////

            }


            //TreatLost
            if( !lost_bp.IsEmpty() ){

                //// don't classify plots based on anything
                //fill_perm_btag( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/BTag", lost_bp);

                fill_disc_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Discr", lost_bp);
                fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen", ttbar);
                fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Reconstruction", lost_bp);
                fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Resolution", lost_bp, ttbar);

                    //fill plots based on Prob value
                if( lost_bp.Prob() < 2 ) disc_cut = "LCut";
                else disc_cut = "GCut";
                fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen/"+disc_cut, ttbar);
                fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Reconstruction/"+disc_cut, lost_bp);
                fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Resolution/"+disc_cut, lost_bp, ttbar);

                //    //fill kin plots for perm objects
                //fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/BHad", lost_bp.BHad());
                //fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/BLep", lost_bp.BLep());
                //if( lost_bp.WJa() ) fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/WJa", lost_bp.WJa());
                //if( lost_bp.WJb() ) fill_jet_kin_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/WJb", lost_bp.WJb());
                /////

                //// classify based on btagging
                //int n_btag = num_btag(lost_bp);
                //if( n_btag < 2 ) btag_cat = "BTag_LCut";
                //else btag_cat = "BTag_GCut";

                //fill_perm_btag( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/BTag", lost_bp);

                //fill_disc_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"", lost_bp);
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/Gen", ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/Reconstruction", lost_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/Resolution", lost_bp, ttbar);

                //    //fill plots based on Prob value
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/Gen/"+disc_cut, ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/Reconstruction/"+disc_cut, lost_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/"+btag_cat+"/Resolution/"+disc_cut, lost_bp, ttbar);
                /////
                
                ///// classify based on gen-level info
                //fill_perm_btag( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/BTag", lost_bp);

                //fill_disc_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"", lost_bp);
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Gen", ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Reconstruction", lost_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Resolution", lost_bp, ttbar);

                //    //fill plots based on Prob value
                //fill_gen_plots( "Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Gen/"+disc_cut, ttbar);
                //fill_reco_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Reconstruction/"+disc_cut, lost_bp);
                //fill_reso_plots("Matched_Perm_Plots/"+event_perm_type+"/TreatLost/Gen_Reqs/"+acceptance+"/"+partial_merge+"/Resolution/"+disc_cut, lost_bp, ttbar);
                /////

            }

        }// end of matched_perm_info()



        void gen_cats( URStreamer &event ){ // create plots for perms that are based on gen-level classifications
            permutator_.reset_3J();

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) return; // skip to next event if it's not semilep
            tracker_.track("Gen_Cats/SEMILEP");


                /// classify based on gen-level info
            string acceptance;
            string partial_merge;
            if( ttbar.all_partons_in_acceptance(20, 2.4) ) acceptance = "In_Acceptance";
            else acceptance = "Not_In_Acceptance";
            if(ttbar.partial_hadronic_merged(0.4) ) partial_merge = "THad_Partial_Merge";
            else partial_merge = "THad_Not_Partial_Merge";
                
            tracker_.track("Gen_Cats/"+acceptance+"/"+partial_merge);        
            fill_gen_plots( "Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Gen", ttbar);

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            string gen_costh_cat;
            if( gen_cosths.first < 0 ) gen_costh_cat = "Costh_Neg";
            else gen_costh_cat = "Costh_Pos";

            //double merged_parton_mass = 0;
            //if( ttbar.only_merged_bhadwja(0.4) ) merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->first).M();
            //if( ttbar.only_merged_bhadwjb(0.4) ) merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->second).M();
            //if( ttbar.only_merged_wjawjb(0.4) )  merged_parton_mass = (*ttbar.had_b() + *ttbar.had_W()->first).M();
            //cout << "merged parton mass: " << merged_parton_mass << endl;

            vector<IDJet*> jets_vector;
            //const vector<Genjet>& genjets = event.genjets();
            //cout << "number of genjets: " << genjets.size() << endl;
            //int i = 1;
            for( auto jet : object_selector_.clean_jets() ){
            //    cout << "jet " << i << endl;
            //    for( auto genjet : genjets ){
            //        double dr = jet->DeltaR(genjet);
            //        cout << "   DR(jet, genjet): " << dr << endl;
            //        if( dr < 0.4 ) cout << "jet pt: " << jet->Pt() << ", genjet pt: " << genjet.Pt() << endl;
            //    }
                jets_vector.push_back(jet);
            //    ++i;
            }
            sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->Pt() > B->Pt() ); });
            IDJet* j1 = jets_vector[0];
            IDJet* j2 = jets_vector[1];
            IDJet* j3 = jets_vector[2];
            fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Before_Matching/"+gen_costh_cat+"/j1", j1);
            fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Before_Matching/"+gen_costh_cat+"/j2", j2);
            fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Before_Matching/"+gen_costh_cat+"/j3", j3);

            // get matched perm from event
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());

            if( mp.IsEmpty() || !(mp.BHad() && mp.BLep()) || !(mp.WJa() || mp.WJb()) ){ //make sure matched_perm has BHad, BLep, and at least one WJet
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Bs_and_WJet/"+gen_costh_cat+"/j1", j1);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Bs_and_WJet/"+gen_costh_cat+"/j2", j2);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Bs_and_WJet/"+gen_costh_cat+"/j3", j3);
            }
            if( mp.unique_matches() < 3 ){ // require perm to have at least 3 distinct objects
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/3_Unique/"+gen_costh_cat+"/j1", j1);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/3_Unique/"+gen_costh_cat+"/j2", j2);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/3_Unique/"+gen_costh_cat+"/j3", j3);
            }
            if( !(mp.Merged_Event() || mp.Lost_Event()) ){ // require event to be merged or lost
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Merged_or_Lost/"+gen_costh_cat+"/j1", j1);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Merged_or_Lost/"+gen_costh_cat+"/j2", j2);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Merged_or_Lost/"+gen_costh_cat+"/j3", j3);
            }
            if( mp.Merged_BHadBLep() ){ // event must have separate bhad and blep
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Merged_Bs/"+gen_costh_cat+"/j1", j1);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Merged_Bs/"+gen_costh_cat+"/j2", j2);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/Merged_Bs/"+gen_costh_cat+"/j3", j3);
            }
            if( mp.IsEmpty() || !(mp.BHad() && mp.BLep()) || !(mp.WJa() || mp.WJb()) ) return; //make sure matched_perm has BHad, BLep, and at least one WJet
            if( mp.unique_matches() < 3 ) return; // require perm to have at least 3 distinct objects
            if( !(mp.Merged_Event() || mp.Lost_Event()) ) return;// require event to be merged or lost
            if( mp.Merged_BHadBLep() ) return; // event must have separate bhad and blep
            //if( !(mp.BHad()->BTagId(IDJet::BTag::CSVMEDIUM) && mp.BLep()->BTagId(IDJet::BTag::CSVMEDIUM)) ) return; // bhad and blep pass btag wp

            /// event is merged or lost with two b's and a wjet after this point

            fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/After_Matching/"+gen_costh_cat+"/j1", j1);
            fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/After_Matching/"+gen_costh_cat+"/j2", j2);
            fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/After_Matching/"+gen_costh_cat+"/j3", j3);

            //// Find best permutations based on Merged and Lost solvers
                // merged perm
            Permutation merged_bp = merged_best_perm(mp);
            if( !merged_bp.IsEmpty() ){
                fill_gen_plots( "Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/TreatMerged/Gen", ttbar);
                fill_reco_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/TreatMerged/Reconstruction", merged_bp);
                fill_reso_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/TreatMerged/Resolution", merged_bp, ttbar);

                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/After_Matching_merged_bp/"+gen_costh_cat+"/j1", j1);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/After_Matching_merged_bp/"+gen_costh_cat+"/j2", j2);
                fill_jet_kin_plots("Gen_Cat_Plots/"+acceptance+"/"+partial_merge+"/Jet_Kin/After_Matching_merged_bp/"+gen_costh_cat+"/j3", j3);
            }

        } // end of gen_cats



        void best_perm_categories( URStreamer &event ){
            permutator_.reset_3J();

            //generator selection
            bool selection = genp_selector_.select(event);
            if( !selection ){
                Logger::log().debug() << "event has no gen selection " << endl;
                return;
            }
            GenTTBar &ttbar = genp_selector_.ttbar_system();

            tracker_.track("Before best perm"); 

            string best_perm_solution;

            // get best perm from event determined by likelihood disc
            std::pair< Permutation, Permutation > best_perms = process_all3J_evt(event); // < merged_bp, lost_bp >
            Permutation merged_bp = best_perms.first;
            Permutation lost_bp = best_perms.second;
            if( merged_bp.IsEmpty() && lost_bp.IsEmpty() ){ // neither best_perm has solution
                tracker_.track("best perm doesn't exist");
                best_perm_solution = "No_BP_Solutions";
                return; // skip to next event if perm is empty
            }
            else if( !merged_bp.IsEmpty() && lost_bp.IsEmpty() ){ // only merged best_perm has solution
                tracker_.track("only merged_bp exists");
                best_perm_solution = "Only_Merged";
                merged_bp_cats( ttbar, merged_bp, best_perm_solution );
                //return; // skip to next event if perm is empty
            }
            else if( merged_bp.IsEmpty() && !lost_bp.IsEmpty() ){ // only lost best_perm has solution
                tracker_.track("only lost_bp exists");
                best_perm_solution = "Only_Lost";
                lost_bp_cats( ttbar, lost_bp, best_perm_solution );
                //return; // skip to next event if perm is empty
            }
            else{ // both best_perms have solution
                tracker_.track("both bp exists");
                best_perm_solution = "Both_BP";
                both_bp_cats( ttbar, merged_bp, lost_bp, best_perm_solution );
                //return; // skip to next event if perm is empty
            }

            if( !merged_bp.IsEmpty() ) merged_bp_cats( ttbar, merged_bp, "Merged_BP" );
            if( !lost_bp.IsEmpty() ) lost_bp_cats( ttbar, lost_bp, "Lost_BP" );

            //cout << best_perm_solution << endl;
        } // end of best_perm_categories



        void both_bp_cats( GenTTBar &ttbar, Permutation &merged_bp, Permutation &lost_bp, string bp_solution ){

            //if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
            //    tracker_.track("Not semilep events");
            //    return;
            //}
            //tracker_.track("semilep");

            string Merged_BP_Class;
            string Lost_BP_Class;

            // values for classifications taken from slide 3 of https://indico.cern.ch/event/714185/contributions/2936743/attachments/1617897/2572311/2018_March15_Update.pdf
            double TreatMerged_equation = merged_bp.Prob()-0.008064*merged_bp.THad().Pt();
            double TreatLost_equation = lost_bp.Prob()+0.016544*lost_bp.THad().Pt();

            if( TreatMerged_equation <= -0.397701 ) Merged_BP_Class = "Merged_BP/Class_Merged";
            else Merged_BP_Class = "Merged_BP/Class_Lost";
            if( TreatLost_equation <= 15.700815 ) Lost_BP_Class = "Lost_BP/Class_Lost";
            else Lost_BP_Class = "Lost_BP/Class_Merged";

            merged_bp_cats( ttbar, merged_bp, "Both_BP/"+Merged_BP_Class );
            lost_bp_cats( ttbar, lost_bp, "Both_BP/"+Lost_BP_Class );

            // event is classified as merged for both discriminants
            if( Merged_BP_Class == "Merged_BP/Class_Merged" && Lost_BP_Class == "Lost_BP/Class_Merged" ) merged_bp_cats( ttbar, merged_bp, "Both_BP/Class_Merged" );

            // event is classified as lost for both discriminants
            else if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Lost" ) lost_bp_cats( ttbar, lost_bp, "Both_BP/Class_Lost" );

            // event is classified as lost for one discriminant and merged for the other
            else{
                merged_bp_cats( ttbar, merged_bp, "Both_BP/Opposite_Class/"+Merged_BP_Class );
                lost_bp_cats( ttbar, lost_bp, "Both_BP/Opposite_Class/"+Lost_BP_Class );

                string merged_bp_thad_mass_cut;
                string lost_bp_thad_mass_cut;
                double thad_mass_cut;
                if( Merged_BP_Class == "Merged_BP/Class_Lost" && Lost_BP_Class == "Lost_BP/Class_Merged" ) thad_mass_cut = 170.;
                else thad_mass_cut = 140.;

                if( merged_bp.THad().M() > thad_mass_cut ) merged_bp_thad_mass_cut = "/THad_Mass/GCut"; // expected to have better resolution
                else merged_bp_thad_mass_cut = "/THad_Mass/LCut"; // expected to have worse resolution
                if( lost_bp.THad().M() > thad_mass_cut ) lost_bp_thad_mass_cut = "/THad_Mass/GCut"; // expected to have worse resolution
                else lost_bp_thad_mass_cut = "/THad_Mass/LCut"; // expected to have better resolution

                merged_bp_cats( ttbar, merged_bp, "Both_BP/Opposite_Class/"+Merged_BP_Class+merged_bp_thad_mass_cut );
                lost_bp_cats( ttbar, lost_bp, "Both_BP/Opposite_Class/"+Lost_BP_Class+lost_bp_thad_mass_cut );

            }
            //cout << "BP cats: " << Merged_BP_Class << ", " << Lost_BP_Class << endl;

        }


        void merged_bp_cats( GenTTBar &ttbar, Permutation &merged_bp, string bp_solution ){

            string merged_perm_status;
            string merged_perm_disc_cut_status;
            if( merged_bp.Prob() < disc_cut_ ){
                merged_perm_disc_cut_status = "LCut";
            }
            else{
                merged_perm_disc_cut_status = "GCut";
            }

            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
                tracker_.track("Not semilep events");
                //merged_perm_status = "NOTSEMILEP";
                //cout << merged_perm_status << ": " << evt_idx_ << endl;
                return;
            }
            tracker_.track("semilep");

            // get matched perm from SEMILEP events
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());


            if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Merged_Event() ) merged_perm_status = "WRONG";
            if( mp.Merged_Event() ){
                if( mp.AreBsSame(merged_bp) && ( merged_bp.WJa() == mp.WJa() || merged_bp.WJa() == mp.WJb() ) ) merged_perm_status = "RIGHT";
                else if( mp.AreBsFlipped(merged_bp) && ( merged_bp.WJa() == mp.WJa() || merged_bp.WJa() == mp.WJb() ) ) merged_perm_status = "MERGED_SWAP";
                else merged_perm_status = "MERGED";
            }
            else merged_perm_status = "WRONG";

            // merged_bp plots

                // filling reco and reso plots in cases where the classification is obvious
            if( bp_solution == "Both_BP/Class_Merged" ){// || bp_solution == "Only_Merged" ){
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Discr", merged_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Merged/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Reconstruction", merged_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Resolution", merged_bp, ttbar );

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Discr/"+merged_perm_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Merged/Gen/"+merged_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Reconstruction/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Merged/Resolution/"+merged_perm_status, merged_bp, ttbar );
            }

                // filling reco and reso plots in for final classifications
            if( bp_solution == "Both_BP/Class_Merged" || bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Merged/THad_Mass/GCut" || bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Lost/THad_Mass/GCut" ){
                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Merged/Discr", merged_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Merged/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Merged/Reconstruction", merged_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Merged/Resolution", merged_bp, ttbar );

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Merged/Discr/"+merged_perm_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Merged/Gen/"+merged_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Merged/Reconstruction/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Merged/Resolution/"+merged_perm_status, merged_bp, ttbar );
            }

                // fill plots for opposite classification cases
            if( bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Lost/THad_Mass/LCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Reconstruction/THad_Mass/LCut/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Resolution/THad_Mass/LCut/"+merged_perm_status, merged_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Reconstruction/THad_Mass/LCut", merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Resolution/THad_Mass/LCut", merged_bp, ttbar );
            }
            else if( bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Lost/THad_Mass/GCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Reconstruction/THad_Mass/GCut/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Resolution/THad_Mass/GCut/"+merged_perm_status, merged_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Reconstruction/THad_Mass/GCut", merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Lost/Resolution/THad_Mass/GCut", merged_bp, ttbar );
            }
            else if( bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Merged/THad_Mass/LCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Reconstruction/THad_Mass/LCut/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Resolution/THad_Mass/LCut/"+merged_perm_status, merged_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Reconstruction/THad_Mass/LCut", merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Resolution/THad_Mass/LCut", merged_bp, ttbar );
            }
            else if( bp_solution == "Both_BP/Opposite_Class/Merged_BP/Class_Merged/THad_Mass/GCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Reconstruction/THad_Mass/GCut/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Resolution/THad_Mass/GCut/"+merged_perm_status, merged_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Reconstruction/THad_Mass/GCut", merged_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Merged_BP/Class_Merged/Resolution/THad_Mass/GCut", merged_bp, ttbar );
            }

            else{
                    // not breaking up by perm comparison
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr", merged_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction", merged_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution", merged_bp, ttbar );

                    // not breaking up by disc val
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr/"+merged_perm_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+merged_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+merged_perm_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+merged_perm_status, merged_bp, ttbar );

                    // splitting up based on disc val
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+merged_perm_disc_cut_status, ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+merged_perm_disc_cut_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+merged_perm_disc_cut_status, merged_bp, ttbar );

                    // splitting up based on perm comp and disc val
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr/"+merged_perm_status+"/"+merged_perm_disc_cut_status, merged_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+merged_perm_status+"/"+merged_perm_disc_cut_status, ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+merged_perm_status+"/"+merged_perm_disc_cut_status, merged_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+merged_perm_status+"/"+merged_perm_disc_cut_status, merged_bp, ttbar );
            }

        } // end of merged_bp_cats


        void lost_bp_cats( GenTTBar &ttbar, Permutation &lost_bp, string bp_solution ){

            string lost_perm_status;
            //string lost_perm_disc_cut_status;
            //if( lost_bp.Prob() < disc_cut_ ){
            //    lost_perm_disc_cut_status = "LCut";
            //}
            //else{
            //    lost_perm_disc_cut_status = "GCut";
            //}


            if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
                tracker_.track("Not semilep events");
                return;
            }
            tracker_.track("semilep");

            // get matched perm from event
            Permutation mp = dr_matcher_.dr_match(
                    genp_selector_.ttbar_final_system(),
                    object_selector_.clean_jets(),
                    object_selector_.lepton(),
                    object_selector_.met(),
                    object_selector_.lepton_charge());


            if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Lost_Event() ) lost_perm_status = "WRONG";
            if( mp.Lost_Event() ){
                if( mp.AreBsSame(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "RIGHT";
                else if( mp.AreBsFlipped(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "LOST_SWAP";
                else lost_perm_status = "LOST";
            }
            else lost_perm_status = "WRONG";

            // lost_bp plots

                // filling reco and reso plots in cases where the classification is obvious
            if( bp_solution == "Both_BP/Class_Lost" ){// || bp_solution == "Only_Lost" ){
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Discr", lost_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Lost/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Reconstruction", lost_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Resolution", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Alpha_Correction", ttbar, lost_bp);

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Discr/"+lost_perm_status, lost_bp );
                fill_gen_plots( "3J_Event_Plots/Clear_Classification/Class_Lost/Gen/"+lost_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Reconstruction/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Resolution/"+lost_perm_status, lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Clear_Classification/Class_Lost/Alpha_Correction/"+lost_perm_status, ttbar, lost_bp);
            }

                // filling reco and reso plots in for final classifications
            if( bp_solution == "Both_BP/Class_Lost" || bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Lost/THad_Mass/LCut" || bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Merged/THad_Mass/LCut" ){
                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Lost/Discr", lost_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Lost/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Lost/Reconstruction", lost_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Lost/Resolution", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Final_Classification/Class_Lost/Alpha_Correction", ttbar, lost_bp );

                    //break up by event categories
                fill_disc_plots("3J_Event_Plots/Final_Classification/Class_Lost/Discr/"+lost_perm_status, lost_bp );
                fill_gen_plots( "3J_Event_Plots/Final_Classification/Class_Lost/Gen/"+lost_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/Final_Classification/Class_Lost/Reconstruction/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Final_Classification/Class_Lost/Resolution/"+lost_perm_status, lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Final_Classification/Class_Lost/Alpha_Correction/"+lost_perm_status, ttbar, lost_bp );
            }

                // fill plots for opposite classification cases
            if( bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Merged/THad_Mass/LCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Reconstruction/THad_Mass/LCut/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Resolution/THad_Mass/LCut/"+lost_perm_status, lost_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Reconstruction/THad_Mass/LCut", lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Resolution/THad_Mass/LCut", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Alpha_Correction/THad_Mass/LCut", ttbar, lost_bp );
            }
            else if( bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Merged/THad_Mass/GCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Reconstruction/THad_Mass/GCut/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Resolution/THad_Mass/GCut/"+lost_perm_status, lost_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Reconstruction/THad_Mass/GCut", lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Resolution/THad_Mass/GCut", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Merged/Alpha_Correction/THad_Mass/GCut", ttbar, lost_bp );
            }
            else if( bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Lost/THad_Mass/LCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Reconstruction/THad_Mass/LCut/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Resolution/THad_Mass/LCut/"+lost_perm_status, lost_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Reconstruction/THad_Mass/LCut", lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Resolution/THad_Mass/LCut", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Alpha_Correction/THad_Mass/LCut", ttbar, lost_bp );
            }
            else if( bp_solution == "Both_BP/Opposite_Class/Lost_BP/Class_Lost/THad_Mass/GCut" ){
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Reconstruction/THad_Mass/GCut/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Resolution/THad_Mass/GCut/"+lost_perm_status, lost_bp, ttbar );
                fill_reco_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Reconstruction/THad_Mass/GCut", lost_bp );
                fill_reso_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Resolution/THad_Mass/GCut", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/Both_BP/Opposite_Class/Lost_BP/Class_Lost/Alpha_Correction/THad_Mass/GCut", ttbar, lost_bp );
            }

            else{
                    // not breaking up by comparison
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr", lost_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen", ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction", lost_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution", lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/"+bp_solution+"/Alpha_Correction", ttbar, lost_bp);


                    // breaking up by perm comp
                fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr/"+lost_perm_status, lost_bp );
                fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+lost_perm_status, ttbar );
                fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+lost_perm_status, lost_bp );
                fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+lost_perm_status, lost_bp, ttbar );
                fill_alpha_correction_plots("3J_Event_Plots/"+bp_solution+"/Alpha_Correction/"+lost_perm_status, ttbar, lost_bp );

                //    // splitting up based on disc val
                //fill_disc_plots("3J_Event_Plots/"+bp_solution+"/Discr/"+lost_perm_status+"/"+lost_perm_disc_cut_status, lost_bp );
                //fill_gen_plots( "3J_Event_Plots/"+bp_solution+"/Gen/"+lost_perm_status+"/"+lost_perm_disc_cut_status, ttbar );
                //fill_reco_plots("3J_Event_Plots/"+bp_solution+"/Reconstruction/"+lost_perm_status+"/"+lost_perm_disc_cut_status, lost_bp );
                //fill_reso_plots("3J_Event_Plots/"+bp_solution+"/Resolution/"+lost_perm_status+"/"+lost_perm_disc_cut_status, lost_bp, ttbar );
            }

        } // end of lost_bp_cats


        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            URStreamer event(tree_);

            while(event.next() /*&& evt_idx_ < 30000*/)
            {
                evt_idx_++;
                if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << endl;
                tracker_.track("All Events");


                //generator selection
                bool selection = genp_selector_.select(event);
                if( !selection ){
                    Logger::log().debug() << "event has no gen selection " << endl;
                    continue;
                }
                tracker_.track("gen selection");
                GenTTBar &ttbar = genp_selector_.ttbar_system();

                string gen_folder;
                if( ttbar.type == GenTTBar::DecayType::FULLHAD ) gen_folder = "FULLHAD";
                else if( ttbar.type == GenTTBar::DecayType::SEMILEP ) gen_folder = "SEMILEP";
                else if( ttbar.type == GenTTBar::DecayType::FULLLEP ) gen_folder = "FULLLEP";
                else gen_folder = "";
                //fill_gen_plots("Gen_Plots/"+gen_folder, ttbar);
                //if( ttbar.had_top() ) cout << "gen thad M: " << ttbar.had_top()->M() << endl;
                //if( ttbar.M() > 700 ) continue;

                int njets = 0;
                if( object_selector_.select(event) ) njets = object_selector_.clean_jets().size();
                if( njets < 3 ) continue;
                //            Logger::log().debug() << "Beginning event " << evt_idx_ << ", njets = " << njets << endl;
                tracker_.track("njet cuts");

                /// 3 jet events
                if( njets == 3 ){ 
                    tracker_.track("njets = 3");
                    matched_perm_info(event);
                    best_perm_categories(event);
                    //merged_bp_cats(event);
                    //lost_bp_cats(event);
                    //gen_cats(event);
                }

            } // end of event loop

            Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
        } // end of analyze()

        //this method is called at the end of the job, by default saves
        //every histogram/tree produced, override it if you need something more
        virtual void end()
        {
            outFile_.Write();
            tracker_.writeTo(outFile_);
            Logger::log().debug() << "End of end() " << evt_idx_ << endl;
        }

        //do you need command-line or cfg options? If so implement this
        //method to book the options you need. CLI parsing is provided
        //by AnalysisFW/interface/URParser.h and uses boost::program_options
        //look here for a quickstart tutorial:
        //http://www.boost.org/doc/libs/1_51_0/doc/html/program_options/tutorial.html

        static void setOptions()
        {
            URParser &parser = URParser::instance();
            opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
            opts.add_options()
                ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                ("report,s", opts::value<int>()->default_value(1), "report every");
        }
};

//make it executable
int main(int argc, char *argv[])
{
    URParser &parser = URParser::instance(argc, argv);
    URDriver<ttbar_reco_3J> test;
    int thing = test.run();
    return thing;
}
