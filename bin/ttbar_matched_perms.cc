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

class ttbar_matched_perms : public AnalyzerBase
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
        ttbar_matched_perms(const std::string output_filename):
            AnalyzerBase("ttbar_matched_perms", output_filename),
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


            vector<string> event_types_lower = {"Merged", "Lost"};
            vector<string> event_treatments = {"TreatMerged", "TreatLost"};
            vector<string> disc_cats = {"LCut", "GCut"};

            for( auto evt_type : event_types_lower ){ // Merged/Lost
                book_disc_solution_type_plots( "Matched_Perm_Plots/"+evt_type );

                for( auto evt_treat : event_treatments ){ // Treat event as Merged or Lost

                    book_disc_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Discr");
                    book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen" );
                    book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Reconstruction" );
                    book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Resolution" );

                    for( auto disc_cat : disc_cats ){
                        book_gen_plots( "Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Gen/"+disc_cat );
                        book_reco_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Reconstruction/"+disc_cat );
                        book_reso_plots("Matched_Perm_Plots/"+evt_type+"/"+evt_treat+"/Resolution/"+disc_cat );
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

            return std::make_pair(reco_thad_cth, reco_tlep_cth);
        }

        //gen
        std::pair< double, double > gen_costh_tops( GenTTBar &ttbar ){

            hyp::TTbar gen_ttang(ttbar);
            auto gen_ttcm = gen_ttang.to_CM();
            double gen_thad_cth = gen_ttang.unit3D().Dot(gen_ttcm.thad().unit3D());
            double gen_tlep_cth = gen_ttang.unit3D().Dot(gen_ttcm.tlep().unit3D());

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
            //string btag_cat;
            //string acceptance;
            //string partial_merge;
            //string event_perm_treatement;

            if( mp.Merged_Event() ) event_perm_type = "Merged";
            if( mp.Lost_Event() ) event_perm_type = "Lost";


            //    /// classify based on gen-level info
            //if( ttbar.all_partons_in_acceptance(20, 2.4) ) acceptance = "In_Acceptance";
            //else acceptance = "Not_In_Acceptance";
            //if(ttbar.partial_hadronic_merged(0.4) ) partial_merge = "THad_Partial_Merge";
            //else partial_merge = "THad_Not_Partial_Merge";
            //            
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

            }


            //TreatLost
            if( !lost_bp.IsEmpty() ){

                    // don't classify plots based on anything
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

            }

        }// end of matched_perm_info()



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
    URDriver<ttbar_matched_perms> test;
    int thing = test.run();
    return thing;
}
