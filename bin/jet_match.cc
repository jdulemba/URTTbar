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
#include "Analyses/URTTbar/interface/TTPermutator.h"
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

#include "TROOT.h"
#include "TTree.h"
//#include <string>

using namespace TMath;
//using namespace std;

class jet_match : public AnalyzerBase
{
    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        double njets_ = 0;

        const char *DRnames_[2] = {"DRP4", "DRP8"};
        double DR_[2] = {0.4, 0.8};
        //            const char *DRnames_[4] = {"DRP4", "DRP5", "DRP6", "DRP8"};
        //	    double DR_[4] = {0.4, 0.5, 0.6, 0.8};
        //            const char *DRnames_[1] = {"DRP4"};
        //	    double DR_[1] = {0.4};

        //histograms
        unordered_map<string, map< string, RObject> > histos_;
        //map<TTNaming, string> naming_;
        CutFlowTracker tracker_; //tracks how many events pass cuts

        int nbins;
        double mass_max_, mass_min_;

        int mass_bins_ = 500;
        int pt_bins_ = 200;
        int DR_bins_ = 200;
        int eta_bins_ = 200;
        int costh_bins_ = 100;

        double DR_min_ = 0.;
        double DR_max_ = 8.;
        double pt_min_ = 0.;
        double pt_max_ = 2000.;
        double eta_min_ = -2.4;
        double eta_max_ = 2.4;
        double costh_min_ = -1.;
        double costh_max_ = 1.;


        //switches
        bool isData_, isTTbar_;

        //selectors and helpers
        TTObjectSelector object_selector_; //selects ttbar objects
        TTGenParticleSelector genp_selector_; //selects generator level objects
        TTBarSolver solver_; //solves ttbar events
        TTGenMatcher matcher_; //matches particles on generator level
        TTPermutator permutator_;

        float evt_weight_;
        MCWeightProducer mc_weights_;
        IDJet::BTag cut_tight_b_ = IDJet::BTag::CSVTIGHT;
        IDJet::BTag cut_medium_b_ = IDJet::BTag::CSVMEDIUM;
        IDJet::BTag cut_loose_b_ = IDJet::BTag::CSVLOOSE;

    public:
        jet_match(const std::string output_filename):
            AnalyzerBase("jet_match", output_filename),
            tracker_(),
            object_selector_(),
            genp_selector_(TTGenParticleSelector::SelMode::LHE),
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
        Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
        if(isData_) {
            if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
            else object_selector_.lepton_type(1);
        }
        if(!isData_) mc_weights_.init(sample);

        //Init Solver
        string filename = "prob_ttJets.permProbComputer.test.root";
        Logger::log().debug() << "solver file: " << filename << endl;
        TFile probfile(DataFile(filename).path().c_str());
        TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
        //			solver_.Init(td, false, true, true, true, true); //probfile, btag, nusolv,massdis,anghad,anglep

        //			cut_tight_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.tightb"));
        //			cut_loose_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.looseb"));
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
        virtual void begin()
        {
            Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
            outFile_.cd();


            opts::variables_map &values = URParser::instance().values();
            string output_file = values["output"].as<std::string>();
            string sample = systematics::get_sample(output_file);
            Logger::log().debug() << "		" << sample << endl;

            if( sample == "ttJetsM0" ){
                nbins = 22;
                //                	m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min_ = 250.;
                mass_max_ = 2000.;
            }
            if( sample == "ttJetsM700" ){
                nbins = 10;
                //                	m_bins[nbins+1] = {700, 800, 900, 1000};
                mass_min_ = 700.;
                mass_max_ = 1000.;
            }
            if( sample == "ttJetsM1000" ){
                nbins = 10;
                //                	m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min_ = 1000.;
                mass_max_ = 2000.;
            }

            // Delta (gen-reco) Plots
            string delt = "Delta_Plots";

            book<TH1F>(delt, "Merged_THad_Delta_mass", "", mass_bins_, -2*mass_max_, 2*mass_max_);
            book<TH1F>(delt, "Merged_THad_Delta_costh", "", costh_bins_, -2*costh_max_, 2*costh_max_);

            book<TH1F>(delt, "Unmerged_THad_Delta_mass", "", mass_bins_, -2*mass_max_, 2*mass_max_);
            book<TH1F>(delt, "Unmerged_THad_Delta_costh", "", costh_bins_, -2*costh_max_, 2*costh_max_);


            //Reco Jets
            for( int i = 0; i < 2; i++ ){
                //Pt
                book<TH1F>(DRnames_[i], "Matched_perm_BHad_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_BLep_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJa_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJb_pt", "", pt_bins_, pt_min_, pt_max_);

                //Eta
                book<TH1F>(DRnames_[i], "Matched_perm_BHad_eta", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_BLep_eta", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJa_eta", "", eta_bins_, eta_min_, eta_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJb_eta", "", eta_bins_, eta_min_, eta_max_);

                //break ttJetsM0 into 700-1000 and > 1000
                book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM700_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM700_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM1000_pt", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM1000_pt", "", pt_bins_, pt_min_, pt_max_);

                book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM700_frac_p", "", 10, 0., 5.);
                book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM700_frac_p", "", 10, 0., 5.);
                book<TH1F>(DRnames_[i], "Matched_perm_WJa_ttbarM1000_frac_p", "", 10, 0., 5.);
                book<TH1F>(DRnames_[i], "Matched_perm_WJb_ttbarM1000_frac_p", "", 10, 0., 5.);

                //Matched objects
                book<TH1F>(DRnames_[i], "Matched_BHadWJa_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_BHadWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_WJaWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_BHadWJa_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_BHadWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_WJaWJb_ptthad_3J", "", pt_bins_, pt_min_, pt_max_);

                book<TH1F>(DRnames_[i], "Matched_BHadWJa_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_BHadWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_WJaWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_BHadWJa_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_BHadWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_WJaWJb_ptthad_4J", "", pt_bins_, pt_min_, pt_max_);

                book<TH1F>(DRnames_[i], "Matched_BHadWJa_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_BHadWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "Matched_WJaWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_BHadWJa_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_BHadWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);
                book<TH1F>(DRnames_[i], "All_Matched_WJaWJb_ptthad_5PJ", "", pt_bins_, pt_min_, pt_max_);

                //BTagging
                book<TH1F>(DRnames_[i], "BTag_less_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BHad_tight_fail", "", nbins, mass_min_, mass_max_);

                book<TH1F>(DRnames_[i], "BTag_less_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_less_BLep_tight_fail", "", nbins, mass_min_, mass_max_);

                book<TH1F>(DRnames_[i], "BTag_great_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BHad_tight_fail", "", nbins, mass_min_, mass_max_);

                book<TH1F>(DRnames_[i], "BTag_great_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "BTag_great_BLep_tight_fail", "", nbins, mass_min_, mass_max_);

                //Merged jets

                //BHad and WJa
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_perm_and_WJb_mass", "", 100, 0., 300.);// invariant mass of WJb and merged BHadWJa
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_perm_mass", "", 100, 0., 300.);// mass of merged BHad & WJa
                book<TH1F>(DRnames_[i], "Merged_WJb_mass", "", 100, 0., 200.);// WJb mass when BHad&WJa merged
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);//DR between gen l,BLep when BHadWJa merged
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_massDivpt", "", 50, 0., 0.5);// mass/pt for merged BHadWJa
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_DRBHadWJb", "", DR_bins_, DR_min_, DR_max_); //DR between merged BHadWJa and WJb

                //BHad and WJb
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_perm_and_WJa_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_perm_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_WJa_mass", "", 100, 0., 200.);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_massDivpt", "", 50, 0., 0.5);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_DRBHadWJa", "", DR_bins_, DR_min_, DR_max_); //DR between merged BHadWJb and WJa

                book<TH1F>(DRnames_[i], "Merged_BHadWJet_l100_pt", "", pt_bins_, pt_min_, pt_max_); //pt for merged bhad and W masses < 100
                book<TH1F>(DRnames_[i], "Merged_BHadWJet_l100_eta", "", eta_bins_, eta_min_, eta_max_); //eta for merged bhad and W masses < 100

                book<TH2D>(DRnames_[i], "Merged_BHadWJet_mass_vs_BHad_mass","",100, 0.,300.,100,0., 300.);

                //BLep and WJa
                book<TH1F>(DRnames_[i], "Merged_BLepWJa_perm_and_WJb_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_BLepWJa_perm_mass", "", 100, 0., 100.);
                book<TH1F>(DRnames_[i], "Merged_BLepWJa_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Merged_BLepWJa_massDivpt", "", 50, 0., 0.5);

                //BLep and WJb
                book<TH1F>(DRnames_[i], "Merged_BLepWJb_perm_and_WJa_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_BLepWJb_perm_mass", "", 100, 0., 100.);
                book<TH1F>(DRnames_[i], "Merged_BLepWJb_perm_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Merged_BLepWJb_massDivpt", "", 50, 0., 0.5);

                book<TH1F>(DRnames_[i], "Merged_BLepWJet_l100_pt", "", pt_bins_, pt_min_, pt_max_); //pt for merged blep and W masses < 100
                book<TH1F>(DRnames_[i], "Merged_BLepWJet_l100_eta", "", eta_bins_, eta_min_, eta_max_); //eta for merged blep and W masses < 100

                book<TH2D>(DRnames_[i], "Merged_BLepWJet_mass_vs_BLep_mass","",100, 0.,300.,100,0., 300.);

                //BHad and BLep
                book<TH1F>(DRnames_[i], "Merged_BHadBLep_perm_and_WJa_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_BHadBLep_perm_and_WJb_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_BHadBLep_perm_mass", "", 100, 0., 300.);
                book<TH1F>(DRnames_[i], "Merged_BHadBLep_massDivpt", "", 50, 0., 0.5);
                book<TH1F>(DRnames_[i], "Merged_BHadBLep_perm_DRBHadBLep", "", DR_bins_, DR_min_, DR_max_);

                book<TH2D>(DRnames_[i], "Merged_BsWJa_mass_vs_Bs_mass", "", 100, 0., 300., 100, 0., 300.);
                book<TH2D>(DRnames_[i], "Merged_BsWJb_mass_vs_Bs_mass", "", 100, 0., 300., 100, 0., 300.);

                //Unmerged jets
                book<TH1F>(DRnames_[i], "Unmerged_BHad_mass", "", 100, 0., 100.);
                book<TH1F>(DRnames_[i], "Unmerged_WJa_mass", "", 100, 0., 100.);
                book<TH1F>(DRnames_[i], "Unmerged_WJb_mass", "", 100, 0., 100.);
                book<TH1F>(DRnames_[i], "Unmerged_BHad_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Unmerged_WJa_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Unmerged_WJb_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Unmerged_BHad_massDivpt", "", 50, 0., 0.5);
                book<TH1F>(DRnames_[i], "Unmerged_WJa_massDivpt", "", 50, 0., 0.5);
                book<TH1F>(DRnames_[i], "Unmerged_WJb_massDivpt", "", 50, 0., 0.5);

                book<TH2D>(DRnames_[i], "Unmerged_BHadWJet_mass_vs_BHad_mass","",100, 0., 300., 100, 0., 300.);
                book<TH2D>(DRnames_[i], "Unmerged_BHadWJet_highest_mass_vs_BHad_mass", "", 100, 0., 300., 100, 0., 300.);

                book<TH1F>(DRnames_[i], "Unmerged_BLep_mass", "", 100, 0., 100.);
                book<TH1F>(DRnames_[i], "Unmerged_BLep_DRLepBLep", "", DR_bins_, DR_min_, DR_max_);
                book<TH1F>(DRnames_[i], "Unmerged_BLep_massDivpt", "", 50, 0., 0.5);

                book<TH2D>(DRnames_[i], "Unmerged_BLepWJet_mass_vs_BLep_mass","",100, 0., 300., 100, 0., 300.);
                book<TH2D>(DRnames_[i], "Unmerged_BLepWJet_highest_mass_vs_BLep_mass", "", 100, 0., 300., 100, 0., 300.);

                // Fraction of merged events
                //3 jets in events 
                //	                book<TH1F>(DRnames_[i], "Merged_BHadBLep_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "Merged_BLepWJa_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "Merged_BLepWJb_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_WJaWJb_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);

                //4 jets in events
                //	                book<TH1F>(DRnames_[i], "Merged_BHadBLep_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "Merged_BLepWJa_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "Merged_BLepWJb_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_WJaWJb_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);

                // 5 or more jets in events
                //	                book<TH1F>(DRnames_[i], "Merged_BHadBLep_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJa_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_BHadWJb_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "Merged_BLepWJa_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "Merged_BLepWJb_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "Merged_WJaWJb_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);

                book<TH1F>(DRnames_[i], "All_BHad_events_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "All_BLep_events_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "All_WJa_events_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "All_WJb_events_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);

                book<TH1F>(DRnames_[i], "All_BHad_events_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "All_BLep_events_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "All_WJa_events_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "All_WJb_events_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);

                book<TH1F>(DRnames_[i], "All_BHad_events_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                //	                book<TH1F>(DRnames_[i], "All_BLep_events_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "All_WJa_events_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);
                book<TH1F>(DRnames_[i], "All_WJb_events_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);

                // Fraction of 3 jets from hadronic side merging
                book<TH1F>(DRnames_[i], "Merged_BHadWJaWJb_vs_mttbar_3J", "", nbins, mass_min_, mass_max_);// 3 jets in events
                book<TH1F>(DRnames_[i], "Merged_BHadWJaWJb_vs_mttbar_4J", "", nbins, mass_min_, mass_max_);// 4 jets in events
                book<TH1F>(DRnames_[i], "Merged_BHadWJaWJb_vs_mttbar_5PJ", "", nbins, mass_min_, mass_max_);// 5 or more jets in events

                book<TH1F>(DRnames_[i], "Merged_BHadWJaWJb_mass", "", 100, 0., 300.);			
                book<TH1F>(DRnames_[i], "Merged_BHadWJaWJb_pt", "", pt_bins_, pt_min_, pt_max_);			
                book<TH1F>(DRnames_[i], "Merged_BHadWJaWJb_eta", "", eta_bins_, eta_min_, eta_max_);			

            }

            Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }

        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

            URStreamer event(tree_);

            while(event.next() /*&& evt_idx_ < 30*/)
            {
                evt_idx_++;
                if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

                auto delt_dir = histos_.find("Delta_Plots");

                //generator selection
                bool selection = genp_selector_.select(event);
                tracker_.track("gen selection");
                if( !selection ){
                    Logger::log().debug() << "event has no selection " << endl;
                    continue;
                }
                GenTTBar &ttbar = genp_selector_.ttbar_system();

                Permutation matched_perm;
                if( object_selector_.select(event) ){
                    if( ttbar.type == GenTTBar::DecayType::SEMILEP ){
                        matched_perm = matcher_.match(
                                genp_selector_.ttbar_final_system(),
                                object_selector_.clean_jets(),
                                object_selector_.loose_electrons(),
                                object_selector_.loose_muons());
                    }
                    matched_perm.SetMET(object_selector_.met());
                }

                if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) continue;
                tracker_.track("semilep");

                njets_ = object_selector_.clean_jets().size();

                //cout << "njets = " << njets_ << endl;

                if( njets_ < 3 ) continue;
                tracker_.track("njet cut");


                //initialize gen jets
                //      			const GenObject* lepton = 0;
                const GenObject* BLep = 0;
                const GenObject* BHad = 0;
                const GenObject* WJa = 0;
                const GenObject* WJb = 0;
                //			const GenTop* tlep = 0;
                const GenTop* thad = 0;

                //Generator level definitions
                //			if( ttbar.lep_top() ) tlep = ttbar.lep_top();
                if( ttbar.had_top() ) thad = ttbar.had_top();
                //			if( ttbar.lepton() ) lepton = ttbar.lepton();
                if( ttbar.lep_b() ) BLep = ttbar.lep_b();
                if( ttbar.had_b() ) BHad = ttbar.had_b();
                if( ttbar.had_W() ){
                    if( ttbar.had_W()->first && !ttbar.had_W()->second ) WJa = ttbar.had_W()->first;
                    if( ttbar.had_W()->second && !ttbar.had_W()->first ) WJb = ttbar.had_W()->second;
                    if( ttbar.had_W()->first && ttbar.had_W()->second ){
                        WJa = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
                        WJb = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;
                    }
                }

                ///////////////Perm Matching
                for( int i = 0; i < 2; i++ ){

                    //				if( njets_ < 3 ) continue;//events must have at least 3 jets

                    auto match_dir = histos_.find(DRnames_[i]);

                    //initialize permutation objects
                    //jets
                    const IDJet* best_BHad = 0;
                    const IDJet* best_BLep = 0;
                    const IDJet* best_WJa = 0;
                    const IDJet* best_WJb = 0;


                    list<IDJet*> bhad_list;
                    float bhad_DR = 1e10;// itialize dr to high number
                    for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
                        if( !(BHad == 0) ){
                            if( (*jets)->DeltaR(*BHad) > DR_[i] ) continue;
                            if( (*jets)->DeltaR(*BHad) < bhad_DR ){
                                bhad_DR = (*jets)->DeltaR(*BHad);
                                bhad_list.push_back(*jets);
                                best_BHad = (bhad_list.back());
                            }
                        }
                    }
                    list<IDJet*> blep_list;
                    float blep_DR = 1e10;// itialize dr to high number
                    for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
                        if( !(BLep == 0) ){
                            if( (*jets)->DeltaR(*BLep) > DR_[i] ) continue;
                            if( (*jets)->DeltaR(*BLep) < blep_DR ){
                                blep_DR = (*jets)->DeltaR(*BLep);
                                blep_list.push_back(*jets);
                                best_BLep = (blep_list.back());
                            }
                        }
                    }
                    list<IDJet*> wja_list;
                    float wja_DR = 1e10;// itialize dr to high number
                    for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
                        if( !(WJa == 0) ){
                            if( (*jets)->DeltaR(*WJa) > DR_[i] ) continue;
                            if( (*jets)->DeltaR(*WJa) < wja_DR ){
                                wja_DR = (*jets)->DeltaR(*WJa);
                                wja_list.push_back(*jets);
                                best_WJa = (wja_list.back());
                            }
                        }
                    }
                    list<IDJet*> wjb_list;
                    float wjb_DR = 1e10;// itialize dr to high number
                    for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
                        if( !(WJb == 0) ){
                            if( (*jets)->DeltaR(*WJb) > DR_[i] ) continue;
                            if( (*jets)->DeltaR(*WJb) < wjb_DR ){
                                wjb_DR = (*jets)->DeltaR(*WJb);
                                wjb_list.push_back(*jets);
                                best_WJb = (wjb_list.back());
                            }
                        }
                    }

                    //Reco Obj. Kinematic Vars
                    if( !(best_BHad == 0) ){//reco BHad exists
                        match_dir->second["Matched_perm_BHad_pt"].fill(best_BHad->Pt());
                        match_dir->second["Matched_perm_BHad_eta"].fill(best_BHad->Eta());

                        if( njets_ == 3) match_dir->second["All_BHad_events_vs_mttbar_3J"].fill(ttbar.M());
                        if( njets_ == 4) match_dir->second["All_BHad_events_vs_mttbar_4J"].fill(ttbar.M());
                        if( njets_ >= 5) match_dir->second["All_BHad_events_vs_mttbar_5PJ"].fill(ttbar.M());

                        if( best_BHad == best_BLep || best_BHad == best_WJa || best_BHad == best_WJb ){//reco BHad merged with something
                            if( best_BHad->BTagId(cut_loose_b_) ) match_dir->second["BTag_less_BHad_loose_pass"].fill(ttbar.M());
                            if( best_BHad->BTagId(cut_medium_b_) ) match_dir->second["BTag_less_BHad_medium_pass"].fill(ttbar.M());	
                            if( best_BHad->BTagId(cut_tight_b_) ) match_dir->second["BTag_less_BHad_tight_pass"].fill(ttbar.M());
                            if( !best_BHad->BTagId(cut_loose_b_) ) match_dir->second["BTag_less_BHad_loose_fail"].fill(ttbar.M());
                            if( !best_BHad->BTagId(cut_medium_b_) ) match_dir->second["BTag_less_BHad_medium_fail"].fill(ttbar.M());	
                            if( !best_BHad->BTagId(cut_tight_b_) ) match_dir->second["BTag_less_BHad_tight_fail"].fill(ttbar.M());
                        }
                    }
                    if( !(best_BLep == 0) ){//reco BLep exists
                        match_dir->second["Matched_perm_BLep_pt"].fill(best_BLep->Pt());
                        match_dir->second["Matched_perm_BLep_eta"].fill(best_BLep->Eta());
                        //	                		match_dir->second["All_BLep_events_vs_mttbar"].fill(ttbar.M());

                        if( best_BLep == best_BHad || best_BLep == best_WJa || best_BLep == best_WJb ){//reco BLep merged with something
                            if( best_BLep->BTagId(cut_loose_b_) ) match_dir->second["BTag_less_BLep_loose_pass"].fill(ttbar.M());
                            if( best_BLep->BTagId(cut_medium_b_) ) match_dir->second["BTag_less_BLep_medium_pass"].fill(ttbar.M());	
                            if( best_BLep->BTagId(cut_tight_b_) ) match_dir->second["BTag_less_BLep_tight_pass"].fill(ttbar.M());
                            if( !best_BLep->BTagId(cut_loose_b_) ) match_dir->second["BTag_less_BLep_loose_fail"].fill(ttbar.M());
                            if( !best_BLep->BTagId(cut_medium_b_) ) match_dir->second["BTag_less_BLep_medium_fail"].fill(ttbar.M());	
                            if( !best_BLep->BTagId(cut_tight_b_) ) match_dir->second["BTag_less_BLep_tight_fail"].fill(ttbar.M());
                        }
                    }
                    if( !(best_WJa == 0) ){//reco WJa exists
                        match_dir->second["Matched_perm_WJa_pt"].fill(best_WJa->Pt());
                        match_dir->second["Matched_perm_WJa_eta"].fill(best_WJa->Eta());

                        if( njets_ == 3 ) match_dir->second["All_WJa_events_vs_mttbar_3J"].fill(ttbar.M());
                        if( njets_ == 4 ) match_dir->second["All_WJa_events_vs_mttbar_4J"].fill(ttbar.M());
                        if( njets_ >= 5 ) match_dir->second["All_WJa_events_vs_mttbar_5PJ"].fill(ttbar.M());

                        if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                            match_dir->second["Matched_perm_WJa_ttbarM700_pt"].fill(best_WJa->Pt());
                            match_dir->second["Matched_perm_WJa_ttbarM700_frac_p"].fill(best_WJa->P()/(ttbar.M()/2));
                        }
                        if( ttbar.M() >= 1000 ){
                            match_dir->second["Matched_perm_WJa_ttbarM1000_pt"].fill(best_WJa->Pt());
                            match_dir->second["Matched_perm_WJa_ttbarM1000_frac_p"].fill(best_WJa->P()/(ttbar.M()/2));
                        }
                    }
                    if( !(best_WJb == 0) ){//reco WJb exists
                        match_dir->second["Matched_perm_WJb_pt"].fill(best_WJb->Pt());
                        match_dir->second["Matched_perm_WJb_eta"].fill(best_WJb->Eta());

                        if( njets_ == 3 ) match_dir->second["All_WJb_events_vs_mttbar_3J"].fill(ttbar.M());
                        if( njets_ == 4 ) match_dir->second["All_WJb_events_vs_mttbar_4J"].fill(ttbar.M());
                        if( njets_ >= 5 ) match_dir->second["All_WJb_events_vs_mttbar_5PJ"].fill(ttbar.M());

                        if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                            match_dir->second["Matched_perm_WJb_ttbarM700_pt"].fill(best_WJb->Pt());
                            match_dir->second["Matched_perm_WJb_ttbarM700_frac_p"].fill(best_WJb->P()/(ttbar.M()/2));
                        }
                        if( ttbar.M() >= 1000 ){
                            match_dir->second["Matched_perm_WJb_ttbarM1000_pt"].fill(best_WJb->Pt());
                            match_dir->second["Matched_perm_WJb_ttbarM1000_frac_p"].fill(best_WJb->P()/(ttbar.M()/2));
                        }
                    }
                    if( !(best_BLep == 0) && !(best_BLep == best_BHad) && !(best_BLep == best_WJa) && !(best_BLep == best_WJb) ){//reco BLep not merged with anything
                        if( best_BLep->BTagId(cut_loose_b_) ) match_dir->second["BTag_great_BLep_loose_pass"].fill(ttbar.M());
                        if( best_BLep->BTagId(cut_medium_b_) ) match_dir->second["BTag_great_BLep_medium_pass"].fill(ttbar.M());	
                        if( best_BLep->BTagId(cut_tight_b_) ) match_dir->second["BTag_great_BLep_tight_pass"].fill(ttbar.M());
                        if( !best_BLep->BTagId(cut_loose_b_) ) match_dir->second["BTag_great_BLep_loose_fail"].fill(ttbar.M());
                        if( !best_BLep->BTagId(cut_medium_b_) ) match_dir->second["BTag_great_BLep_medium_fail"].fill(ttbar.M());	
                        if( !best_BLep->BTagId(cut_tight_b_) ) match_dir->second["BTag_great_BLep_tight_fail"].fill(ttbar.M());

                        match_dir->second["Unmerged_BLep_mass"].fill(best_BLep->M());
                        match_dir->second["Unmerged_BLep_massDivpt"].fill(best_BLep->M()/best_BLep->Pt());
                        if( !(best_WJa == 0) ){
                            match_dir->second["Unmerged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJa).M());
                        }
                        if( !(best_WJb == 0) ){
                            match_dir->second["Unmerged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJb).M());
                        }
                        if( !(best_WJa == 0) && !(best_WJb == 0) ){
                            if( best_WJa->M() > best_WJb->M() ) match_dir->second["Unmerged_BLepWJet_highest_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJa).M());
                            if( best_WJa->M() < best_WJb->M() ) match_dir->second["Unmerged_BLepWJet_highest_mass_vs_BLep_mass"].fill(best_BLep->M(),(*best_BLep+*best_WJb).M());
                        }


                    }
                    if( !(best_BHad == 0) && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) && !(best_BHad == best_BLep) ){//reco BHad not merged with anything
                        if( best_BHad->BTagId(cut_loose_b_) ) match_dir->second["BTag_great_BHad_loose_pass"].fill(ttbar.M());
                        if( best_BHad->BTagId(cut_medium_b_) ) match_dir->second["BTag_great_BHad_medium_pass"].fill(ttbar.M());	
                        if( best_BHad->BTagId(cut_tight_b_) ) match_dir->second["BTag_great_BHad_tight_pass"].fill(ttbar.M());
                        if( !best_BHad->BTagId(cut_loose_b_) ) match_dir->second["BTag_great_BHad_loose_fail"].fill(ttbar.M());
                        if( !best_BHad->BTagId(cut_medium_b_) ) match_dir->second["BTag_great_BHad_medium_fail"].fill(ttbar.M());	
                        if( !best_BHad->BTagId(cut_tight_b_) ) match_dir->second["BTag_great_BHad_tight_fail"].fill(ttbar.M());

                        if( thad && best_WJa && best_WJb && !(best_WJa == best_WJb) && !(best_WJa == best_BLep) && !(best_WJb == best_BLep) ){//reco WJa adn WJb exist but aren't merged to anything
                            delt_dir->second["Unmerged_THad_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa+*best_WJb).M());
                            delt_dir->second["Unmerged_THad_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa+*best_WJb).CosTheta());
                        }

                        match_dir->second["Unmerged_BHad_mass"].fill(best_BHad->M());
                        match_dir->second["Unmerged_BHad_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
                        if( best_WJa ){
                            match_dir->second["Unmerged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJa).M());
                        }

                        if( best_WJb ){
                            match_dir->second["Unmerged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJb).M());
                        }

                        if( !(best_WJa == 0) && !(best_WJb == 0) ){
                            if( best_WJa->M() > best_WJb->M() ) match_dir->second["Unmerged_BHadWJet_highest_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJa).M());
                            if( best_WJa->M() < best_WJb->M() ) match_dir->second["Unmerged_BHadWJet_highest_mass_vs_BHad_mass"].fill(best_BHad->M(),(*best_BHad+*best_WJb).M());
                        }
                    }
                    if( !(best_WJa == 0) && !(best_WJa == best_BHad) && !(best_WJa == best_WJb) && !(best_WJa == best_BLep) ){//reco WJa not merged with anything

                        match_dir->second["Unmerged_WJa_massDivpt"].fill(best_WJa->M()/best_WJa->Pt());
                        match_dir->second["Unmerged_WJa_mass"].fill(best_WJa->M());
                    }
                    if( !(best_WJb == 0) && !(best_WJb == best_BHad) && !(best_WJb == best_WJa) && !(best_WJb == best_BLep) ){//reco WJb not merged with anything

                        match_dir->second["Unmerged_WJb_massDivpt"].fill(best_WJb->M()/best_WJb->Pt());
                        match_dir->second["Unmerged_WJb_mass"].fill(best_WJb->M());
                    }
                    if( !(best_BHad == 0) && !(best_WJa == 0) ){ //reco BHad and WJa exist
                        if( njets_ == 3 ) match_dir->second["All_Matched_BHadWJa_ptthad_3J"].fill(thad->Pt());
                        if( njets_ == 4 ) match_dir->second["All_Matched_BHadWJa_ptthad_4J"].fill(thad->Pt());
                        if( njets_ > 4 ) match_dir->second["All_Matched_BHadWJa_ptthad_5PJ"].fill(thad->Pt());
                    }
                    if( !(best_BHad == 0) && !(best_WJb == 0) ){ //reco BHad and WJb exist
                        if( njets_ == 3 ) match_dir->second["All_Matched_BHadWJb_ptthad_3J"].fill(thad->Pt());
                        if( njets_ == 4 ) match_dir->second["All_Matched_BHadWJb_ptthad_4J"].fill(thad->Pt());
                        if( njets_ > 4 ) match_dir->second["All_Matched_BHadWJb_ptthad_5PJ"].fill(thad->Pt());
                    }
                    if( !(best_WJa == 0) && !(best_WJb == 0) ){ //reco WJa and WJb exist
                        if( njets_ == 3 ) match_dir->second["All_Matched_WJaWJb_ptthad_3J"].fill(thad->Pt());
                        if( njets_ == 4 ) match_dir->second["All_Matched_WJaWJb_ptthad_4J"].fill(thad->Pt());
                        if( njets_ > 4 ) match_dir->second["All_Matched_WJaWJb_ptthad_5PJ"].fill(thad->Pt());
                    }
                    if( best_BHad == best_WJa && !(best_BHad == 0) && !(best_BHad == best_BLep) && !(best_BHad == best_WJb) ){//only reco BHad and WJa merged

                        if( njets_ == 3 ){
                            match_dir->second["Merged_BHadWJa_vs_mttbar_3J"].fill(ttbar.M());
                            match_dir->second["Matched_BHadWJa_ptthad_3J"].fill(thad->Pt());
                        }
                        if( njets_ == 4 ){
                            match_dir->second["Merged_BHadWJa_vs_mttbar_4J"].fill(ttbar.M());
                            match_dir->second["Matched_BHadWJa_ptthad_4J"].fill(thad->Pt());
                        }
                        if( njets_ >= 5 ){
                            match_dir->second["Merged_BHadWJa_vs_mttbar_5PJ"].fill(ttbar.M());
                            match_dir->second["Matched_BHadWJa_ptthad_5PJ"].fill(thad->Pt());
                        }

                        match_dir->second["Merged_BHadWJa_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
                        match_dir->second["Merged_BHadWJa_perm_mass"].fill(best_BHad->M());//merged b and wjet mass
                        if( !(best_WJb == 0) && !(best_WJb == best_BHad) && !(best_WJb == best_BLep) ){//reco WJb exists and not merged
                            match_dir->second["Merged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(), (*best_BHad+*best_WJb).M());
                            match_dir->second["Merged_BHadWJa_perm_and_WJb_mass"].fill((*best_BHad+*best_WJb).M());//comb invariant mass of merged jet and other wjet

                            match_dir->second["Merged_BHadWJa_DRBHadWJb"].fill(best_BHad->DeltaR(*best_WJb));

                            if( (*best_BHad+*best_WJb).M() < 100 ){
                                match_dir->second["Merged_BHadWJet_l100_pt"].fill(best_BHad->Pt());
                                match_dir->second["Merged_BHadWJet_l100_eta"].fill(best_BHad->Eta());
                            }

                            match_dir->second["Merged_WJb_mass"].fill(best_WJb->M());

                            if( thad ){
                                delt_dir->second["Merged_THad_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJb).M());
                                delt_dir->second["Merged_THad_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJb).CosTheta());
                            }
                        }
                    }
                    if( best_BHad == best_WJb && !(best_BHad == 0) && !(best_BHad == best_BLep) && !(best_BHad == best_WJa) ){//only reco BHad and WJb merged

                        if( njets_ == 3 ){
                            match_dir->second["Merged_BHadWJb_vs_mttbar_3J"].fill(ttbar.M());
                            match_dir->second["Matched_BHadWJb_ptthad_3J"].fill(thad->Pt());
                        }
                        if( njets_ == 4 ){
                            match_dir->second["Merged_BHadWJb_vs_mttbar_4J"].fill(ttbar.M());
                            match_dir->second["Matched_BHadWJb_ptthad_4J"].fill(thad->Pt());
                        }
                        if( njets_ >= 5 ){
                            match_dir->second["Merged_BHadWJb_vs_mttbar_5PJ"].fill(ttbar.M());
                            match_dir->second["Matched_BHadWJb_ptthad_5PJ"].fill(thad->Pt());
                        }

                        match_dir->second["Merged_BHadWJb_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
                        match_dir->second["Merged_BHadWJb_perm_mass"].fill(best_BHad->M());//merged b and wjet mass
                        if( !(best_WJa == 0) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) ){//reco WJa exists and not merged
                            match_dir->second["Merged_BHadWJet_mass_vs_BHad_mass"].fill(best_BHad->M(), (*best_BHad+*best_WJa).M());
                            match_dir->second["Merged_BHadWJb_perm_and_WJa_mass"].fill((*best_BHad+*best_WJa).M());//comb mass of merged jet and other wjet

                            match_dir->second["Merged_BHadWJb_DRBHadWJa"].fill(best_BHad->DeltaR(*best_WJa));

                            if( (*best_BHad+*best_WJa).M() < 100 ){
                                match_dir->second["Merged_BHadWJet_l100_pt"].fill(best_BHad->Pt());
                                match_dir->second["Merged_BHadWJet_l100_eta"].fill(best_BHad->Eta());
                            }

                            match_dir->second["Merged_WJa_mass"].fill(best_WJa->M());

                            if( thad ){
                                delt_dir->second["Merged_THad_Delta_mass"].fill(thad->M()-(*best_BHad+*best_WJa).M());
                                delt_dir->second["Merged_THad_Delta_costh"].fill(thad->CosTheta()-(*best_BHad+*best_WJa).CosTheta());
                            }
                        }
                    }
                    if( best_WJa == best_WJb && !(best_WJa == 0) && !(best_WJa == best_BHad) && !(best_WJa == best_BLep) ){//only reco WJa and WJb merged

                        if( njets_ == 3 ){
                            match_dir->second["Matched_WJaWJb_ptthad_3J"].fill(thad->Pt());
                            match_dir->second["Merged_WJaWJb_vs_mttbar_3J"].fill(ttbar.M());
                        }
                        if( njets_ == 4 ){
                            match_dir->second["Matched_WJaWJb_ptthad_4J"].fill(thad->Pt());
                            match_dir->second["Merged_WJaWJb_vs_mttbar_4J"].fill(ttbar.M());
                        }
                        if( njets_ > 4 ){
                            match_dir->second["Matched_WJaWJb_ptthad_5PJ"].fill(thad->Pt());
                            match_dir->second["Merged_WJaWJb_vs_mttbar_5PJ"].fill(ttbar.M());
                        }

                        if( thad && best_BHad ){
                            delt_dir->second["Merged_THad_Delta_mass"].fill(thad->M()-(*best_WJa+*best_BHad).M());
                            delt_dir->second["Merged_THad_Delta_costh"].fill(thad->CosTheta()-(*best_WJa+*best_BHad).CosTheta());
                        }
                    }

                    if( best_BLep == best_WJa && !(best_BLep == 0) && !(best_BHad == best_BLep) && !(best_BLep == best_WJb) ){//only reco BLep and WJa merged

                        //	                		match_dir->second["Merged_BLepWJa_vs_mttbar"].fill(ttbar.M());

                        match_dir->second["Merged_BLepWJa_massDivpt"].fill(best_BLep->M()/best_BLep->Pt());
                        match_dir->second["Merged_BLepWJa_perm_mass"].fill(best_BLep->M());//merged b and wjet mass
                        if( !(best_WJb == 0) && !(best_WJb == best_BLep) && !(best_WJb == best_BHad) ){//reco WJb exists and not merged
                            match_dir->second["Merged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(), (*best_BLep+*best_WJb).M());
                            match_dir->second["Merged_BLepWJa_perm_and_WJb_mass"].fill((*best_BLep+*best_WJb).M());//comb invariant mass of merged jet and other wjet
                            if( (*best_BLep+*best_WJb).M() < 100 ){
                                match_dir->second["Merged_BLepWJet_l100_pt"].fill(best_BLep->Pt());
                                match_dir->second["Merged_BLepWJet_l100_eta"].fill(best_BLep->Eta());
                            }

                            //match_dir->second["Merged_WJb_mass"].fill(best_WJb_DRP4->M());
                        }
                        //if( object_selector_.clean_jets().size() == 3 ) match_dir->second["Matched_BHadWJa_ptthad_3J"].fill(thad->Pt());
                        //if( object_selector_.clean_jets().size() == 4 ) match_dir->second["Matched_BHadWJa_ptthad_4J"].fill(thad->Pt());
                        //if( object_selector_.clean_jets().size() > 4 ) match_dir->second["Matched_BHadWJa_ptthad_5PJ"].fill(thad->Pt());
                    }
                    if( best_BLep == best_WJb && !(best_BLep == 0) && !(best_BHad == best_BLep) && !(best_BLep == best_WJa) ){//only reco BLep and WJb merged

                        //	                		match_dir->second["Merged_BLepWJb_vs_mttbar"].fill(ttbar.M());

                        match_dir->second["Merged_BLepWJb_massDivpt"].fill(best_BLep->M()/best_BLep->Pt());
                        match_dir->second["Merged_BLepWJb_perm_mass"].fill(best_BLep->M());//merged b and wjet mass
                        if( !(best_WJa == 0) && !(best_WJa == best_BLep) && !(best_WJa == best_BHad) ){//reco WJa exists and not merged
                            match_dir->second["Merged_BLepWJet_mass_vs_BLep_mass"].fill(best_BLep->M(), (*best_BLep+*best_WJa).M());
                            match_dir->second["Merged_BLepWJb_perm_and_WJa_mass"].fill((*best_BLep+*best_WJa).M());//comb mass of merged jet and other wjet
                            if( (*best_BLep+*best_WJa).M() < 100 ){
                                match_dir->second["Merged_BLepWJet_l100_pt"].fill(best_BLep->Pt());
                                match_dir->second["Merged_BLepWJet_l100_eta"].fill(best_BLep->Eta());
                            }

                            //	match_dir->second["Merged_WJa_mass"].fill(best_WJa_DRP4->M());
                        }
                        //if( object_selector_.clean_jets().size() == 3 ) match_dir->second["Matched_BHadWJb_ptthad_3J"].fill(thad->Pt());
                        //if( object_selector_.clean_jets().size() == 4 ) match_dir->second["Matched_BHadWJb_ptthad_4J"].fill(thad->Pt());
                        //if( object_selector_.clean_jets().size() > 4 ) match_dir->second["Matched_BHadWJb_ptthad_5PJ"].fill(thad->Pt());
                    }
                    if( (best_BHad == best_BLep) && !(best_BHad == 0) && !(best_BHad == best_WJa) && !(best_BHad == best_WJb) ){//only reco BHad and BLep merged

                        //	                		match_dir->second["Merged_BHadBLep_vs_mttbar"].fill(ttbar.M());

                        match_dir->second["Merged_BHadBLep_massDivpt"].fill(best_BHad->M()/best_BHad->Pt());
                        match_dir->second["Merged_BHadBLep_perm_mass"].fill(best_BHad->M());
                        if( !(best_WJa == 0) && !(best_WJa == best_WJb) && !(best_WJa == best_BHad) ){//reco WJa exists and not merged
                            match_dir->second["Merged_BsWJa_mass_vs_Bs_mass"].fill(best_BHad->M(), (*best_BHad+*best_WJa).M());
                            match_dir->second["Merged_BHadBLep_perm_and_WJa_mass"].fill((*best_BHad+*best_WJa).M());
                        }
                        if( !(best_WJb == 0) && !(best_WJb == best_WJa) && !(best_WJb == best_BHad) ){//reco WJb exists and not merged
                            match_dir->second["Merged_BsWJb_mass_vs_Bs_mass"].fill(best_BHad->M(), (*best_BHad+*best_WJb).M());
                            match_dir->second["Merged_BHadBLep_perm_and_WJb_mass"].fill((*best_BHad+*best_WJb).M());
                        }
                    }

                    // BHad, WJa, and WJb are all merged, but not with BLep
                    if( !(best_BHad == 0) && (best_BHad == best_WJa) && (best_BHad == best_WJb) && !(best_BHad == best_BLep) ){
                        if( njets_ == 3 ) match_dir->second["Merged_BHadWJaWJb_vs_mttbar_3J"].fill(ttbar.M());
                        if( njets_ == 4 ) match_dir->second["Merged_BHadWJaWJb_vs_mttbar_4J"].fill(ttbar.M());
                        if( njets_ >= 5 ) match_dir->second["Merged_BHadWJaWJb_vs_mttbar_5PJ"].fill(ttbar.M());

                        match_dir->second["Merged_BHadWJaWJb_mass"].fill(best_BHad->M());
                        match_dir->second["Merged_BHadWJaWJb_pt"].fill(best_BHad->Pt());
                        match_dir->second["Merged_BHadWJaWJb_eta"].fill(best_BHad->Eta());
                    }

                }

            }

            Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
        }

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
    URDriver<jet_match> test;
    int thing = test.run();
    return thing;
}
