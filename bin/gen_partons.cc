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

class gen_partons : public AnalyzerBase
{
    private:
        //counters
        unsigned long evt_idx_ = 0; //event index
        double njets_ = 0;

        const char *DRnames_[2] = {"DRP4", "DRP8"};
        double DR_[2] = {0.4, 0.8};
        const char *GenObjNames_[5] = {"Lep", "BLep", "BHad", "WJa", "WJb"};
        const char *GenJetNames_[4] = {"BLep", "BHad", "WJa", "WJb"};
        const char *GenTopNames_[2] = {"THad", "TLep"};
        //            const char *DRnames_[4] = {"DRP4", "DRP5", "DRP6", "DRP8"};
        //	    double DR_[4] = {0.4, 0.5, 0.6, 0.8};
        //            const char *DRnames_[1] = {"DRP4"};
        //	    double DR_[1] = {0.4};

//        decltype map< string, GenObject* > GenObjs_ = {
//            {"lepton", GenObject* lepton},
//            {"BHad", GenObject* BHad},
//            {"BLep", GenObject* BLep},
//            {"WJa", GenObject* WJa},
//            {"WJb", GenObject* WJb}
//        };
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
        int deltphi_bins_ = 100;

        double DR_min_ = 0.;
        double DR_max_ = 8.;
        double pt_min_ = 0.;
        double pt_max_ = 2000.;
        double eta_min_ = -5.0;
        double eta_max_ = 5.0;
        double costh_min_ = -1.;
        double costh_max_ = 1.;
        double deltphi_min_ = -4*Pi();
        double deltphi_max_ = 4*Pi();


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


        float cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_;


    public:
        gen_partons(const std::string output_filename):
            AnalyzerBase("gen_partons", output_filename),
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

        // get parameters from cfg file
        URParser& parser = URParser::instance();

        parser.addCfgParameter<float>("gen_jets", "ptmin", "minimum pt");
        parser.addCfgParameter<float>("gen_jets", "etamax", "maximum eta");
        parser.addCfgParameter<float>("gen_jets", "lead_ptmin", "minimum leading jet pt");
        parser.parseArguments();
        
        cut_jet_ptmin_ = parser.getCfgPar<float>("gen_jets", "ptmin" );
        cut_jet_etamax_ = parser.getCfgPar<float>("gen_jets", "etamax");
        cut_leadjet_ptmin_ = parser.getCfgPar<float>("gen_jets", "lead_ptmin" );

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
            else if( sample == "ttJetsM700" ){
                nbins = 10;
                //                	m_bins[nbins+1] = {700, 800, 900, 1000};
                mass_min_ = 700.;
                mass_max_ = 1000.;
            }
            else if( sample == "ttJetsM1000" ){
                nbins = 10;
                //                	m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
                mass_min_ = 1000.;
                mass_max_ = 2000.;
            }
            else{
                nbins = 20;
                mass_min_ = 250.;
                mass_max_ = 1000.;
            }

            string folder = "Gen_Plots";

            //Gen Objects
            	//DR
            //book<TH2D>(folder, "DR_LepBHad_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_LepBLep_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_LepWJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_LepWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_BHadBLep_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_BHadWJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_BHadWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_BLepWJa_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_BLepWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DR_WJaWJb_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //
            //book<TH1F>(folder, "DR_LepBHad", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_LepBLep", "", DR_bins_, DR_min_, DR_max_) ;
            //book<TH1F>(folder, "DR_LepWJa", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_LepWJb", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_BHadBLep", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_BHadWJa", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_BHadWJb", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_BLepWJa", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_BLepWJb", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DR_WJaWJb", "", DR_bins_, DR_min_, DR_max_);
            //
            //book<TH2D>(folder, "DRmin_THad_vs_mttbar", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DRmin_TLep_vs_mttbar", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DRmin_THad_vs_ptTHad", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DRmin_TLep_vs_ptTLep", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);
            //
            //book<TH1F>(folder, "DRmin_THad", "", DR_bins_, DR_min_, DR_max_);
            //book<TH1F>(folder, "DRmin_TLep", "", DR_bins_, DR_min_, DR_max_);
            //
            //book<TH2D>(folder, "DRmax_THad_vs_mttbar", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
            //book<TH2D>(folder, "DRmax_THad_vs_ptTHad", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);
            //
            //book<TH1F>(folder, "DRmax_THad", "", DR_bins_, DR_min_, DR_max_);
            
            for( int i = 0; i < 2; i++ ){
            //			book<TH1F>(folder, std::string("DRmin_THad_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            //	        book<TH1F>(folder, std::string("DRmin_THad_l")+DRnames_[i]+"_vs_ptTHad", "", nbins, pt_min_, pt_max_);
            //	        book<TH1F>(folder, std::string("DRmin_TLep_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            //	        book<TH1F>(folder, std::string("DRmin_TLep_l")+DRnames_[i]+"_vs_ptTLep", "", nbins, pt_min_, pt_max_);
            //	
            //	        book<TH1F>(folder, std::string("DRmin_THad_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            //	        book<TH1F>(folder, std::string("DRmin_THad_g")+DRnames_[i]+"_vs_ptTHad", "", nbins, pt_min_, pt_max_);
            //	        book<TH1F>(folder, std::string("DRmin_TLep_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            //	        book<TH1F>(folder, std::string("DRmin_TLep_g")+DRnames_[i]+"_vs_ptTLep", "", nbins, pt_min_, pt_max_);
            //
            //	        book<TH1F>(folder, std::string("DRmax_THad_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            //	        book<TH1F>(folder, std::string("DRmax_THad_l")+DRnames_[i]+"_vs_ptTHad", "", nbins, pt_min_, pt_max_);
            //	
            //	        book<TH1F>(folder, std::string("DRmax_THad_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            //	        book<TH1F>(folder, std::string("DRmax_THad_g")+DRnames_[i]+"_vs_ptTHad", "", nbins, pt_min_, pt_max_);
            
                book<TH1F>(folder, std::string("DR_LepBHad_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_LepBLep_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_LepWJa_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_LepWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BHadBLep_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BHadWJa_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BHadWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BLepWJa_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BLepWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_WJaWJb_l")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            
                book<TH1F>(folder, std::string("DR_LepBHad_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_LepBLep_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_LepWJa_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_LepWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BHadBLep_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BHadWJa_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BHadWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BLepWJa_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_BLepWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
                book<TH1F>(folder, std::string("DR_WJaWJb_g")+DRnames_[i]+"_vs_mttbar", "", nbins, mass_min_, mass_max_);
            }


            string pt_plots = "Pt_Plots";
            string eta_plots = "Eta_Plots";
            string costh_plots = "Costh_Plots";
            string mass_plots = "Mass_Plots";

                // Plots for gen objects (lep, bhad, blep ...)
            for( const auto& genobj : GenObjNames_ ){ 
                    //Pt
                book<TH2D>(pt_plots, std::string("Pt_")+genobj+"_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                book<TH1F>(pt_plots, std::string("Pt_")+genobj, "", pt_bins_, pt_min_, pt_max_);

                    //Eta
                book<TH2D>(eta_plots, std::string("Eta_")+genobj+"_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                book<TH1F>(eta_plots, std::string("Eta_")+genobj, "", eta_bins_, eta_min_, eta_max_);

                    //Costh
                book<TH2D>(costh_plots, std::string("Costh_")+genobj+"_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, costh_bins_, costh_min_, costh_max_);
                book<TH1F>(costh_plots, std::string("Costh_")+genobj, "", costh_bins_, costh_min_, costh_max_);

                    //Mass
                book<TH1F>(mass_plots, std::string("Mass_")+genobj, "", mass_bins_, 0, 30);
            }

            string TwoD_plots = "2D_Plots";

                // Plots for only gen jets (blep,bhad, wja,wjb)
            for( const auto& genjet : GenJetNames_ ){
                book<TH2D>(TwoD_plots, std::string("3J_")+genjet+std::string("_Eta_vs_")+genjet+std::string("_Pt"), "", 200, 0., 1000., 100, 0., 5.0);
            }

                // Plots for gen tops (thad, tlep)
            for( const auto& gentop : GenTopNames_ ){
                    //Pt
                book<TH2D>(pt_plots, std::string("Pt_")+gentop+"_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                book<TH1F>(pt_plots, std::string("Pt_")+gentop, "", pt_bins_, pt_min_, pt_max_);

                    //Eta
                book<TH2D>(eta_plots, std::string("Eta_")+gentop+"_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                book<TH1F>(eta_plots, std::string("Eta_")+gentop, "", eta_bins_, eta_min_, eta_max_);

                    //Costh
                book<TH2D>(costh_plots, std::string("Costh_")+gentop+"_vs_Mtt", "", mass_bins_, mass_min_, mass_max_, costh_bins_, costh_min_, costh_max_);
                book<TH1F>(costh_plots, std::string("Costh_")+gentop, "", costh_bins_, costh_min_, costh_max_);

                    //Mass
                book<TH1F>(mass_plots, std::string("Mass_")+gentop, "", mass_bins_, 125., 250.);
            }
            
                // ttbar plots
            book<TH1F>(pt_plots, "Pt_ttbar", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(eta_plots, "Eta_ttbar", "", eta_bins_, eta_min_, eta_max_);
            book<TH1F>(costh_plots, "Costh_ttbar", "", costh_bins_, costh_min_, costh_max_);
            book<TH1F>(mass_plots, "Mass_ttbar", "", mass_bins_, mass_min_, mass_max_);


            string deltphi_plots = "DeltPhi_Plots";
            // Delta Phi between bs
            book<TH1F>(deltphi_plots, "DeltPhi_BHadBLep_4J", "", deltphi_bins_, deltphi_min_, deltphi_max_);

            string sys_plots = "Sys_Plots";
            //System Plots
            book<TH1F>(sys_plots, "nJets", "", 15, 0.5, 15.5);

            // plots for different DR values
            for( int i = 0; i < 2; i++ ){
                // all numbers of jets allowed
                book<TH1F>(DRnames_[i], "Gen_Had_Resolved_vs_mttbar", "", nbins, mass_min_, mass_max_); // 3 hadronic jets resolved
                book<TH1F>(DRnames_[i], "Gen_Merged_BHadWJet_vs_mttbar", "", nbins, mass_min_, mass_max_); // BHad and WJa or WJb merged
                book<TH1F>(DRnames_[i], "Gen_Merged_WJets_vs_mttbar", "", nbins, mass_min_, mass_max_); // jets from W merged
                book<TH1F>(DRnames_[i], "Gen_Merged_THad_vs_mttbar", "", nbins, mass_min_, mass_max_); // all jets on hadronic side merged
                //                book<TH1F>(DRnames_[i], "Gen_Non_Reconstructable_vs_mttbar", "", nbins, mass_min_, mass_max_); // 1 or more partons not matched
                book<TH1F>(DRnames_[i], "Gen_Had_Resolved_THadpt", "", pt_bins_, pt_min_, pt_max_); // 3 hadronic jets resolved
                book<TH1F>(DRnames_[i], "Gen_Merged_BHadWJet_THadpt", "", pt_bins_, pt_min_, pt_max_); // BHad and WJa or WJb merged
                book<TH1F>(DRnames_[i], "Gen_Merged_WJets_THadpt", "", pt_bins_, pt_min_, pt_max_); // jets from W merged
                book<TH1F>(DRnames_[i], "Gen_Merged_THad_THadpt", "", pt_bins_, pt_min_,pt_max_); // all jets on hadronic side merged
            }

            string had_comp = "Had_comp";

            // all jets
            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_Mtt", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_Mtt", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_Mtt", "", mass_bins_, mass_min_, mass_max_);

            book<TH2D>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Merged_THad_DRP8_THadpt_vs_Mttbar", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);

            // 3 jets
            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_3J", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_THadpt_3J", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_3J", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_Mtt_3J", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_Mtt_3J", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_Mtt_3J", "", mass_bins_, mass_min_, mass_max_);

            book<TH2D>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar_3J", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Merged_THad_DRP8_THadpt_vs_Mttbar_3J", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar_3J", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);

            // 4 jets
            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_4J", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_THadpt_4J", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_4J", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_Mtt_4J", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_Mtt_4J", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_Mtt_4J", "", mass_bins_, mass_min_, mass_max_);

            book<TH2D>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar_4J", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Merged_THad_DRP8_THadpt_vs_Mttbar_4J", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar_4J", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);

            // 5+ jets
            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_5PJ", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_THadpt_5PJ", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_5PJ", "", pt_bins_, pt_min_, pt_max_);

            book<TH1F>(had_comp, "Gen_Had_Resolved_DRP4_Mtt_5PJ", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Merged_THad_DRP8_Mtt_5PJ", "", mass_bins_, mass_min_, mass_max_);
            book<TH1F>(had_comp, "Gen_Partially_Merged_DRP8_Mtt_5PJ", "", mass_bins_, mass_min_, mass_max_);

            book<TH2D>(had_comp, "Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar_5PJ", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Merged_THad_DRP8_THadpt_vs_Mttbar_5PJ", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
            book<TH2D>(had_comp, "Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar_5PJ", "", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);

            //    // diff DR for partially merged
            //            // all jets
            //            book<TH1F>(had_comp, "Gen_Merged_BHadWJet_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_WJets_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            //
            //            // 3 jets
            //            book<TH1F>(had_comp, "Gen_Merged_BHadWJet_DRP4_THadpt_3J", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_WJets_DRP4_THadpt_3J", "", pt_bins_, pt_min_, pt_max_);
            //
            //            // 4 jets
            //            book<TH1F>(had_comp, "Gen_Merged_BHadWJet_DRP4_THadpt_4J", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_WJets_DRP4_THadpt_4J", "", pt_bins_, pt_min_, pt_max_);
            //
            //            // 5+ jets
            //            book<TH1F>(had_comp, "Gen_Merged_BHadWJet_DRP4_THadpt_5PJ", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_WJets_DRP4_THadpt_5PJ", "", pt_bins_, pt_min_, pt_max_);

            // one parton lost from acceptance
            book<TH1F>(had_comp, "Three_Parton_Gen_Had_Resolved_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_Gen_Merged_THad_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_Gen_Merged_BHadWJet_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_Gen_Merged_WJets_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);

            // blep missing
            book<TH1F>(had_comp, "Three_Parton_BLep_Missing_Gen_Had_Resolved_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_BLep_Missing_Gen_Merged_THad_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_BLep_Missing_Gen_Merged_BHadWJet_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_BLep_Missing_Gen_Merged_WJets_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);

            // bhad missing
            book<TH1F>(had_comp, "Three_Parton_BHad_Missing_Gen_Had_Resolved_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_BHad_Missing_Gen_Merged_THad_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_BHad_Missing_Gen_Merged_BHadWJet_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_BHad_Missing_Gen_Merged_WJets_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);

            // wja missing
            book<TH1F>(had_comp, "Three_Parton_WJa_Missing_Gen_Had_Resolved_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_WJa_Missing_Gen_Merged_THad_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_WJa_Missing_Gen_Merged_BHadWJet_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_WJa_Missing_Gen_Merged_WJets_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);

            // wjb missing
            book<TH1F>(had_comp, "Three_Parton_WJb_Missing_Gen_Had_Resolved_DRP4_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_WJb_Missing_Gen_Merged_THad_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_WJb_Missing_Gen_Merged_BHadWJet_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);
            book<TH1F>(had_comp, "Three_Parton_WJb_Missing_Gen_Merged_WJets_DRP8_THadpt", "", pt_bins_, pt_min_, pt_max_);


            //            book<TH1F>(had_comp, "Gen_Had_Resolved_THadpt", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_BHadWJet_THadpt", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_WJets_THadpt", "", pt_bins_, pt_min_, pt_max_);
            //            book<TH1F>(had_comp, "Gen_Merged_THad_THadpt", "", pt_bins_, pt_min_, pt_max_);

            Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }

        //This method is called once every file, contains the event loop
        ///run your proper analysis here
        virtual void analyze()
        {
            Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;
        cout << "jet pt min: " << cut_jet_ptmin_ << ", max eta: " << cut_jet_etamax_ << ", lead jet pt: " << cut_leadjet_ptmin_ << endl;

            URStreamer event(tree_);

            while(event.next() /*&& evt_idx_ < 30*/)
            {
                evt_idx_++;
                if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

                auto gen_dir = histos_.find("Gen_Plots");
                auto pt_dir = histos_.find("Pt_Plots");
                auto eta_dir = histos_.find("Eta_Plots");
                auto costh_dir = histos_.find("Costh_Plots");
                auto mass_dir = histos_.find("Mass_Plots");
                auto dphi_dir = histos_.find("DeltPhi_Plots");
                auto sys_dir = histos_.find("Sys_Plots");
                auto hc_dir = histos_.find("Had_comp");
                auto twoD_dir = histos_.find("2D_Plots");

                //generator selection
                bool selection = genp_selector_.select(event);
                if( !selection ){
                    Logger::log().debug() << "event has no selection " << endl;
                    continue;
                }
                tracker_.track("gen selection");
                GenTTBar &ttbar = genp_selector_.ttbar_system();


                //            njets_ = 0;
                if( object_selector_.select(event) ) njets_ = object_selector_.clean_jets().size();
                //            if( njets_ < 3 ) continue;

                if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ) continue;
                sys_dir->second["nJets"].fill(njets_);

                vector<GenObject*> GenObjs;
                vector<GenObject*> GenJets;
                vector<GenTop*> GenTops;
                //initialize gen partons
                GenObject* lepton = ttbar.lepton();
                GenObject* BLep = ttbar.lep_b();
                GenObject* BHad = ttbar.had_b();
                GenObject* WJa = (ttbar.had_W()->first->Pt() > ttbar.had_W()->second->Pt() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
                GenObject* WJb = (ttbar.had_W()->first->Pt() > ttbar.had_W()->second->Pt() ) ? ttbar.had_W()->second : ttbar.had_W()->first;
                GenTop* THad = ttbar.had_top();
                GenTop* TLep = ttbar.lep_top();

                //GenTop* leading_top = (THad->Pt() > TLep->Pt()) ? THad : TLep;
                //GenTop* subleading_top = (THad->Pt() > TLep->Pt()) ? TLep : THad;

                if( !(BHad && BLep && WJa && WJb && THad && TLep && lepton) ){
                    cout << "not all partons" << endl;
                    throw 42;
                }

                GenObjs.push_back(lepton);
                GenObjs.push_back(BLep);
                GenObjs.push_back(BHad);
                GenObjs.push_back(WJa);
                GenObjs.push_back(WJb);
                GenJets.push_back(BLep);
                GenJets.push_back(BHad);
                GenJets.push_back(WJa);
                GenJets.push_back(WJb);
                GenTops.push_back(THad);
                GenTops.push_back(TLep);

                if( !genp_selector_.is_in_acceptance(GenTTBar::DecayType::SEMILEP) ){ // not all partons in acceptance
                    if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){ // one parton falls outside the acceptance (pt=cut_jet_ptmin_, eta=cut_jet_etamax_)
                        /// any of the partons can fall outside
                        if( ttbar.resolved_had_partons(0.4) ){ // hadronic patons resolved at DR=0.4
                            hc_dir->second["Three_Parton_Gen_Had_Resolved_DRP4_THadpt"].fill(THad->Pt());
                        }


                        else{ // partons not resolved at DR = 0.4
                            if( ttbar.merged_had_partons(0.8) ){ // all 3 had partons merged at DR = 0.8
                                hc_dir->second["Three_Parton_Gen_Merged_THad_DRP8_THadpt"].fill(THad->Pt());
                            }

                            else{ // all 3 had partons not merged at DR = 0.8
                                if( ttbar.only_merged_wjawjb(0.8) ){// W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_Gen_Merged_WJets_DRP8_THadpt"].fill(THad->Pt());
                                }

                                if( ttbar.only_merged_bhadw_partons(0.8) ){ // BHad and one W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_Gen_Merged_BHadWJet_DRP8_THadpt"].fill(THad->Pt());
                                }
                            }
                        }
                    }

                    if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) && !ttbar.is_bhad_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){ // only bhad parton falls outside the acceptance (pt=cut_jet_ptmin_, eta=cut_jet_etamax_)
                        if( ttbar.resolved_had_partons(0.4) ){ // hadronic partons resolved at DR = 0.4
                            hc_dir->second["Three_Parton_BHad_Missing_Gen_Had_Resolved_DRP4_THadpt"].fill(THad->Pt());
                        }

                        else{ // partons not resolved at DR = 0.4
                            if( ttbar.merged_had_partons(0.8) ){// all 3 had partons merged at DR = 0.8
                                hc_dir->second["Three_Parton_BHad_Missing_Gen_Merged_THad_DRP8_THadpt"].fill(THad->Pt());
                            }

                            else{ // all 3 had partons not merged at DR = 0.8
                                if( ttbar.only_merged_wjawjb(0.8) ){ // W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_BHad_Missing_Gen_Merged_WJets_DRP8_THadpt"].fill(THad->Pt());
                                }

                                if( ttbar.only_merged_bhadw_partons(0.8) ){ // BHad and one W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_BHad_Missing_Gen_Merged_BHadWJet_DRP8_THadpt"].fill(THad->Pt());
                                }
                            }
                        }
                    }

                    if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) && !ttbar.is_blep_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){ // only blep parton falls outside the acceptance (pt=cut_jet_ptmin_, eta=cut_jet_etamax_)
                        if( ttbar.resolved_had_partons(0.4) ){ // hadronic partons resolved at DR = 0.4
                            hc_dir->second["Three_Parton_BLep_Missing_Gen_Had_Resolved_DRP4_THadpt"].fill(THad->Pt());
                        }

                        else{ // partons not resolved at DR = 0.4
                            if( ttbar.merged_had_partons(0.8) ){// all 3 had partons merged at DR = 0.8
                                hc_dir->second["Three_Parton_BLep_Missing_Gen_Merged_THad_DRP8_THadpt"].fill(THad->Pt());
                            }

                            else{ // all 3 had partons not merged at DR = 0.8
                                if( ttbar.only_merged_wjawjb(0.8) ){ // W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_BLep_Missing_Gen_Merged_WJets_DRP8_THadpt"].fill(THad->Pt());
                                }

                                if( ttbar.only_merged_bhadw_partons(0.8) ){ // BHad and one W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_BLep_Missing_Gen_Merged_BHadWJet_DRP8_THadpt"].fill(THad->Pt());
                                }
                            }
                        }
                    }

                    if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) && !ttbar.is_wja_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){ // only wja parton falls outside the acceptance (pt=cut_jet_ptmin_, eta=cut_jet_etamax_)
                        if( ttbar.resolved_had_partons(0.4) ){ // hadronic partons resolved at DR = 0.4
                            hc_dir->second["Three_Parton_WJa_Missing_Gen_Had_Resolved_DRP4_THadpt"].fill(THad->Pt());
                        }

                        else{ // partons not resolved at DR = 0.4
                            if( ttbar.merged_had_partons(0.8) ){// all 3 had partons merged at DR = 0.8
                                hc_dir->second["Three_Parton_WJa_Missing_Gen_Merged_THad_DRP8_THadpt"].fill(THad->Pt());
                            }

                            else{ // all 3 had partons not merged at DR = 0.8
                                if( ttbar.only_merged_wjawjb(0.8) ){ // W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_WJa_Missing_Gen_Merged_WJets_DRP8_THadpt"].fill(THad->Pt());
                                }

                                if( ttbar.only_merged_bhadw_partons(0.8) ){ // BHad and one W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_WJa_Missing_Gen_Merged_BHadWJet_DRP8_THadpt"].fill(THad->Pt());
                                }
                            }
                        }
                    }

                    if( ttbar.three_partons_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) && !ttbar.is_wjb_in_acceptance(cut_jet_ptmin_, cut_jet_etamax_, cut_leadjet_ptmin_) ){ // only wjb parton falls outside the acceptance (pt=cut_jet_ptmin_, eta=cut_jet_etamax_)
                        if( ttbar.resolved_had_partons(0.4) ){ // hadronic partons resolved at DR = 0.4
                            hc_dir->second["Three_Parton_WJb_Missing_Gen_Had_Resolved_DRP4_THadpt"].fill(THad->Pt());
                        }

                        else{ // partons not resolved at DR = 0.4
                            if( ttbar.merged_had_partons(0.8) ){// all 3 had partons merged at DR = 0.8
                                hc_dir->second["Three_Parton_WJb_Missing_Gen_Merged_THad_DRP8_THadpt"].fill(THad->Pt());
                            }

                            else{ // all 3 had partons not merged at DR = 0.8
                                if( ttbar.only_merged_wjawjb(0.8) ){ // W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_WJb_Missing_Gen_Merged_WJets_DRP8_THadpt"].fill(THad->Pt());
                                }

                                if( ttbar.only_merged_bhadw_partons(0.8) ){ // BHad and one W partons merged at DR = 0.8
                                    hc_dir->second["Three_Parton_WJb_Missing_Gen_Merged_BHadWJet_DRP8_THadpt"].fill(THad->Pt());
                                }
                            }
                        }
                    }

                    continue;
                } // end of not all partons in acceptance

                //            njets_ = object_selector_.clean_jets().size();

                //			if( njets_ < 3 ) continue;


                // Kinematic variables
                    // Gen object hists
                int i = 0;
                for( auto& genobj : GenObjs ){
                        //Pt
                    pt_dir->second[std::string("Pt_")+GenObjNames_[i]+"_vs_Mtt"].fill(ttbar.M(), genobj->Pt());
                    pt_dir->second[std::string("Pt_")+GenObjNames_[i]].fill(genobj->Pt());

                        //Eta
                    eta_dir->second[std::string("Eta_")+GenObjNames_[i]+"_vs_Mtt"].fill(ttbar.M(), genobj->Eta());
                    eta_dir->second[std::string("Eta_")+GenObjNames_[i]].fill(genobj->Eta());

                        //Costh
                    costh_dir->second[std::string("Costh_")+GenObjNames_[i]+"_vs_Mtt"].fill(ttbar.M(), genobj->CosTheta());
                    costh_dir->second[std::string("Costh_")+GenObjNames_[i]].fill(genobj->CosTheta());

                        //Mass
                    mass_dir->second[std::string("Mass_")+GenObjNames_[i]].fill(genobj->M());

                    ++i;
                }

                    // 2D Gen jet hists
                if( njets_ == 3 ){
                    int j = 0;
                    for( auto& genjet : GenJets ){
                        twoD_dir->second[std::string("3J_")+GenJetNames_[j]+std::string("_Eta_vs_")+GenJetNames_[j]+std::string("_Pt")].fill(genjet->Pt(), Abs(genjet->Eta()));
                        //cout << GenJetNames_[j] << " Pt: " << genjet->Pt() << endl;
                        //cout << GenJetNames_[j] << " Eta: " << Abs(genjet->Eta()) << endl;
                        ++j;
                    }
                }


                hyp::TTbar gen_ttang(ttbar);
                auto gen_ttcm = gen_ttang.to_CM();
                double gen_thad_cth = gen_ttang.unit3D().Dot(gen_ttcm.thad().unit3D());// thad costheta* value
                double gen_tlep_cth = gen_ttang.unit3D().Dot(gen_ttcm.tlep().unit3D());


                    // Gen top hists
                int j = 0;
                for( auto& gentop : GenTops ){
                        //Pt
                    pt_dir->second[std::string("Pt_")+GenTopNames_[j]+"_vs_Mtt"].fill(ttbar.M(), gentop->Pt());
                    pt_dir->second[std::string("Pt_")+GenTopNames_[j]].fill(gentop->Pt());

                        //Eta
                    eta_dir->second[std::string("Eta_")+GenTopNames_[j]+"_vs_Mtt"].fill(ttbar.M(), gentop->Eta());
                    eta_dir->second[std::string("Eta_")+GenTopNames_[j]].fill(gentop->Eta());

                        //Costh
                    if( strcmp( GenTopNames_[j], "THad" ) == 0 ){
                        costh_dir->second[std::string("Costh_")+GenTopNames_[j]+"_vs_Mtt"].fill(ttbar.M(), gen_thad_cth);
                        costh_dir->second[std::string("Costh_")+GenTopNames_[j]].fill(gen_thad_cth);
                    }
                    if( strcmp( GenTopNames_[j], "TLep" ) == 0 ){
                        costh_dir->second[std::string("Costh_")+GenTopNames_[j]+"_vs_Mtt"].fill(ttbar.M(), gen_tlep_cth);
                        costh_dir->second[std::string("Costh_")+GenTopNames_[j]].fill(gen_tlep_cth);
                    }

                        //Mass
                    mass_dir->second[std::string("Mass_")+GenTopNames_[j]].fill(gentop->M());

                    ++j;
                }

                    // ttbar hists
                pt_dir->second["Pt_ttbar"].fill(ttbar.Pt());
                eta_dir->second["Eta_ttbar"].fill(ttbar.Eta());
                costh_dir->second["Costh_ttbar"].fill(ttbar.CosTheta());
                mass_dir->second["Mass_ttbar"].fill(ttbar.M());


                //initialize kinematic vars
                //			double DR_LepBHad = -1;
                //			double DR_LepBLep = -1;
                //			double DR_LepWJa = -1;
                //			double DR_LepWJb = -1;
                //			double DR_BHadBLep = -1;
                //			double DR_BHadWJa = -1;
                //			double DR_BHadWJb = -1;
                //			double DR_BLepWJa = -1;
                //			double DR_BLepWJb = -1;
                //			double DR_WJaWJb = -1;


                if( njets_ == 4 ) dphi_dir->second["DeltPhi_BHadBLep_4J"].fill(BHad->Phi()-BLep->Phi());

                for( int i = 0; i < 2; i++ ){

                    auto DR_dir = histos_.find(DRnames_[i]);

                    /// parton level hadronic event categories for DR = 0.4, 0.8
                    /// all 3 had jets merged
                    if( ttbar.merged_had_partons(DR_[i]) ){
                        DR_dir->second["Gen_Merged_THad_vs_mttbar"].fill(ttbar.M());
                        DR_dir->second["Gen_Merged_THad_THadpt"].fill(THad->Pt());
                    }

                    /// W jets merged
                    if( ttbar.only_merged_wjawjb(DR_[i]) ){
                        DR_dir->second["Gen_Merged_WJets_vs_mttbar"].fill(ttbar.M());
                        DR_dir->second["Gen_Merged_WJets_THadpt"].fill(THad->Pt());
                    }

                    ///BHad and one W jet merged
                    if( ttbar.only_merged_bhadw_partons(DR_[i]) ){
                        DR_dir->second["Gen_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                        DR_dir->second["Gen_Merged_BHadWJet_THadpt"].fill(THad->Pt());
                    }

                    //		    /// BHad and WJa merged
                    //			    if( ttbar.merged_bhadwja_partons(DR_[i]) ){
                    //				    DR_dir->second["Gen_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                    //				    DR_dir->second["Gen_Merged_BHadWJet_THadpt"].fill(THad->Pt());
                    //				}
                    //				
                    //		    /// BHad and WJb merged
                    //			    if( ttbar.merged_bhadwjb_partons(DR_[i]) ){
                    //				    DR_dir->second["Gen_Merged_BHadWJet_vs_mttbar"].fill(ttbar.M());
                    //				    DR_dir->second["Gen_Merged_BHadWJet_THadpt"].fill(THad->Pt());
                    //				}

                    /// hadronic jets resolved
                    if( ttbar.resolved_had_partons(DR_[i]) ){
                        DR_dir->second["Gen_Had_Resolved_vs_mttbar"].fill(ttbar.M());
                        DR_dir->second["Gen_Had_Resolved_THadpt"].fill(THad->Pt());
                    }
                }


                /// hadronic jets resolved at DR = 0.4
                if( ttbar.resolved_had_partons(0.4) ){
                    hc_dir->second["Gen_Had_Resolved_DRP4_THadpt"].fill(THad->Pt());
                    hc_dir->second["Gen_Had_Resolved_DRP4_Mtt"].fill(ttbar.M());
                    hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar"].fill(ttbar.M(), THad->Pt());
                    if( njets_ == 3 ){
                        hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_3J"].fill(THad->Pt());
                        hc_dir->second["Gen_Had_Resolved_DRP4_Mtt_3J"].fill(ttbar.M());
                        hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar_3J"].fill(ttbar.M(), THad->Pt());
                    }
                    if( njets_ == 4 ){
                        hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_4J"].fill(THad->Pt());
                        hc_dir->second["Gen_Had_Resolved_DRP4_Mtt_4J"].fill(ttbar.M());
                        hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar_4J"].fill(ttbar.M(), THad->Pt());
                    }
                    if( njets_ > 4 ){
                        hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_5PJ"].fill(THad->Pt());
                        hc_dir->second["Gen_Had_Resolved_DRP4_Mtt_5PJ"].fill(ttbar.M());
                        hc_dir->second["Gen_Had_Resolved_DRP4_THadpt_vs_Mttbar_5PJ"].fill(ttbar.M(), THad->Pt());
                    }
                }

                /// jets not resolved at DR = 0.4
                else if( !(ttbar.resolved_had_partons(0.4)) ){
                    /// all 3 had jets merged at DR = 0.8
                    if( ttbar.merged_had_partons(0.8) ){
                        hc_dir->second["Gen_Merged_THad_DRP8_THadpt"].fill(THad->Pt());
                        hc_dir->second["Gen_Merged_THad_DRP8_Mtt"].fill(ttbar.M());
                        hc_dir->second["Gen_Merged_THad_DRP8_THadpt_vs_Mttbar"].fill(ttbar.M(), THad->Pt());
                        if( njets_ == 3 ){
                            hc_dir->second["Gen_Merged_THad_DRP8_THadpt_3J"].fill(THad->Pt());
                            hc_dir->second["Gen_Merged_THad_DRP8_Mtt_3J"].fill(ttbar.M());
                            hc_dir->second["Gen_Merged_THad_DRP8_THadpt_vs_Mttbar_3J"].fill(ttbar.M(), THad->Pt());
                        }
                        if( njets_ == 4 ){
                            hc_dir->second["Gen_Merged_THad_DRP8_THadpt_4J"].fill(THad->Pt());
                            hc_dir->second["Gen_Merged_THad_DRP8_Mtt_4J"].fill(ttbar.M());
                            hc_dir->second["Gen_Merged_THad_DRP8_THadpt_vs_Mttbar_4J"].fill(ttbar.M(), THad->Pt());
                        }
                        if( njets_ > 4 ){
                            hc_dir->second["Gen_Merged_THad_DRP8_THadpt_5PJ"].fill(THad->Pt());
                            hc_dir->second["Gen_Merged_THad_DRP8_Mtt_5PJ"].fill(ttbar.M());
                            hc_dir->second["Gen_Merged_THad_DRP8_THadpt_vs_Mttbar_5PJ"].fill(ttbar.M(), THad->Pt());
                        }
                    }

                    /// jets not all merged at DR = 0.8
                    else if( !(ttbar.merged_had_partons(0.8)) ){
                        hc_dir->second["Gen_Partially_Merged_DRP8_THadpt"].fill(THad->Pt());
                        hc_dir->second["Gen_Partially_Merged_DRP8_Mtt"].fill(ttbar.M());
                        hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar"].fill(ttbar.M(), THad->Pt());
                        if( njets_ == 3 ){
                            hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_3J"].fill(THad->Pt());
                            hc_dir->second["Gen_Partially_Merged_DRP8_Mtt_3J"].fill(ttbar.M());
                            hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar_3J"].fill(ttbar.M(), THad->Pt());
                        }
                        if( njets_ == 4 ){
                            hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_4J"].fill(THad->Pt());
                            hc_dir->second["Gen_Partially_Merged_DRP8_Mtt_4J"].fill(ttbar.M());
                            hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar_4J"].fill(ttbar.M(), THad->Pt());
                        }
                        if( njets_ > 4 ){
                            hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_5PJ"].fill(THad->Pt());
                            hc_dir->second["Gen_Partially_Merged_DRP8_Mtt_5PJ"].fill(ttbar.M());
                            hc_dir->second["Gen_Partially_Merged_DRP8_THadpt_vs_Mttbar_5PJ"].fill(ttbar.M(), THad->Pt());
                        }
                    }
                } // end of partons not resolved

                //			if( njets_ < 3 ) continue;

                for( int i=0; i<2; i++ ){
                        // Lep and BLep
                    if( ttbar.only_merged_lepblep(DR_[i]) ){ // only lep and blep merged
                        gen_dir->second[std::string("DR_LepBLep_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_lepblep(DR_[i]) ){ // lep and blep not merged
                        gen_dir->second[std::string("DR_LepBLep_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // Lep and BHad
                    if( ttbar.only_merged_lepbhad(DR_[i]) ){ // only lep and bhad merged
                        gen_dir->second[std::string("DR_LepBHad_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_lepbhad(DR_[i]) ){ // lep and bhad not merged
                        gen_dir->second[std::string("DR_LepBHad_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // Lep and WJa
                    if( ttbar.only_merged_lepwja(DR_[i]) ){ // only lep and wja merged
                        gen_dir->second[std::string("DR_LepWJa_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_lepwja(DR_[i]) ){ // lep and wja not merged
                        gen_dir->second[std::string("DR_LepWJa_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // Lep and WJb
                    if( ttbar.only_merged_lepwjb(DR_[i]) ){ // only lep and wjb merged
                        gen_dir->second[std::string("DR_LepWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_lepwjb(DR_[i]) ){ // lep and wjb not merged
                        gen_dir->second[std::string("DR_LepWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // BLep and BHad
                    if( ttbar.only_merged_blepbhad(DR_[i]) ){ // only blep and bhad merged
                        gen_dir->second[std::string("DR_BHadBLep_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_blepbhad(DR_[i]) ){ // blep and bhad not merged
                        gen_dir->second[std::string("DR_BHadBLep_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // BLep and WJa
                    if( ttbar.only_merged_blepwja(DR_[i]) ){ // only blep and wja merged
                        gen_dir->second[std::string("DR_BLepWJa_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_blepwja(DR_[i]) ){ // blep and wja not merged
                        gen_dir->second[std::string("DR_BLepWJa_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // BLep and WJb
                    if( ttbar.only_merged_blepwjb(DR_[i]) ){ // only blep and wjb merged
                        gen_dir->second[std::string("DR_BLepWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_blepwjb(DR_[i]) ){ // blep and wjb not merged
                        gen_dir->second[std::string("DR_BLepWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // BHad and WJa
                    if( ttbar.only_merged_bhadwja(DR_[i]) ){ // only bhad and wja merged
                        gen_dir->second[std::string("DR_BHadWJa_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_bhadwja(DR_[i]) ){ // bhad and wja not merged
                        gen_dir->second[std::string("DR_BHadWJa_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // BHad and WJb
                    if( ttbar.only_merged_bhadwjb(DR_[i]) ){ // only bhad and wjb merged
                        gen_dir->second[std::string("DR_BHadWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_bhadwjb(DR_[i]) ){ // bhad and wjb not merged
                        gen_dir->second[std::string("DR_BHadWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }

                        // WJa and WJb
                    if( ttbar.only_merged_wjawjb(DR_[i]) ){ // only wja and wjb merged
                        gen_dir->second[std::string("DR_WJaWJb_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                    else if( !ttbar.only_merged_wjawjb(DR_[i]) ){ // wja and wjb not merged
                        gen_dir->second[std::string("DR_WJaWJb_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                    }
                } // end of loop to make merged fractions hists

                //
                //			//find min value of DR for THad objects
                //			list<double> THad_DR_list;
                //			list<double>::iterator THad_DR_it;
                //			if( !(BHad == 0) && !(WJa == 0) ) THad_DR_list.push_back(DR_BHadWJa);
                //			if( !(BHad == 0) && !(WJb == 0) ) THad_DR_list.push_back(DR_BHadWJb);
                //			if( !(WJa == 0) && !(WJb == 0) ) THad_DR_list.push_back(DR_WJaWJb);
                //			double THad_DRmin = 1e2;
                //			double THad_DRmax = 0;
                //			for( THad_DR_it = THad_DR_list.begin(); THad_DR_it != THad_DR_list.end(); ++THad_DR_it ){
                //				if( *THad_DR_it < THad_DRmin ){
                //					THad_DRmin = *THad_DR_it;
                //				}
                //				if( *THad_DR_it > THad_DRmax ){
                //					THad_DRmax = *THad_DR_it;
                //				}
                //			}
                //			if( THad_DRmin < 0.4 ){
                //				for( THad_DR_it = THad_DR_list.begin(); THad_DR_it != THad_DR_list.end(); ++THad_DR_it ){
                //					if( *THad_DR_it > THad_DRmax ){
                //						THad_DRmax = *THad_DR_it;
                //					}
                //				}
                //
                //			}
                //			if( !(THad_DRmax == 0) ){
                //				gen_dir->second["DRmax_THad_vs_mttbar"].fill(ttbar.M(), THad_DRmax);
                //				gen_dir->second["DRmax_THad_vs_ptTHad"].fill(THad->Pt(), THad_DRmax);
                //				gen_dir->second["DRmax_THad"].fill(THad_DRmax);
                //
                //				for( int i = 0; i < 2; i++ ){
                //					if( THad_DRmax < DR_[i] ){
                //						gen_dir->second[std::string("DRmax_THad_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
                //						gen_dir->second[std::string("DRmax_THad_l")+DRnames_[i]+"_vs_ptTHad"].fill(THad->Pt());
                //					}
                //					if( THad_DRmax >= DR_[i] ){
                //						gen_dir->second[std::string("DRmax_THad_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                //						gen_dir->second[std::string("DRmax_THad_g")+DRnames_[i]+"_vs_ptTHad"].fill(THad->Pt());
                //					}
                //				}
                //			}
                //			if( !(THad_DRmin == 1e2) ){
                //				gen_dir->second["DRmin_THad_vs_mttbar"].fill(ttbar.M(), THad_DRmin);
                //				gen_dir->second["DRmin_THad_vs_ptTHad"].fill(THad->Pt(), THad_DRmin);
                //				gen_dir->second["DRmin_THad"].fill(THad_DRmin);
                //
                //				for( int i = 0; i < 2; i++ ){
                //					if( THad_DRmin < DR_[i] ){
                //						gen_dir->second[std::string("DRmin_THad_l")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M()); 
                //						gen_dir->second[std::string("DRmin_THad_l")+DRnames_[i]+"_vs_ptTHad"].fill(THad->Pt());
                //					}
                //					if( THad_DRmin >= DR_[i] ){
                //						gen_dir->second[std::string("DRmin_THad_g")+DRnames_[i]+"_vs_mttbar"].fill(ttbar.M());
                //						gen_dir->second[std::string("DRmin_THad_g")+DRnames_[i]+"_vs_ptTHad"].fill(THad->Pt());
                //					}
                //				}
                //			}

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
    URDriver<gen_partons> test;
    int thing = test.run();
    return thing;
}
