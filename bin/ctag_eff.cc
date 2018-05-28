#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include <map>
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include <algorithm>
//#include "NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include <list>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include <boost/algorithm/string/predicate.hpp>
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "TMath.h"
#include "TRandom3.h"
#include "Analyses/URTTbar/interface/Permutation.h"
#include <set>
#include "Analyses/URTTbar/interface/IDMet.h"
#include "TUUID.h"   
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"

#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"

using namespace TMath;

static map<string, IDJet::BTag> available_bjet_id = {
    {"none"     , IDJet::BTag::NONE},
    {"csvMedium", IDJet::BTag::CSVMEDIUM},
    {"csvLoose" , IDJet::BTag::CSVLOOSE},
    {"csvTight" , IDJet::BTag::CSVTIGHT},

    {"ctagMedium", IDJet::BTag::CTAGMEDIUM},
    {"ctagLoose" , IDJet::BTag::CTAGLOOSE},
    {"ctagTight" , IDJet::BTag::CTAGTIGHT},

    {"DeepCSVMedium", IDJet::BTag::DEEPCSVMEDIUM},
    {"DeepCSVLoose" , IDJet::BTag::DEEPCSVLOOSE},
    {"DeepCSVTight" , IDJet::BTag::DEEPCSVTIGHT},

    //{"DeepctagMedium", IDJet::BTag::DEEPCTAGMEDIUM},
    //{"DeepctagLoose" , IDJet::BTag::DEEPCTAGLOOSE},
    //{"DeepctagTight" , IDJet::BTag::DEEPCTAGTIGHT},

    //{"mvaMedium", IDJet::BTag::MVAMEDIUM},
    //{"mvaLoose" , IDJet::BTag::MVALOOSE},
    //{"mvaTight" , IDJet::BTag::MVATIGHT},
};	

static map<string, std::function<bool(const Permutation &, const Permutation &)> > available_ordering = {
    {"full_discriminant", [](const Permutation &one, const Permutation &two) {return  one.Prob()  < two.Prob() ;}         },
    {"nu_chisq"         , [](const Permutation &one, const Permutation &two) {return  one.NuChisq() 	 < two.NuChisq() ;}	},
    {"nu_discriminant"	 , [](const Permutation &one, const Permutation &two) {return  one.NuDiscr() 	 < two.NuDiscr() ;}	},
    {"btag_discriminant", [](const Permutation &one, const Permutation &two) {return  one.BDiscr()  	 < two.BDiscr()  ;}	},
    {"mass_discriminant", [](const Permutation &one, const Permutation &two) {return  one.MassDiscr() < two.MassDiscr();} },
};


class ctag_eff : public AnalyzerBase
{
    public:
        enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
    private:
        //histograms and helpers
        map<string, map<string, RObject> > histos_;
        map<TTNaming, string> naming_;
        CutFlowTracker tracker_;

        //switches
        bool isData_, isTTbar_, pdfs_;

        //selectors and helpers
        TTGenParticleSelector genp_selector_;
        TTObjectSelector object_selector_;
        TTPermutator permutator_;
        TTGenMatcher matcher_;
        TTBarSolver solver_;
        float evt_weight_;
        TRandom3 randomizer_;// = TRandom3(98765);
        MCWeightProducer mc_weights_;
        BTagSFProducer btag_sf_;
        PDFuncertainty pdf_uncs_;

        vector<systematics::SysShifts> systematics_;
        string cut_ordering_;
        std::function<bool(const Permutation &, const Permutation &)> ordering_fcn_;
        map<string, std::function<bool(const IDJet*)> > working_points_;
        map<string, BTagSFProducer*> wp_SFs_;  

        //Scale factors
        LeptonSF electron_sf_, muon_sf_;

        IDJet::BTag cut_tight_b_=IDJet::BTag::NONE;
        IDJet::BTag cut_loose_b_=IDJet::BTag::NONE;



        unsigned long evt_idx_ = 0;
        vector<float> pt_binning_ = {25.,35.,50.,70.,200.};

    public:
        std::vector<systematics::SysShifts> get_systematics(std::string outname) {
            std::string sample = systematics::get_sample(outname);
            std::vector<systematics::SysShifts> full_sys = {
                systematics::SysShifts::NOSYS, systematics::SysShifts::JES_UP, systematics::SysShifts::JES_DW,
                systematics::SysShifts::JER_UP, systematics::SysShifts::JER_DW,
                systematics::SysShifts::PU_UP, systematics::SysShifts::PU_DW,
                systematics::SysShifts::BTAG_UP, systematics::SysShifts::BTAG_DW, //Systematic on B-Jet selection
                systematics::SysShifts::BTAG_B_UP, systematics::SysShifts::BTAG_B_DW, //systematics on B/C-Tagging
                systematics::SysShifts::BTAG_L_UP, systematics::SysShifts::BTAG_L_DW, //systematics on B/C-Tagging
                systematics::SysShifts::BTAG_C_UP, systematics::SysShifts::BTAG_C_DW, //systematics on B/C-Tagging
            };
            std::vector<systematics::SysShifts> nosys = {systematics::SysShifts::NOSYS};

            if(boost::starts_with(sample, "ttJets")) {
                if(sample.find("_") != std::string::npos) {
                    //sys shifted sample!
                    return nosys;
                }
                else return full_sys;
            }
            else if(!boost::starts_with(sample, "data")) {
                return full_sys;
            }
            else return nosys;
        }

        ctag_eff(const std::string output_filename):
            AnalyzerBase("ctag_eff", output_filename), 
            tracker_(),
            working_points_(),
            object_selector_(),
            permutator_(),
            matcher_(),
            solver_(true),
            mc_weights_(),
            evt_weight_(1.),
            electron_sf_("electron_sf", false),
            muon_sf_("muon_sf"),
            randomizer_(),
            btag_sf_("best_permutation.tightb", "best_permutation.looseb", 0.5, -1, -1),
            //pdf_uncs_(108) // value for 2017
            pdf_uncs_(254) // value for 2016
    {
        btag_sf_.ignore_partial_shifts();
        cut_tight_b_ = btag_sf_.tight_cut();
        cut_loose_b_ = btag_sf_.loose_cut();
        // cout << "bcuts " << cut_tight_b_ << " " << cut_loose_b_ << 
        //   " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;

        URParser &parser = URParser::instance();

        TUUID id;  
        randomizer_.SetSeed(id.Hash());   

        //set tracker
        tracker_.use_weight(&evt_weight_);
        object_selector_.set_tracker(&tracker_);
        permutator_.set_tracker(&tracker_);

        //find out which sample are we running on
        opts::variables_map &values = parser.values();
        string output_file = values["output"].as<std::string>();
        string sample = systematics::get_sample(output_file);
        isData_  = boost::starts_with(sample, "data");
        isTTbar_ = boost::starts_with(sample, "ttJets");
        Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
        systematics_ = get_systematics(output_file);    
        int nosys = values["nosys"].as<int>();
        int nopdf = values["nopdf"].as<int>();
        if(nosys == 1) {
            systematics_ = {systematics::SysShifts::NOSYS};
            Logger::log().info() << "DISABLING SYSTEMATICS!" << endl;
        }
        if(isData_) {
            if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
            else object_selector_.lepton_type(1);
        }
        pdfs_ = isTTbar_ && (systematics_.size() > 1) && (nopdf == 0);

        if(!isData_) mc_weights_.init(sample);

        //SET CUTS FROM CFG
        cut_ordering_ = parser.getCfgPar<string>("permutations.ordering");
        ordering_fcn_ = available_ordering[cut_ordering_];

        working_points_["notag"]     = [](const IDJet* jet) {return false;};
        working_points_["csvLoose"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVLOOSE);}; 
        working_points_["csvTight"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVTIGHT);};
        working_points_["csvMedium"] = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::CSVMEDIUM);};
        working_points_["ctagLoose"]  = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::CTAGLOOSE);};
        working_points_["ctagTight"]  = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::CTAGTIGHT);};
        working_points_["ctagMedium"] = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::CTAGMEDIUM);};
        //working_points_["cmvaLoose"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::MVALOOSE);};
        //working_points_["cmvaTight"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::MVATIGHT);};
        //working_points_["cmvaMedium"] = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::MVAMEDIUM);};
        working_points_["DeepCSVLoose"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::DEEPCSVLOOSE);};
        working_points_["DeepCSVMedium"] = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::DEEPCSVMEDIUM);};
        working_points_["DeepCSVTight"]  = [](const IDJet* jet) {return jet->BTagId(IDJet::BTag::DEEPCSVTIGHT);};
        //working_points_["DeepctagLoose"]  = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::DEEPCTAGLOOSE);};
        //working_points_["DeepctagTight"]  = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::DEEPCTAGTIGHT);};
        //working_points_["DeepctagMedium"] = [](const IDJet* jet) {return jet->CTagId(IDJet::BTag::DEEPCTAGMEDIUM);};

        //Get appropriate SFs for the probe working points
        DataFile wjet_efficiency(parser.getCfgPar<string>("general.wjets_efficiencies"));
        DataFile csv_sfs(parser.getCfgPar<string>("general.csv_sffile"));
        DataFile deepcsv_sfs(parser.getCfgPar<string>("general.deepcsv_sffile"));
        DataFile ctag_sfs(parser.getCfgPar<string>("general.ctag_sffile"));
        //DataFile deepctag_sfs(parser.getCfgPar<string>("general.deepctag_sffile"));
        //DataFile cmva_sfs(parser.getCfgPar<string>("general.cmva_sffile"));

        wp_SFs_["csvLoose" ] 	= new BTagSFProducer(csv_sfs, wjet_efficiency, IDJet::BTag::CSVLOOSE , IDJet::BTag::NONE, 0.5, -1, -1);
        wp_SFs_["csvTight" ] 	= new BTagSFProducer(csv_sfs, wjet_efficiency, IDJet::BTag::CSVTIGHT , IDJet::BTag::NONE, 0.5, -1, -1); 
        wp_SFs_["csvMedium"] 	= new BTagSFProducer(csv_sfs, wjet_efficiency, IDJet::BTag::CSVMEDIUM, IDJet::BTag::NONE, 0.5, -1, -1); 

        wp_SFs_["DeepCSVLoose" ] = new BTagSFProducer(deepcsv_sfs, wjet_efficiency, IDJet::BTag::DEEPCSVLOOSE , IDJet::BTag::NONE, 0.5, -1, -1);
        wp_SFs_["DeepCSVTight" ] = new BTagSFProducer(deepcsv_sfs, wjet_efficiency, IDJet::BTag::DEEPCSVTIGHT , IDJet::BTag::NONE, 0.5, -1, -1); 
        wp_SFs_["DeepCSVMedium"] = new BTagSFProducer(deepcsv_sfs, wjet_efficiency, IDJet::BTag::DEEPCSVMEDIUM, IDJet::BTag::NONE, 0.5, -1, -1); 

        //wp_SFs_["cmvaLoose" ] = new BTagSFProducer(cmva_sfs, wjet_efficiency, IDJet::BTag::MVALOOSE , IDJet::BTag::NONE, 0.5, -1, -1); 
        //wp_SFs_["cmvaTight" ] = new BTagSFProducer(cmva_sfs, wjet_efficiency, IDJet::BTag::MVATIGHT , IDJet::BTag::NONE, 0.5, -1, -1); 
        //wp_SFs_["cmvaMedium"] = new BTagSFProducer(cmva_sfs, wjet_efficiency, IDJet::BTag::MVAMEDIUM, IDJet::BTag::NONE, 0.5, -1, -1); 

        wp_SFs_["ctagLoose" ] = new BTagSFProducer(ctag_sfs, wjet_efficiency, IDJet::BTag::CTAGLOOSE , IDJet::BTag::NONE, 0.5, -1, 0.5); 
        wp_SFs_["ctagTight" ] = new BTagSFProducer(ctag_sfs, wjet_efficiency, IDJet::BTag::CTAGTIGHT , IDJet::BTag::NONE, 0.5, -1, 0.5); 
        wp_SFs_["ctagMedium"] = new BTagSFProducer(ctag_sfs, wjet_efficiency, IDJet::BTag::CTAGMEDIUM, IDJet::BTag::NONE, 0.5, -1, 0.5); 

        //wp_SFs_["DeepctagLoose" ] = new BTagSFProducer(deepctag_sfs, wjet_efficiency, IDJet::BTag::DEEPCTAGLOOSE , IDJet::BTag::NONE, 0.5, -1, 1); 
        //wp_SFs_["DeepctagTight" ] = new BTagSFProducer(deepctag_sfs, wjet_efficiency, IDJet::BTag::DEEPCTAGTIGHT , IDJet::BTag::NONE, 0.5, -1, 1); 
        //wp_SFs_["DeepctagMedium"] = new BTagSFProducer(deepctag_sfs, wjet_efficiency, IDJet::BTag::DEEPCTAGMEDIUM, IDJet::BTag::NONE, 0.5, -1, 1); 

        wp_SFs_["csvLoose" ]->ignore_general_shifts();
        wp_SFs_["csvTight" ]->ignore_general_shifts();
        wp_SFs_["csvMedium"]->ignore_general_shifts();

        wp_SFs_["ctagLoose" ]->ignore_general_shifts(); 
        wp_SFs_["ctagTight" ]->ignore_general_shifts(); 
        wp_SFs_["ctagMedium"]->ignore_general_shifts(); 

        //wp_SFs_["cmvaLoose" ]->ignore_general_shifts(); 
        //wp_SFs_["cmvaTight" ]->ignore_general_shifts(); 
        //wp_SFs_["cmvaMedium"]->ignore_general_shifts(); 

        wp_SFs_["DeepCSVLoose" ]->ignore_general_shifts();
        wp_SFs_["DeepCSVTight" ]->ignore_general_shifts(); 
        wp_SFs_["DeepCSVMedium"]->ignore_general_shifts(); 

        //wp_SFs_["DeepctagLoose" ]->ignore_general_shifts(); 
        //wp_SFs_["DeepctagTight" ]->ignore_general_shifts(); 
        //wp_SFs_["DeepctagMedium"]->ignore_general_shifts(); 
        // working_points_[] = [](const Jet* jet) {};

        naming_[TTNaming::RIGHT ] = "semilep_visible_right";
        naming_[TTNaming::RIGHT_THAD ] = "semilep_right_thad";
        naming_[TTNaming::RIGHT_WHAD ] = "semilep_right_whad";
        naming_[TTNaming::RIGHT_TLEP ] = "semilep_right_tlep";
        naming_[TTNaming::WRONG ] = "semilep_wrong" 			 ;
        naming_[TTNaming::OTHER ] = 	"other"              ;

        Logger::log().debug() << "Constructor completed" << std::endl;
    }

        ~ctag_eff() {
            wp_SFs_.erase(wp_SFs_.begin(), wp_SFs_.end());
        }

        // int get_bin(IDJet *j1, IDJet *j2) {
        // 	pt_binning_
        // }

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
        void book_combo_plots(string folder){
            book<TH1F>(folder, "mass_discriminant", "", 20,   0., 20.);

            book<TH2F>(folder, "flav_map", ";leading;subleading", 3, 0., 3., 3, 0., 3.);

            book<TH1F>(folder, "tmasshad", "", 50, 0., 500);			
            book<TH1F>(folder, "Wmasshad", "", 50, 0., 500);			
            book<TH1F>(folder, "tpthad", "", 50, 0., 500);			
            //book<TH1F>(folder, "evt_weight", "", 100, 0., 50.);
        }

        void fill_combo_plots(string folder, const Permutation &hyp){
            auto dir = histos_.find(folder);
            dir->second["mass_discriminant"].fill(hyp.MassDiscr(), evt_weight_);
            //dir->second["evt_weight"].fill(evt_weight_);
            float lflav = -1;
            switch(hyp.WJa()->hadronFlavour()) {
                case 5: lflav = 2.5; break;
                case 4: lflav = 1.5; break;
                default: lflav = 0.5; break;      
            }
            float sflav = -1;
            switch(hyp.WJb()->hadronFlavour()) {
                case 5:  sflav = 2.5; break;
                case 4:  sflav = 1.5; break;
                default: sflav = 0.5; break;      
            }
            dir->second["flav_map"].fill(lflav, sflav, evt_weight_);

            double whad_mass = hyp.WHad().M();
            double thad_mass = hyp.THad().M();
            dir->second["Wmasshad"].fill(whad_mass, evt_weight_);
            dir->second["tmasshad"].fill(thad_mass, evt_weight_);
            dir->second["tpthad"].fill(hyp.THad().Pt(), evt_weight_);
        }

        void book_pdf_plots(string folder) {
            pdf_uncs_.book_replicas<TH1F>(folder, "mass_discriminant", 20,   0., 20.);
        }
        void fill_pdf_plots(string folder, const Permutation &hyp, URStreamer& streamer) {
            pdf_uncs_.fill_replicas(folder, "mass_discriminant", hyp.MassDiscr(), evt_weight_, streamer);
        }

        void book_presel_plots(string folder) {
            book<TH1F>(folder, "nvtx_noweight", "", 100, 0., 100.);
            book<TH1F>(folder, "nvtx", "", 100, 0., 100.);
            book<TH1F>(folder, "rho", "", 100, 0., 100.);
            book<TH1F>(folder, "weight", "", 100, 0., 20.);
            book<TH1F>(folder, "lep_pt"   , ";p_{T}(l) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "lep_eta"  , ";#eta(l) (GeV)", 300, -3, 3);
            book<TH1F>(folder, "jets_pt"  , ";p_{T}(j) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "jets_eta" , ";#eta(j) (GeV)",  300, -3, 3);
            book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
            book<TH1F>(folder, "jets_CvsL" , "", 55, -1., 1.1);
            book<TH1F>(folder, "jets_CvsB" , "", 55, -1., 1.1);
            book<TH1F>(folder, "jets_CSV"  , "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_cMVA" , "", 55, -1., 1.1);
            //book<TH1F>(folder, "jets_DeepCSVCvsL" , "", 55,  0., 1.1);
            //book<TH1F>(folder, "jets_DeepCSVCvsB" , "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_DeepCSVb" , "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_DeepCSVl" , "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_DeepCSVbb", "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_DeepCSVc" , "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_DeepCSVcc", "", 55,  0., 1.1);
            book<TH1F>(folder, "jets_DeepCSVbD", "", 55,  0., 1.1);

            // book<TH1F>(folder, "jets_hflav_CvsL_B" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "jets_hflav_CvsL_C" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "jets_hflav_CvsL_S" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "jets_hflav_CvsL_L" , "", 55, -1., 1.1);
        }

        void fill_presel_plots(string folder, URStreamer &event) {
            auto dir = histos_.find(folder);
            dir->second["nvtx"].fill(event.vertexs().size(), evt_weight_);
            dir->second["nvtx_noweight"].fill(event.vertexs().size());
            dir->second["rho"].fill(event.rho().value(), evt_weight_);
            dir->second["weight"].fill(evt_weight_);
            dir->second["lep_pt"].fill(object_selector_.lepton()->Pt(), evt_weight_);
            dir->second["lep_eta"].fill(object_selector_.lepton()->Eta(), evt_weight_);
            for(IDJet* jet : object_selector_.clean_jets()) {
                dir->second["jets_pt"].fill(jet->Pt(), evt_weight_);
                dir->second["jets_eta"].fill(jet->Eta(), evt_weight_);
                dir->second["jets_CvsL"].fill(jet->CvsLtag(), evt_weight_);
                dir->second["jets_CvsB"].fill(jet->CvsBtag(), evt_weight_);
                dir->second["jets_CSV" ].fill(jet->csvIncl(), evt_weight_);
                dir->second["jets_cMVA"].fill(jet->CombinedMVA(), evt_weight_);


                //dir->second["jets_DeepCSVCvsL" ].fill(jet->DeepCSVCvsLtag(), evt_weight_);
                //dir->second["jets_DeepCSVCvsB" ].fill(jet->DeepCSVCvsBtag(), evt_weight_);

                dir->second["jets_DeepCSVb" ].fill(jet->DeepCSVProbB(), evt_weight_);
                dir->second["jets_DeepCSVl" ].fill(jet->DeepCSVProbUDSG(), evt_weight_);
                dir->second["jets_DeepCSVbb"].fill(jet->DeepCSVProbBB(), evt_weight_);
                dir->second["jets_DeepCSVc" ].fill(jet->DeepCSVProbC(), evt_weight_);
                dir->second["jets_DeepCSVcc"].fill(jet->DeepCSVProbCC(), evt_weight_);
                dir->second["jets_DeepCSVbD"].fill(jet->DeepCSVProbB()+jet->DeepCSVProbBB(), evt_weight_);
                // int hflav = fabs(jet->hadronFlavour());
                // int pflav = fabs(jet->partonFlavour());
                // string hstr;
                // if(hflav == 5) hstr="B";
                // else if(hflav == 4) hstr="C";
                // else if(pflav == 3) hstr="S";
                // else hstr="L"; 
                // dir->second["jets_hflav_CvsL_"+hstr].fill(jet->CvsLtag(), evt_weight_);

        //cout << "DeepCSV: "  << jet->DeepCSVProbB()+jet->DeepCSVProbBB() << endl;
        //cout << "DeepCSV CvsL: " << jet->DeepCSVCvsLtag() << endl;
        //cout << "DeepCSV CvsB: " << jet->DeepCSVCvsBtag() << endl;
            }
            dir->second["njets" ].fill(object_selector_.clean_jets().size(), evt_weight_);
        }

        void book_notag_plots(string folder){
            book<TH1F>(folder, "evt_weight", "", 100, 0., 50.);
            book<TH1F>(folder, "btag_sf", "", 50, 0.8, 1.4);
            book<TH1F>(folder, "btag_sf_B", "", 50, 0.8, 1.4);
            book<TH1F>(folder, "btag_sf_C", "", 50, 0.8, 1.4);
            book<TH1F>(folder, "btag_sf_L", "", 50, 0.8, 1.4);
            book<TH2F>(folder, "btag_sf_vs_pt_B", ";p_{T}(b) (GeV)", 100, 0., 500., 50, 0.8, 1.4);
            book<TH2F>(folder, "btag_sf_vs_pt_C", ";p_{T}(c) (GeV)", 100, 0., 500., 50, 0.8, 1.4);
            book<TH2F>(folder, "btag_sf_vs_pt_L", ";p_{T}(l) (GeV)", 100, 0., 500., 50, 0.8, 1.4);
            book<TH1F>(folder, "muon_sf", "", 50, 0.9, 1.1);
            book<TH2F>(folder, "muon_sf_vs_mu_pt", ";p_{T}(#mu) (GeV)", 500, 0., 500., 50, 0.9, 1.1);
            book<TH2F>(folder, "muon_sf_vs_mu_eta", ";#eta(#mu) (GeV)", 300, -3., 3., 50, 0.9, 1.1);
            book<TH1F>(folder, "njets"    , "", 50, 0., 50.);
            book<TH1F>(folder, "lep_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
            book<TH1F>(folder, "had_b_pt" , ";p_{T}(b) (GeV)", 100, 0., 500.);
            book<TH1F>(folder, "lep_pt"   , ";p_{T}(#ell) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "Whad_mass", ";m_{W}(had) (GeV)", 28, 0., 140.);
            book<TH1F>(folder, "Whad_DR"  , ";m_{W}(had) (GeV)", 100, 0., 10.);
            book<TH1F>(folder, "Whad_pt"  , ";m_{W}(had) (GeV)", 100, 0., 500.);
            book<TH1F>(folder, "thad_mass", ";m_{t}(had) (GeV)", 60, 100., 400.);
            book<TH1F>(folder, "thad_pt", ";p_{T}(t_{had}) (GeV)", 100, 0., 500.);
            book<TH2F>(folder, "Whad_jet_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);
            book<TH2F>(folder, "BJet_jet_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);
            book<TH2F>(folder, "leadB_leadW_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);
            book<TH2F>(folder, "subB_subW_pts" , ";lead pT; sublead pT", 50, 0., 500., 50, 0., 500.);

            //info by flavor    
            book<TH1F>(folder, "Wjets_hflav_CvsL_B" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Wjets_hflav_CvsL_C" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Wjets_hflav_CvsL_L" , "", 55, -1., 1.1);

            book<TH1F>(folder, "Wjets_hflav_CvsB_B" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Wjets_hflav_CvsB_C" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Wjets_hflav_CvsB_L" , "", 55, -1., 1.1);

            //DeepCSV Disc values
            //book<TH1F>(folder, "Wjets_hflav_DeepCSVCvsL_B" , "", 55, 0., 1.1);
            //book<TH1F>(folder, "Wjets_hflav_DeepCSVCvsL_C" , "", 55, 0., 1.1);
            //book<TH1F>(folder, "Wjets_hflav_DeepCSVCvsL_L" , "", 55, 0., 1.1);

            //book<TH1F>(folder, "Wjets_hflav_DeepCSVCvsB_B" , "", 55, 0., 1.1);
            //book<TH1F>(folder, "Wjets_hflav_DeepCSVCvsB_C" , "", 55, 0., 1.1);
            //book<TH1F>(folder, "Wjets_hflav_DeepCSVCvsB_L" , "", 55, 0., 1.1);

            book<TH1F>(folder, "Wjets_hflav_DeepCSVbD_B" , "", 55, 0., 1.1);
            book<TH1F>(folder, "Wjets_hflav_DeepCSVbD_C" , "", 55, 0., 1.1);
            book<TH1F>(folder, "Wjets_hflav_DeepCSVbD_L" , "", 55, 0., 1.1);

            book<TH1F>(folder, "Wjets_hflav_jpt_C" , "", 100, 0., 500.);    
            book<TH1F>(folder, "Wjets_hflav_jpt_L" , "", 100, 0., 500.);    
            // book<TH1F>(folder, "Wjets_hflavJP_jpt_C" , "", 20, 25., 100.);    
            // book<TH1F>(folder, "Wjets_hflavJP_jpt_L" , "", 20, 25., 100.);    
            // book<TH2F>(folder, "Wjets_hflav_jpt_jeta_C" , "", 20, 25., 100., 16, -2.5, 2.5);    
            // book<TH2F>(folder, "Wjets_hflav_jpt_jeta_L" , "", 20, 25., 100., 16, -2.5, 2.5);    


            // book<TH1F>(folder, "Wja_hflav_CvsL_B" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "Wja_hflav_CvsL_C" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "Wja_hflav_CvsL_L" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "Wjb_hflav_CvsL_B" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "Wjb_hflav_CvsL_C" , "", 55, -1., 1.1);
            // book<TH1F>(folder, "Wjb_hflav_CvsL_L" , "", 55, -1., 1.1);

            book<TH1F>(folder, "Wja_CvsL" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Wjb_CvsL" , "", 55, -1., 1.1);

            book<TH1F>(folder, "Wjets_CvsL" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Wjets_CvsB" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Bjets_CvsL" , "", 55, -1., 1.1);
            book<TH1F>(folder, "Bjets_CvsB" , "", 55, -1., 1.1);
            book<TH1F>(folder, "WjetCSV"   , "", 40, -20., 20.);
            book<TH1F>(folder, "Wjets_CMVA" , "", 55, -1., 1.1);
            //book<TH1F>(folder, "Wjets_DeepCSVCvsL" , "", 55,  0., 1.1);
            //book<TH1F>(folder, "Wjets_DeepCSVCvsB" , "", 55,  0., 1.1);
            book<TH1F>(folder, "Wjets_DeepCSVb" , "", 55,  0., 1.1);
            book<TH1F>(folder, "Wjets_DeepCSVl" , "", 55,  0., 1.1);
            book<TH1F>(folder, "Wjets_DeepCSVbb", "", 55,  0., 1.1);
            book<TH1F>(folder, "Wjets_DeepCSVc" , "", 55,  0., 1.1);
            book<TH1F>(folder, "Wjets_DeepCSVcc", "", 55,  0., 1.1);
            book<TH1F>(folder, "Wjets_DeepCSVbD", "", 55,  0., 1.1);
        }

        void fill_notag_plots(string folder, Permutation &hyp, systematics::SysShifts shift){
            auto dir = histos_.find(folder);
            dir->second["evt_weight"].fill(evt_weight_);
            dir->second["btag_sf"].fill(btag_sf_.scale_factor({hyp.BHad(), hyp.BLep()}, shift), evt_weight_);
            dir->second["muon_sf"].fill( muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta()), evt_weight_ );
            dir->second["muon_sf_vs_mu_pt"].fill( object_selector_.muon()->Pt(), muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta()), evt_weight_ );
            dir->second["muon_sf_vs_mu_eta"].fill( object_selector_.muon()->Eta(), muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta()), evt_weight_ );
            dir->second["njets"    ].fill(object_selector_.clean_jets().size(), evt_weight_);
            dir->second["lep_b_pt" ].fill(hyp.BLep()->Pt(), evt_weight_);
            dir->second["had_b_pt" ].fill(hyp.BHad()->Pt(), evt_weight_);
            dir->second["lep_pt"   ].fill(hyp.L()->Pt(), evt_weight_);
            dir->second["Whad_mass"].fill(hyp.WHad().M(), evt_weight_);
            dir->second["Whad_DR"  ].fill(hyp.WJa()->DeltaR(*hyp.WJb()), evt_weight_);
            dir->second["thad_mass"].fill(hyp.THad().M() , evt_weight_);
            dir->second["thad_pt"].fill(hyp.THad().Pt() , evt_weight_);
            dir->second["WjetCSV"  ].fill(hyp.WJa()->csvIncl(), evt_weight_);
            dir->second["WjetCSV"  ].fill(hyp.WJb()->csvIncl(), evt_weight_);

            const IDJet *lb = (hyp.BHad()->Pt() > hyp.BLep()->Pt()) ? hyp.BHad() : hyp.BLep();
            const IDJet *sb = (hyp.BHad()->Pt() > hyp.BLep()->Pt()) ? hyp.BLep() : hyp.BHad();
            const IDJet *lj = (hyp.WJa()->Pt()  > hyp.WJb()->Pt()) ? hyp.WJa() : hyp.WJb();
            const IDJet *sj = (hyp.WJa()->Pt()  > hyp.WJb()->Pt()) ? hyp.WJb() : hyp.WJa();

            dir->second["Whad_jet_pts" ].fill(lj->Pt(), sj->Pt()  , evt_weight_);		
            dir->second["BJet_jet_pts" ].fill(lb->Pt(), sb->Pt()  , evt_weight_);		
            dir->second["leadB_leadW_pts"].fill(lb->Pt(), lj->Pt(), evt_weight_);
            dir->second["subB_subW_pts"  ].fill(sb->Pt(), sj->Pt(), evt_weight_);

            string wj="Wja";
            for(IDJet* jet : {hyp.WJa(), hyp.WJb()}){
                int hflav = fabs(jet->hadronFlavour());
                int pflav = fabs(jet->partonFlavour());
                string hstr;
                if(hflav == 5) hstr="B";
                else if(hflav == 4) {
                    hstr="C";
                    dir->second["Wjets_hflav_jpt_C"].fill(jet->Pt(), evt_weight_);
                    // dir->second["Wjets_hflavJP_jpt_C"].fill(jet->Pt(), evt_weight_);
                    // dir->second["Wjets_hflav_jpt_jeta_C"].fill(jet->Pt(), jet->Eta(), evt_weight_);
                }
                // else if(pflav == 3) {
                //   hstr="S";
                //   dir->second["Wjets_hflav_jpt_LS"].fill(jet->Pt(), evt_weight_);
                // }
                else {
                    hstr="L"; 
                    dir->second["Wjets_hflav_jpt_L"].fill(jet->Pt(), evt_weight_);
                    // dir->second["Wjets_hflavJP_jpt_L"].fill(jet->Pt(), evt_weight_);
                    // dir->second["Wjets_hflav_jpt_jeta_L"].fill(jet->Pt(), jet->Eta(), evt_weight_);
                }

                dir->second["btag_sf_"+hstr].fill(btag_sf_.scale_factor({hyp.BHad(), hyp.BLep()}, shift), evt_weight_);
                dir->second["btag_sf_vs_pt_"+hstr].fill(jet->Pt(), btag_sf_.scale_factor({hyp.BHad(), hyp.BLep()}, shift), evt_weight_);
                dir->second["Wjets_hflav_CvsL_"+hstr].fill(jet->CvsLtag(), evt_weight_);
                dir->second["Wjets_hflav_CvsB_"+hstr].fill(jet->CvsBtag(), evt_weight_);

                // DeepCSV Discriminant values for each wjet flavor
                //dir->second["Wjets_hflav_DeepCSVCvsL_"+hstr].fill(jet->DeepCSVCvsLtag(), evt_weight_);
                //dir->second["Wjets_hflav_DeepCSVCvsB_"+hstr].fill(jet->DeepCSVCvsBtag(), evt_weight_);
                dir->second["Wjets_hflav_DeepCSVbD_"+hstr].fill(jet->DeepCSVProbB()+jet->DeepCSVProbBB(), evt_weight_);
                // dir->second[wj+"_hflav_CvsL_"+hstr].fill(jet->CvsLtag(), evt_weight_);
            }

            dir->second["Wja_CvsL"].fill(hyp.WJa()->CvsLtag() , evt_weight_);
            dir->second["Wjb_CvsL"].fill(hyp.WJb()->CvsLtag() , evt_weight_);

            dir->second["Wjets_CvsL"].fill(hyp.WJa()->CvsLtag() , evt_weight_);
            dir->second["Wjets_CvsL"].fill(hyp.WJb()->CvsLtag() , evt_weight_);
            dir->second["Wjets_CvsB"].fill(hyp.WJa()->CvsBtag() , evt_weight_);
            dir->second["Wjets_CvsB"].fill(hyp.WJb()->CvsBtag() , evt_weight_);
            dir->second["Bjets_CvsL"].fill(hyp.BHad()->CvsLtag(), evt_weight_);
            dir->second["Bjets_CvsL"].fill(hyp.BLep()->CvsLtag(), evt_weight_);
            dir->second["Bjets_CvsB"].fill(hyp.BHad()->CvsBtag(), evt_weight_);
            dir->second["Bjets_CvsB"].fill(hyp.BLep()->CvsBtag(), evt_weight_);
            dir->second["Wjets_CMVA"].fill(hyp.WJa()->CombinedMVA(), evt_weight_);
            dir->second["Wjets_CMVA"].fill(hyp.WJb()->CombinedMVA(), evt_weight_);

            //dir->second["Wjets_DeepCSVCvsL" ].fill(hyp.WJa()->DeepCSVCvsLtag(), evt_weight_);
            //dir->second["Wjets_DeepCSVCvsB" ].fill(hyp.WJa()->DeepCSVCvsBtag(), evt_weight_);
            dir->second["Wjets_DeepCSVb" ].fill(hyp.WJa()->DeepCSVProbB(), evt_weight_);
            dir->second["Wjets_DeepCSVl" ].fill(hyp.WJa()->DeepCSVProbUDSG(), evt_weight_);
            dir->second["Wjets_DeepCSVbb"].fill(hyp.WJa()->DeepCSVProbBB(), evt_weight_);
            dir->second["Wjets_DeepCSVc" ].fill(hyp.WJa()->DeepCSVProbC(), evt_weight_);
            dir->second["Wjets_DeepCSVcc"].fill(hyp.WJa()->DeepCSVProbCC(), evt_weight_);
            dir->second["Wjets_DeepCSVbD"].fill(hyp.WJa()->DeepCSVProbB()+hyp.WJa()->DeepCSVProbBB(), evt_weight_);

            //dir->second["Wjets_DeepCSVCvsL" ].fill(hyp.WJb()->DeepCSVCvsLtag(), evt_weight_);
            //dir->second["Wjets_DeepCSVCvsB" ].fill(hyp.WJb()->DeepCSVCvsBtag(), evt_weight_);
            dir->second["Wjets_DeepCSVb" ].fill(hyp.WJb()->DeepCSVProbB(), evt_weight_);
            dir->second["Wjets_DeepCSVl" ].fill(hyp.WJb()->DeepCSVProbUDSG(), evt_weight_);
            dir->second["Wjets_DeepCSVbb"].fill(hyp.WJb()->DeepCSVProbBB(), evt_weight_);
            dir->second["Wjets_DeepCSVc" ].fill(hyp.WJb()->DeepCSVProbC(), evt_weight_);
            dir->second["Wjets_DeepCSVcc"].fill(hyp.WJb()->DeepCSVProbCC(), evt_weight_);
            dir->second["Wjets_DeepCSVbD"].fill(hyp.WJb()->DeepCSVProbB()+hyp.WJb()->DeepCSVProbBB(), evt_weight_);
        }

        void book_jet_plots(string folder, bool extended){
            book<TH1F>(folder, "hadronflav","hflav", 28, -0.5, 27.5);
            if(!extended) return;
            book<TH1F>(folder, "eta"	,"eta"	, 100, -5, 5);
            book<TH1F>(folder, "pt" 	,"pt" 	, 100, 0 , 500);
            book<TH1F>(folder, "phi"	,"phi"	, 100, -4, 4);
            book<TH1F>(folder, "pflav_smart","pflav", 55, -27.5, 27.5);
            book<TH1F>(folder, "abs_pflav_smart","pflav", 28, -0.5, 27.5);
            book<TH1F>(folder, "energy", ";E_{jet} (GeV)", 100, 0., 500.);
            book<TH1F>(folder, "ncharged", "", 50, 0., 50.);						
            book<TH1F>(folder, "nneutral", "", 50, 0., 50.);						
            book<TH1F>(folder, "ntotal"  , "", 50, 0., 50.);						
        }

        void fill_jet_plots(string folder, const IDJet* jet, bool extended){//, const Genparticle* genp=0){
            auto dir = histos_.find(folder);
            dir->second["hadronflav"].fill(fabs(jet->hadronFlavour()), evt_weight_);
            if(!extended) return;
            dir->second["eta"	 ].fill(jet->Eta(), evt_weight_);
            dir->second["pt" 	 ].fill(jet->Pt() , evt_weight_);
            dir->second["phi"	 ].fill(jet->Phi(), evt_weight_);
            dir->second["energy"].fill(jet->E(), evt_weight_);

            dir->second["ncharged"].fill(jet->numChargedHadrons(), evt_weight_);
            dir->second["nneutral"].fill(jet->numNeutralHadrons(), evt_weight_);
            dir->second["ntotal"  ].fill(jet->numChargedHadrons()+jet->numNeutralHadrons(), evt_weight_);

            dir->second["pflav_smart"].fill(jet->flavor(), evt_weight_);
            dir->second["abs_pflav_smart"].fill(fabs(jet->flavor()), evt_weight_);
        }

        virtual void begin()
        {
            Logger::log().debug() << "ctag_eff::begin" << endl;
            outFile_.cd();
            vector<string> folders;

            //FIXME use folders as root dir, makes much more sense!
            //Would also be nice to have nothing instead of "all" for the others
            //semilep_visible_right becomes a SubdirectoryView in the plotter
            if(isTTbar_) folders = {"semilep_visible_right", "semilep_right_thad",
                "semilep_right_tlep", "semilep_right_whad", "semilep_wrong", "other"};
            else folders = {""};
            // string folders[] = {"all", "semilep_visible_right", "semilep_right_thad", 
            // 										"semilep_right_tlep", "semilep_right_whad", "semilep_wrong", "other"};
            string wjet_folders[] = {"leading", "subleading"};
            string tagging[] = {"lead_tagged", "sublead_tagged", "both_tagged", "both_untagged"};

            //TH1::AddDirectory(true);
            //if(isTTbar_) book_hyp_plots("gen");
            for(auto& genCategory : folders){			
                for(auto& sys : systematics_){
                    string gcategory;
                    if(!genCategory.empty()) gcategory  = genCategory +"/";
                    string sys_name = systematics::shift_to_name.at(sys);
                    if(genCategory.empty() || genCategory == "semilep_visible_right") {
                        book_presel_plots(gcategory+sys_name+"/preselection");
                        book_combo_plots(gcategory+sys_name+"/permutations");
                    }
                    string criterion = cut_ordering_;

                    for(auto& wp_item : working_points_) {
                        Logger::log().debug()	<< "Booking plots for " << wp_item.first << endl;
                        string working_point = wp_item.first;
                        string base;
                        if(!genCategory.empty()) base  = genCategory +"/";
                        base += sys_name + "/" + criterion + "/" + working_point;
                        book<TH1F>(base, "dummy"  , "", 1, 0., 500.);
                        for(auto& tag : tagging){
                            if(working_point == "notag" && tag != "both_untagged") continue;
                            string folder = base + "/" + tag;

                            book_combo_plots(folder);
                            if(pdfs_ && sys == systematics::SysShifts::NOSYS) {
                                book_pdf_plots(folder);
                            }
                            bool notag = (working_point == "notag");
                            for(auto& w_folder : wjet_folders){
                                string current = folder + "/" + w_folder;
                                book_jet_plots(current, notag);
                            }
                            if(notag && sys == systematics::SysShifts::NOSYS) {
                                book_notag_plots(folder);
                            } //if(working_point == "notag" and sytematic==NOSYS)
                            //if(notag) {
                            //	book_notag_plots(folder);
                            //} //if(working_point == "notag")
                        }//for(auto& tag : tagging)
                    }//for(auto& wp_item : working_points_)
                }//for(auto& sys : systematics){
            }//for(auto& genCategory : folders)
            }

            void fill(string folder, Permutation &hyp, bool extended){//, TTbarHypothesis *genHyp=0) {
                auto dir = histos_.find(folder);
                if(dir == histos_.end()) {
                    Logger::log().error() << "fill: histogram folder: " << folder <<
                        " does not exists!" << endl;
                    throw 40;
                }

                fill_combo_plots(folder, hyp);

                //FILL JET INFORMATION
                const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
                const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();			
                fill_jet_plots(folder + "/leading",       leading, extended);
                fill_jet_plots(folder + "/subleading", subleading, extended);
            }

            /*void fill_other_jet_plots(string folder, Permutation &hyp, std::function<bool(const IDJet*)>& fcn, float weight) {
              vector<IDJet*>& jets = object_selector_.clean_jets();
            //folder += sys_name + "/" + criterion + "/" + working_point;
            auto dir = histos_.find(folder);
            if(dir == histos_.end()) {
            Logger::log().error() << "fill_other_jet_plots: histogram folder: " << folder <<
            " does not exists!" << endl;
            throw 40;
            }
            const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
            const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();
            // dir->second["csvL_csvS" ].fill(leading->csvIncl(), subleading->csvIncl(), weight);

            set<IDJet*> hypjets;
            hypjets.insert(hyp.WJa());
            hypjets.insert(hyp.WJb());
            hypjets.insert(hyp.BLep());
            hypjets.insert(hyp.BHad());
            }*/

            string get_wjet_category(Permutation &hyp, std::function<bool(const IDJet*)>& fcn) {
                const IDJet *leading    = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJa() : hyp.WJb();
                const IDJet *subleading = (hyp.WJa()->E() > hyp.WJb()->E()) ? hyp.WJb() : hyp.WJa();
                bool lead_tag = fcn(leading   );
                bool sub_tag  = fcn(subleading);
                if(lead_tag && sub_tag) return "/both_tagged";
                else if(lead_tag) return "/lead_tagged";
                else if(sub_tag)  return "/sublead_tagged";
                return "/both_untagged";
            }

            TTNaming get_ttdir_name(Permutation &gen, Permutation &reco) {
                if(reco.IsCorrect(gen)) return TTNaming::RIGHT;
                else if(reco.IsTHadCorrect(gen)) return TTNaming::RIGHT_THAD; 
                else if(gen.IsWHadComplete() && reco.IsWHadCorrect(gen)) return TTNaming::RIGHT_WHAD;
                else if(reco.IsBLepCorrect(gen)) return TTNaming::RIGHT_TLEP;
                else if(genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) return TTNaming::WRONG;
                else return TTNaming::OTHER;
            }


            // top pT reweighting formula based on section 5.5 of AN-16-272
            double top_pt_reweighting( string sys_var ){
                double top_pt1, top_pt2;

                if( sys_var == "nominal" ){top_pt1 = 0.; top_pt2 = 0.;}
                else if( sys_var == "up_nom" ){top_pt1 = 1.; top_pt2 = 0.;}
                else if( sys_var == "down_nom" ){top_pt1 = -1.; top_pt2 = 0.;}
                else if( sys_var == "nom_up" ){top_pt1 = 0.; top_pt2 = 1.;}
                else if( sys_var == "nom_down" ){top_pt1 = 0.; top_pt2 = -1.;}
                else{
                    Logger::log().error() << sys_var << " isn't a valid option for the top pt reweighting!" << endl;
                    throw 40;
                }

                double p0 =  6.15025*pow(10., -2.) + 3.243*pow(10., -2.)*top_pt1 - 4.353*pow(10., -7.)*top_pt2;
                double p1 = -5.17833*pow(10., -4.) - 1.404*pow(10., -4.)*top_pt1 - 1.005*pow(10., -4.)*top_pt2;

                double weight = exp(p0+p1*( (genp_selector_.top()->Pt()+genp_selector_.tbar()->Pt())/2 ));
                return weight;
            }

            void process_evt(systematics::SysShifts shift, URStreamer &event){
                tracker_.track("start");
                //float weight = 1.;

                //select reco objects
                if( !object_selector_.select(event, shift) ) return;
                tracker_.track("object selection");

                //MC Weight for lepton selection
                if(!isData_) { 
                    if(object_selector_.tight_muons().size() == 1)
                        evt_weight_ *= muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta());
                    if(object_selector_.tight_electrons().size() == 1)
                        evt_weight_ *= electron_sf_.get_sf(object_selector_.electron()->Pt(), object_selector_.electron()->etaSC());
                }

                //tracker_.track("after lepton sf applied");

                string sys_name = systematics::shift_to_name.at(shift);
                string presel_dir = sys_name;
                if(isTTbar_){
                    presel_dir = naming_.at(TTNaming::RIGHT) + "/" + sys_name;
                }
                fill_presel_plots(presel_dir+"/preselection", event);
                //tracker_.track("after presel plots");

                if( !permutator_.preselection(
                            object_selector_.clean_jets(), object_selector_.lepton(), 
                            object_selector_.met() ) ) return;
                tracker_.track("perm preselection");


                //Find best permutation
                Permutation best_permutation;
                for(auto test_perm : permutator_.pemutations()) {
                    solver_.Solve(test_perm);
                    double bjet_lpt = Max(test_perm.BHad()->Pt(), test_perm.BLep()->Pt());
                    fill_combo_plots(presel_dir+"/permutations", test_perm);
                    //if(bjet_lpt < test_perm.WJa()->Pt()) continue;
                    if(ordering_fcn_(test_perm, best_permutation)){
                        best_permutation = test_perm;
                    }
                }
                if(!best_permutation.IsComplete() || best_permutation.Prob() > 1E9) return; //FIXME, is right??? best_permutation.Prob() > 1E9
                tracker_.track("best perm");

                //NOW, cut on btagging
                if(!(best_permutation.BHad()->BTagId(cut_tight_b_) || best_permutation.BLep()->BTagId(cut_tight_b_))) {return;}
                if(!(best_permutation.BHad()->BTagId(cut_loose_b_) && best_permutation.BLep()->BTagId(cut_loose_b_))) {return;}
                tracker_.track("perm btag");

                size_t capped_jets_size = permutator_.capped_jets().size();
                best_permutation.permutating_jets(capped_jets_size);

                //find mc weight for btag
                if(!isData_){
                    evt_weight_ *= btag_sf_.scale_factor({best_permutation.BHad(), best_permutation.BLep()}, shift);
                }
                //find mc weight for top pt for ttbar events
                if( isTTbar_ ){
                    double top_pt_weight = top_pt_reweighting("nominal");
                    evt_weight_ *= top_pt_weight;
                    //cout << "evt weight: " << evt_weight_ << endl;
                }

                //Gen matching
                Permutation matched_perm;
                // if(isTTbar_) Logger::log().debug() << " --  GEN  --" << genp_selector_.ttbar_final_system() << std::endl;
                if(isTTbar_ && genp_selector_.ttbar_system().type == GenTTBar::DecayType::SEMILEP) {
                    matched_perm = matcher_.match(
                            genp_selector_.ttbar_final_system(),
                            object_selector_.clean_jets(), 
                            object_selector_.veto_electrons(),
                            object_selector_.veto_muons()
                            );
                    matched_perm.SetMET(object_selector_.met());
                }
                tracker_.track("matched perm");

                string sys_dir = sys_name;

                double weight = evt_weight_;
                string ttsubdir = "";
                if(isTTbar_) { 			
                    TTNaming dir_id = get_ttdir_name(matched_perm, best_permutation);
                    ttsubdir = naming_.at(dir_id) + "/";
                }

                for(auto& wpoint : working_points_){        
                    auto wp_sf = wp_SFs_.find(wpoint.first);
                    if(!isData_ && wp_sf != wp_SFs_.end()) evt_weight_ = weight*wp_sf->second->scale_factor({best_permutation.WJa(), best_permutation.WJb()}, shift);
                    else evt_weight_ = weight;

                    string jet_category = get_wjet_category(best_permutation, wpoint.second);
                    string folder = ttsubdir+sys_dir+"/"+cut_ordering_+"/"+wpoint.first;
                    //fill_other_jet_plots(folder, best_permutation, wpoint.second, evt_weight_);
                    folder += jet_category;
                    //Logger::log().debug() << "filling: " << folder << endl;
                    bool notag = (wpoint.first == "notag");
                    if(notag && shift == systematics::SysShifts::NOSYS) fill_notag_plots(folder, best_permutation, shift); //fill notag plots only for NOSYS
                    //if(notag) fill_notag_plots(folder, best_permutation);
                    fill(folder, best_permutation, notag);
                    if(pdfs_ && shift == systematics::SysShifts::NOSYS) fill_pdf_plots(folder, best_permutation, event);
                }
            }

            //This method is called once every file, contains the event loop
            //run your proper analysis here
            virtual void analyze()
            {
                opts::variables_map &values = URParser::instance().values();
                int limit = values["limit"].as<int>();
                int skip  = values["skip"].as<int>();
                int report = values["report"].as<int>();

                if(evt_idx_ >= limit) return;
                Logger::log().debug() << "ctag_eff::analyze" << endl;
                URStreamer event(tree_);

                string output_file = values["output"].as<std::string>();
                string sample = systematics::get_sample(output_file);
                //if( boost::starts_with(sample, "data") && sample.find("E") != std::string::npos ) Logger::log().debug() << endl << endl << "You're running on " << sample << ". MAKE SURE YOU ARE USING THE RIGHT RUN SELECTION!!" << endl << endl;

                tracker_.deactivate();
                while(event.next())
                {

                    //if( boost::starts_with(sample, "data") && sample.find("E") != std::string::npos ){
                    //    //cout << "Run number: " << event.run << endl;
                    //    //if( event.run > 304671 ) continue; // for combining C to E data runs
                    //    //if( event.run <= 304671 ) continue; // for combining E to F data runs
                    //}
                    // if(evt_idx_ % 1000 == 0) Logger::log().warning() << "Beginning event: " <<
                    //                           evt_idx_ << endl;
                    if(limit > 0 && evt_idx_ > limit) {
                        return;
                    }
                    evt_idx_++;
                    if(skip > 0 && evt_idx_ < skip) {
                        continue;
                    }
                    if(evt_idx_ % report == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

                    //long and time consuming
                    if(isTTbar_){
                        genp_selector_.select(event);			
                    }

                    for(auto shift : systematics_){

                        evt_weight_ = (isData_) ? 1. : mc_weights_.evt_weight(event, shift);

                        //cout << "Ev weight: " << evt_weight_ << endl;
                        if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
                        //Logger::log().debug() << "processing: " << shift << endl;
                        process_evt(shift, event);
                        if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
                    }

                } //while(event.next())

                Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
            }

            //this method is called at the end of the job, by default saves
            //every histogram/tree produced, override it if you need something more
            virtual void end(){
                outFile_.Write();
                tracker_.writeTo(outFile_);
            }

            //do you need command-line or cfg options? If so implement this 
            //method to book the options you need. CLI parsing is provided
            //by AnalysisFW/interface/URParser.h and uses boost::program_options
            //look here for a quickstart tutorial: 
            //http://www.boost.org/doc/libs/1_51_0/doc/html/program_options/tutorial.html
            static void setOptions() {
                URParser &parser = URParser::instance();
                opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
                opts.add_options()
                    ("nosys", opts::value<int>()->default_value(0), "do not run systematics")
                    ("nopdf", opts::value<int>()->default_value(0), "do not run pdf uncertainties")
                    ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                    ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                    ("report", opts::value<int>()->default_value(10000), "report every");

                parser.addCfgParameter<std::string>("general", "ctag_sffile", "");
                //parser.addCfgParameter<std::string>("general", "deepctag_sffile", "");
                //parser.addCfgParameter<std::string>("general", "cmva_sffile", "");
                parser.addCfgParameter<std::string>("general", "csv_sffile", "");
                parser.addCfgParameter<std::string>("general", "deepcsv_sffile", "");
                parser.addCfgParameter<std::string>("general", "wjets_efficiencies", "");

                parser.addCfgParameter<std::string>("permutations", "ordering", "ID to be applied");
                //("permutations.ordering", opts::value<string>()->default_value("mass_discriminant"));
            }
            };

            //make it executable
            int main(int argc, char *argv[])
            {
                URParser &parser = URParser::instance(argc, argv);
                URDriver<ctag_eff> test;
                return test.run();
            }
