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
//#include "Analyses/URTTbar/interface/JetScaler.h"
#include "TUUID.h"   
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"

#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "Analyses/URTTbar/interface/DR_TTGenMatcher.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"
#include <sstream>
#include <unordered_map>
#include "URAnalysis/AnalysisFW/interface/EventList.h"
#include "Analyses/URTTbar/interface/Hypotheses.h"
#include "TROOT.h"
//#include <fstream>
//#include "JetMETCorrections/Modules/interface/JetResolution.h"
//#include "Analyses/URTTbar/interface/JERFile.h"

using namespace TMath;
using namespace systematics;
typedef SysShifts Sys;

class htt_simple : public AnalyzerBase
{
    public:
        //enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_WHAD, RIGHT_TLEP, WRONG, OTHER};
        struct SyncInfo {
            int Run;
            long int LumiSection;
            long int Event; 
            bool RecoSuccess; 
            float MassTT;
            bool hasMuon;
        };

    private:
        //histograms and helpers
        unordered_map<string, unordered_map<string, RObject> > histos_;
        CutFlowTracker tracker_;

        //switches
        bool isTTbar_, isSignal_, isData_, runsys_, has_pdfs_, optim_;

        //
        bool isTTJetsM1000_, isTTJetsM700_;
        //

        //selectors and helpers
        TTGenParticleSelector genp_selector_;
        TTGenMatcher matcher_;
        DR_TTGenMatcher dr_matcher_;
        TTPermutator permutator_;
        TTObjectSelector object_selector_;
        float evt_weight_;
        TRandom3 randomizer_;// = TRandom3(98765);
        MCWeightProducer mc_weights_;
        BTagSFProducer btag_sf_;
        TTBarSolver solver_;
        PDFuncertainty pdf_uncs_;

        //Scale factors
        LeptonSF electron_sf_, muon_sf_;

        //    JME::JetResolution pt_jet_res_; // initialize jet energy pt resolution
        //    JME::JetResolutionScaleFactor jet_res_sf_; // initialize jet energy resolution scale factors


        unsigned long evt_idx_ = 0;
        vector<systematics::SysShifts> systematics_;

        //cuts
        IDJet::BTag cut_tight_b_=IDJet::BTag::NONE, cut_loose_b_=IDJet::BTag::NONE;

        //sync
        bool sync_;
        TTree *sync_tree_;
        SyncInfo sync_info_;

    public:
        inline double MT(TLorentzVector *l, TLorentzVector *met) {
            return sqrt(pow(l->Pt() + met->Pt(), 2) - pow(l->Px() + met->Px(), 2) - pow(l->Py() + met->Py(), 2));
        }

        htt_simple(const std::string output_filename):
            AnalyzerBase("htt_simple", output_filename), 
            tracker_(),
            genp_selector_(TTGenParticleSelector::SelMode::LHE),
            //genp_selector_(),
            matcher_(),
            permutator_(),
            object_selector_(),
            mc_weights_(),
            evt_weight_(1.),
            electron_sf_("electron_sf", false),
            muon_sf_("muon_sf"),
            randomizer_(),    
            btag_sf_("permutations.tightb", "permutations.looseb"),
            solver_(true),
            pdf_uncs_(108),
            systematics_(),
            sync_(false),
            sync_tree_(0),
            sync_info_(){
                object_selector_.allow_loose();//allow loose but not tight lepton events
                Logger::log().debug() << "htt_simple ctor" << endl;
                cut_tight_b_ = btag_sf_.tight_cut();
                cut_loose_b_ = btag_sf_.loose_cut();
                //cout << "bcuts " << cut_tight_b_ << " " << cut_loose_b_ << endl;
                //  " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;

                URParser &parser = URParser::instance();

                //    parser.addCfgParameter<string>("JERC", "JER_SF","");
                //    parser.addCfgParameter<string>("JERC", "PT_JER","");
                //    parser.parseArguments();
                //
                //    DataFile jer_sf_fname_(parser.getCfgPar<string>("JERC","JER_SF"));
                //    DataFile pt_jer_fname_(parser.getCfgPar<string>("JERC","PT_JER"));
                //
                //    jet_res_sf_ = JME::JetResolutionScaleFactor(jer_sf_fname_.path());
                //    pt_jet_res_ = JME::JetResolution(pt_jer_fname_.path());

                TUUID id;  
                randomizer_.SetSeed(id.Hash());   

                //find out which sample are we running on
                opts::variables_map &values = parser.values();
                string output_file = values["output"].as<std::string>();
                string sample = systematics::get_sample(output_file);
                isSignal_ = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
                isTTbar_ = boost::starts_with(sample, "ttJets");
                isData_  = boost::starts_with(sample, "data");

                //
                isTTJetsM1000_ = boost::starts_with(sample, "ttJetsM1000");
                isTTJetsM700_ = boost::starts_with(sample, "ttJetsM700");
                //
                //set tracker
                if(!values.count("noweights")) tracker_.use_weight(&evt_weight_);
                object_selector_.set_tracker(&tracker_);
                permutator_.set_tracker(&tracker_);

                if(isData_){
                    if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
                    else object_selector_.lepton_type(1);
                }
                sync_ = values.count("sync");
                optim_ = values.count("optimization");
                has_pdfs_ = !(
                        boost::starts_with(sample, "QCD") || 
                        boost::starts_with(sample, "singletbar_tW") || boost::starts_with(sample, "singlet_tW") ||
                        boost::starts_with(sample, "WZ") || boost::starts_with(sample, "WW") || boost::starts_with(sample, "ZZ")
                        );
                bool isTTShift = boost::starts_with(sample, "ttJets_");
                if(sync_) tracker_.verbose(true);
                if(!isData_) mc_weights_.init(sample, has_pdfs_);

                Logger::log().debug() << "  NOSYS: " << values.count("nosys") << endl;

                runsys_ = !(values.count("nosys") || isData_ || isTTShift || sync_);
                if(!runsys_)
                    systematics_ = {systematics::SysShifts::NOSYS};
                else {
                    systematics_ = {
                        Sys::NOSYS, 
                        Sys::JES_UP,  Sys::JES_DW, 
                        Sys::JER_UP,  Sys::JER_DW, 
                        Sys::MET_UP,  Sys::MET_DW,
                        Sys::PU_UP,   Sys::PU_DW,
                        Sys::BEFF_UP, Sys::BEFF_DW, 
                        Sys::BFAKE_UP, Sys::BFAKE_DW, 
                        Sys::LEPEFF_UP, Sys::LEPEFF_DW,
                    };
                    if(isTTbar_) {
                        systematics_.push_back(Sys::HDAMP_UP);
                        systematics_.push_back(Sys::HDAMP_DW);
                        systematics_.push_back(Sys::RENORM_UP);	
                        systematics_.push_back(Sys::RENORM_DW); 
                        systematics_.push_back(Sys::FACTOR_UP);	
                        systematics_.push_back(Sys::FACTOR_DW);
                        systematics_.push_back(Sys::RENFACTOR_UP);
                        systematics_.push_back(Sys::RENFACTOR_DW);
                    }
                }
            }

        TDirectory* getDir(std::string path){
            TDirectory* out = (TDirectory*) outFile_.Get(path.c_str());
            if(out) return out;
            outFile_.mkdir(path.c_str());
            return (TDirectory*) outFile_.Get(path.c_str());
        }

        template <class H, typename ... Args>
            void book(std::string folder, std::string name, Args ... args)
            {
                getDir(folder)->cd();
                histos_[folder][name] = RObject::book<H>(name.c_str(), args ...);
            }


    ////////////////////// book and fill 3 jet plots
            // reco
        void book_3J_reco_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", 180, 200., 2000.);
            book<TH1F>(folder+"/Mass", "THad", "", 500, 0., 500.);
            book<TH1F>(folder+"/Costh", "THad", "", 200, -1., 1.);
            book<TH1F>(folder+"/Costh", "TLep", "", 200, -1., 1.);
            book<TH1F>(folder+"/Pt", "THad", "", 100, 0., 1000.);
            book<TH1F>(folder+"/Pt", "TLep", "", 100, 0., 1000.);
            book<TH1F>(folder+"/Pt", "TTbar", "",100, 0., 1000.);
            book<TH1F>(folder+"/Eta", "THad", "", 100, -2.5, 2.5);
            book<TH1F>(folder+"/Eta", "TLep", "", 100, -2.5, 2.5);
        }
        void fill_3J_pre_alpha_reco_plots( string folder, Permutation &perm ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(perm.LVect().M(), evt_weight_);
            pt_dir->second["TTbar"].fill(perm.LVect().Pt(), evt_weight_);

            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(perm.THad().M(), evt_weight_);
            pt_dir->second["THad"].fill(perm.THad().Pt(), evt_weight_);
            eta_dir->second["THad"].fill(perm.THad().Eta(), evt_weight_);
            costh_dir->second["THad"].fill(reco_cosths.first, evt_weight_);
            pt_dir->second["TLep"].fill(perm.TLep().Pt(), evt_weight_);
            eta_dir->second["TLep"].fill(perm.TLep().Eta(), evt_weight_);
            costh_dir->second["TLep"].fill(reco_cosths.second, evt_weight_);
        }
        void fill_3J_post_alpha_reco_plots( string folder, Permutation &perm ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

                // costh
            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

                // corrected thad
            std::pair< TLorentzVector, double > corrected_thad = alpha_thad(perm); // corrected thad tvec, costh* value

            mass_dir->second["THad"].fill( corrected_thad.first.M(), evt_weight_ );
            mass_dir->second["TTbar"].fill( (corrected_thad.first + perm.TLep()).M(), evt_weight_ );

            costh_dir->second["THad"].fill( corrected_thad.second, evt_weight_ );
            costh_dir->second["TLep"].fill(reco_cosths.second, evt_weight_);

            pt_dir->second["THad"].fill( corrected_thad.first.Pt(), evt_weight_ );
            pt_dir->second["TLep"].fill(perm.TLep().Pt(), evt_weight_);
            pt_dir->second["TTbar"].fill( (corrected_thad.first + perm.TLep()).Pt(), evt_weight_ );

            eta_dir->second["THad"].fill( corrected_thad.first.Eta(), evt_weight_ );
            eta_dir->second["TLep"].fill(perm.TLep().Eta(), evt_weight_);
        }

            // resolution
        void book_3J_reso_plots( string folder ){
            book<TH1F>(folder+"/Mass", "TTbar", "", 500, -1000., 1000.);
            book<TH1F>(folder+"/Mass", "THad", "", 500, -1000., 500.);
            book<TH1F>(folder+"/Costh", "THad", "", 200, -2., 2.);
            book<TH1F>(folder+"/Costh", "TLep", "", 200, -2., 2.);
            book<TH1F>(folder+"/Pt", "THad", "", 200, -500., 500.);
            book<TH1F>(folder+"/Pt", "TLep", "", 200, -500., 500.);
            book<TH1F>(folder+"/Pt", "TTbar", "",200, -500., 500.);
            book<TH1F>(folder+"/Eta", "THad", "", 200, -5., 5.);
            book<TH1F>(folder+"/Eta", "TLep", "", 200, -5., 5.);
        }
        void fill_3J_pre_alpha_reso_plots( string folder, Permutation &perm, GenTTBar &ttbar ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

            mass_dir->second["TTbar"].fill(ttbar.M() - perm.LVect().M(), evt_weight_);
            pt_dir->second["TTbar"].fill(ttbar.Pt() - perm.LVect().Pt(), evt_weight_);

            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

            mass_dir->second["THad"].fill(ttbar.had_top()->M() - perm.THad().M(), evt_weight_);
            pt_dir->second["THad"].fill(ttbar.had_top()->Pt() - perm.THad().Pt(), evt_weight_);
            eta_dir->second["THad"].fill(ttbar.had_top()->Eta() - perm.THad().Eta(), evt_weight_);
            costh_dir->second["THad"].fill(gen_cosths.first - reco_cosths.first, evt_weight_);
            pt_dir->second["TLep"].fill(ttbar.lep_top()->Pt() - perm.TLep().Pt(), evt_weight_);
            eta_dir->second["TLep"].fill(ttbar.lep_top()->Eta() - perm.TLep().Eta(), evt_weight_);
            costh_dir->second["TLep"].fill(gen_cosths.second - reco_cosths.second, evt_weight_);
        }
        void fill_3J_post_alpha_reso_plots( string folder, Permutation &perm, GenTTBar &ttbar ){
            auto mass_dir = histos_.find(folder+"/Mass");
            auto costh_dir = histos_.find(folder+"/Costh");
            auto pt_dir = histos_.find(folder+"/Pt");
            auto eta_dir = histos_.find(folder+"/Eta");

                // costh
            std::pair< double, double > gen_cosths = gen_costh_tops(ttbar); // < gen thad, tlep costh >
            std::pair< double, double > reco_cosths = reco_costh_tops(perm); // < reco thad, tlep costh >

                // corrected thad
            std::pair< TLorentzVector, double > corrected_thad = alpha_thad(perm); // corrected thad tvec, costh* value


            mass_dir->second["TTbar"].fill( ttbar.M() - (corrected_thad.first + perm.TLep()).M(), evt_weight_ );
            pt_dir->second["TTbar"].fill( ttbar.Pt() - (corrected_thad.first + perm.TLep()).Pt(), evt_weight_ );

            mass_dir->second["THad"].fill( ttbar.had_top()->M() - corrected_thad.first.M(), evt_weight_ );
            pt_dir->second["THad"].fill( ttbar.had_top()->Pt() - corrected_thad.first.Pt(), evt_weight_ );
            eta_dir->second["THad"].fill( ttbar.had_top()->Eta() - corrected_thad.first.Eta(), evt_weight_ );
            costh_dir->second["THad"].fill( gen_cosths.first - corrected_thad.second, evt_weight_ );
            pt_dir->second["TLep"].fill( ttbar.lep_top()->Pt() - perm.TLep().Pt(), evt_weight_ );
            eta_dir->second["TLep"].fill( ttbar.lep_top()->Eta() - perm.TLep().Eta(), evt_weight_ );
            costh_dir->second["TLep"].fill( gen_cosths.second - reco_cosths.second, evt_weight_ );
        }


        //// selection plots

        void fill_3J_combo_plots(string folder, Permutation &hyp){
            if(optim_) return;
            auto dir = histos_.find(folder);
            dir->second["mass_discriminant"].fill(hyp.MassDiscr(), evt_weight_);
            dir->second["nu_discriminant"  ].fill(hyp.NuDiscr(), evt_weight_);
            dir->second["full_discriminant"].fill(hyp.Prob(), evt_weight_);

                // corrected thad
            std::pair< TLorentzVector, double > corrected_thad = alpha_thad(hyp); // corrected thad tvec, costh* value

            //double whad_mass = hyp.WHad().M();
            double thad_mass = corrected_thad.first.M();
            //dir->second["Wmasshad"].fill(whad_mass, evt_weight_);
            dir->second["tmasshad"].fill(thad_mass, evt_weight_);
            dir->second["lbratio" ].fill(hyp.L()->Pt()/hyp.BLep()->Pt(), evt_weight_);
            //double minpt = ((hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJb() : hyp.WJa())->Pt();
            double minpt = hyp.WJa()->Pt();
            dir->second["j2bratio"].fill(minpt/hyp.BHad()->Pt(), evt_weight_);
            dir->second["bjets_pt"].fill(hyp.BHad()->Pt(), evt_weight_);
            dir->second["bjets_pt"].fill(hyp.BLep()->Pt(), evt_weight_);
            dir->second["wjets_pt"].fill(hyp.WJa()->Pt(), evt_weight_);		
            //dir->second["wjets_pt"].fill(hyp.WJb()->Pt(), evt_weight_);		
        }


        void fill_3J_selection_plots(string folder, URStreamer &event, Permutation &hyp, bool pdf){
            fill_presel_plots(folder, event);
            fill_3J_combo_plots(folder, hyp);
            auto dir = histos_.find(folder);


            /// get tlep costh, corrected thad tvec and costh
                // costh
            std::pair< double, double > reco_cosths = reco_costh_tops(hyp); // < reco thad, tlep costh >

                // corrected thad
            std::pair< TLorentzVector, double > corrected_thad = alpha_thad(hyp); // corrected thad tvec, costh* value
            ///

            ////PDF uncertainties
            //if(has_pdfs_) {
            //	const vector<Mcweight>& ws =  event.MCWeights();
            //	for(size_t h = 0 ; h < ws.size() ; ++h) {
            //		dir->second["PDFYields"].fill(h, evt_weight_*ws[h].weights()/ws[0].weights());
            //	}
            //}

            dir->second["pt_thad" ].fill( corrected_thad.first.Pt() , evt_weight_ );
            dir->second["eta_thad"].fill( corrected_thad.first.Eta(), evt_weight_ );
            dir->second["pt_tlep" ].fill( hyp.TLep().Pt() , evt_weight_ );
            dir->second["eta_tlep"].fill( hyp.TLep().Eta(), evt_weight_ );
            dir->second["pt_tt"   ].fill( (corrected_thad.first + hyp.TLep()).Pt() , evt_weight_ );
            dir->second["eta_tt"  ].fill( (corrected_thad.first + hyp.TLep()).Eta(), evt_weight_ );
            dir->second["m_tt"].fill( (corrected_thad.first + hyp.TLep()).M(), evt_weight_ );
            if(pdf) pdf_uncs_.fill_replicas(folder, "m_tt", (corrected_thad.first + hyp.TLep()).M(), evt_weight_, event);

            if(object_selector_.clean_jets().size() == 3)
                dir->second["full_discriminant_3j" ].fill(hyp.Prob(), evt_weight_);
            else if(object_selector_.clean_jets().size() == 4)
                dir->second["full_discriminant_4j" ].fill(hyp.Prob(), evt_weight_);
            else if(object_selector_.clean_jets().size() == 5)
                dir->second["full_discriminant_5j" ].fill(hyp.Prob(), evt_weight_);
            else
                dir->second["full_discriminant_6Pj"].fill(hyp.Prob(), evt_weight_);

            ////choose up/down
            //TLorentzVector whad(hyp.WHad());
            //TLorentzVector j1(*hyp.WJa());
            //TLorentzVector j2(*hyp.WJb());
            //j1.Boost(whad.BoostVector()*-1);
            //j2.Boost(whad.BoostVector()*-1);
            //if(j2.E() > j1.E()) hyp.SwapWJets();

            ////Angular variables
            //hyp::TTbar ttang(hyp);
            //auto ttcm = ttang.to_CM(); 
            //auto tlepcm = ttang.tlep().to_CM();
            //auto thadcm = ttang.thad().to_CM();

            ////top angles
            double tlep_ctstar = min(fabs(reco_cosths.second), 0.99999);
            double thad_ctstar = min(fabs(corrected_thad.second), 0.99999);
            dir->second["tlep_ctstar"].fill(tlep_ctstar, evt_weight_);
            dir->second["thad_ctstar"].fill(thad_ctstar, evt_weight_);
            dir->second["mtt_tlep_ctstar"].fill( (corrected_thad.first + hyp.TLep()).M(), tlep_ctstar, evt_weight_ );

            if(pdf){
                    pdf_uncs_.fill_replicas2D(
                    folder, "mtt_tlep_ctstar", 
                    (corrected_thad.first + hyp.TLep()).M(), tlep_ctstar,
                    //hyp.LVect().M(), fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())),
                    evt_weight_, event
                    );
            }

            //dir->second["cdelta_ld"].fill(tlepcm.W().l().unit3D().Dot(thadcm.W().down().unit3D()), evt_weight_);

            ////helicity frame
            //dir->second["hframe_ctheta_d"].fill(
            //        thadcm.W().down().unit3D().Dot(ttcm.thad().unit3D()), 
            //        evt_weight_
            //        );
        } // fill_3J_selection_plots
    ////////////////////// book and fill 3 jet plots





        void book_combo_plots(string folder){
            if(optim_) return;
            book<TH1F>(folder, "mass_discriminant", "", 20,   0., 20.);
            book<TH1F>(folder, "nu_discriminant", "", 20,   0., 20.);
            book<TH1F>(folder, "full_discriminant", "", 40,   0., 40.);

            book<TH1F>(folder, "tmasshad", "", 100, 0., 500);			
            book<TH1F>(folder, "Wmasshad", "", 100, 0., 500);			
            book<TH1F>(folder, "lbratio" , "", 100, 0., 10.);			
            book<TH1F>(folder, "j2bratio", "", 100, 0., 10.);			
            book<TH1F>(folder, "bjets_pt", "", 500, 0., 500.);			
            book<TH1F>(folder, "wjets_pt", "", 500, 0., 500.);			
        }

        void fill_combo_plots(string folder, const Permutation &hyp){
            if(optim_) return;
            auto dir = histos_.find(folder);
            dir->second["mass_discriminant"].fill(hyp.MassDiscr(), evt_weight_);
            dir->second["nu_discriminant"  ].fill(hyp.NuDiscr(), evt_weight_);
            dir->second["full_discriminant"].fill(hyp.Prob(), evt_weight_);

            double whad_mass = hyp.WHad().M();
            double thad_mass = hyp.THad().M();
            dir->second["Wmasshad"].fill(whad_mass, evt_weight_);
            dir->second["tmasshad"].fill(thad_mass, evt_weight_);
            dir->second["lbratio" ].fill(hyp.L()->Pt()/hyp.BLep()->Pt(), evt_weight_);
            double minpt = ((hyp.WJa()->Pt() > hyp.WJb()->Pt()) ? hyp.WJb() : hyp.WJa())->Pt();
            dir->second["j2bratio"].fill(minpt/hyp.BHad()->Pt(), evt_weight_);
            dir->second["bjets_pt"].fill(hyp.BHad()->Pt(), evt_weight_);
            dir->second["bjets_pt"].fill(hyp.BLep()->Pt(), evt_weight_);
            dir->second["wjets_pt"].fill(hyp.WJa()->Pt(), evt_weight_);		
            dir->second["wjets_pt"].fill(hyp.WJb()->Pt(), evt_weight_);		
        }

        void book_presel_plots(string folder) {
            if(optim_) return;
            book<TH1F>(folder, "nvtx_noweight", "", 100, 0., 100.);
            book<TH1F>(folder, "nvtx", "", 100, 0., 100.);
            book<TH1F>(folder, "rho", "", 100, 0., 100.);
            book<TH1F>(folder, "weight", "", 100, 0., 5.);
            book<TH1F>(folder, "lep_pt"   , ";p_{T}(l) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "lep_eta"  , ";#eta(l) (GeV)", 300, -3, 3);
            book<TH1F>(folder, "jets_pt"  , ";p_{T}(j) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "jets_eta" , ";#eta(j) (GeV)",  300, -3, 3);
            book<TH1F>(folder, "lead_jet_pt"  , ";p_{T}(j) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "lead_jet_eta" , ";#eta(j) (GeV)",  300, -3, 3);
            book<TH1F>(folder, "jets_CSV" , ";#eta(j) (GeV)",  200, -1, 1);
            book<TH1F>(folder, "max_jets_CSV" , ";#eta(j) (GeV)",  200, -1, 1);
            book<TH1F>(folder, "lep_iso", ";#eta(j) (GeV)",  200, 0, 10);
            book<TH1F>(folder, "lep_wp" , ";#eta(j) (GeV)",  4, 0, 4);
            book<TH1F>(folder, "MT" , ";#eta(j) (GeV)",  500, 0, 500);
            book<TH1F>(folder, "MET" , ";#eta(j) (GeV)",  500, 0, 1000);
            book<TH1F>(folder, "METPhi" , ";#eta(j) (GeV)",  314, -1*Pi(), Pi());

            book<TH1F>(folder, "csv" , "",  100, -0.1, 1.1);
            book<TH1F>(folder, "csv_p11" , "",  100, -0.1, 1.1);

            book<TH1F>(folder, "njets"    , "", 50, 0., 50.);

            book<TH2F>(folder, "MT_iso" , ";#eta(j) (GeV)"  ,  10, 0, 100, 10, 0, 1);
            book<TH2F>(folder, "MT_btag" , ";#eta(j) (GeV)" ,  10, 0, 100, 10, 0, 1);
            book<TH2F>(folder, "iso_btag" , ";#eta(j) (GeV)",  10, 0, 1,   10, 0, 1);
            book<TH2D>(folder, "jets_csv_WP", "", 4, 0., 4., 4, 0., 4.);
        }

        int btag_idval(const IDJet* jet) {
            if(jet->BTagId(IDJet::BTag::CSVTIGHT)) return 3;
            else if(jet->BTagId(IDJet::BTag::CSVMEDIUM)) return 2;
            else if(jet->BTagId(IDJet::BTag::CSVLOOSE) ) return 1;
            return 0;
        }

        void fill_presel_plots(string folder, URStreamer &event){
            if(optim_) return;
            auto dir = histos_.find(folder);
            if(dir == histos_.end()) {
                Logger::log().error() << "could not find: " << folder << endl;
                throw 42;
            }

            dir->second["nvtx"].fill(event.vertexs().size(), evt_weight_);
            dir->second["nvtx_noweight"].fill(event.vertexs().size());
            dir->second["rho"].fill(event.rho().value(), evt_weight_);
            dir->second["weight"].fill(evt_weight_);
            dir->second["lep_pt"].fill(object_selector_.lepton()->Pt(), evt_weight_);
            dir->second["lep_eta"].fill(object_selector_.lepton()->Eta(), evt_weight_);
            dir->second["njets" ].fill(object_selector_.clean_jets().size(), evt_weight_);
            double max_csv = -10000;
            double max_pt = -1;
            double max_eta = 0;
            for(IDJet* jet : object_selector_.clean_jets()){
                dir->second["jets_pt"].fill(jet->Pt(), evt_weight_);
                dir->second["jets_eta"].fill(jet->Eta(), evt_weight_);
                dir->second["jets_CSV"].fill(jet->csvIncl(), evt_weight_);
                if(jet->csvIncl() > max_csv) max_csv = jet->csvIncl();
                dir->second["csv"    ].fill(jet->csvIncl(), evt_weight_);
                dir->second["csv_p11"].fill(pow(jet->csvIncl(), 11), evt_weight_);
                if(jet->Pt() > max_pt) {
                    max_pt = jet->Pt();
                    max_eta = jet->Eta();
                }
            }
            dir->second["lead_jet_pt" ].fill(max_pt , evt_weight_);
            dir->second["lead_jet_eta"].fill(max_eta, evt_weight_);

            dir->second["MET"   ].fill(object_selector_.met()->Et() , evt_weight_);
            dir->second["METPhi"].fill(object_selector_.met()->Phi(), evt_weight_);
            double mt = MT(object_selector_.lepton(), object_selector_.met());
            dir->second["MT"].fill(mt, evt_weight_);
            dir->second["max_jets_CSV"].fill(max_csv, evt_weight_);

            auto &clean_jets = object_selector_.clean_jets();
            sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->csvIncl() > B->csvIncl());});		
            dir->second["jets_csv_WP"].fill(btag_idval(clean_jets[0]), btag_idval(clean_jets[1]), evt_weight_);

            double iso = -1.;
            if(object_selector_.lepton_type() == 1) {
                iso = object_selector_.muon()->RelPFIsoDb();
            }
            else {
                auto electron = object_selector_.electron();			
                iso = electron->PFIsolationRho2015();
                int elwp = -1;
                if(electron->TightID25ns()) elwp = 3;
                else if(electron->MediumID25ns()) elwp =2;
                else if(electron->LooseID25ns()) elwp =1;
                else if(electron->VetoID25ns()) elwp =0;
                dir->second["lep_wp" ].fill(elwp, evt_weight_);
            }
            dir->second["lep_iso"].fill(iso, evt_weight_);
            dir->second["MT_iso"  ].fill(mt, iso, evt_weight_);
            dir->second["MT_btag" ].fill(mt, max_csv, evt_weight_);
            dir->second["iso_btag"].fill(iso, max_csv, evt_weight_);
        }

        void book_selection_plots(string folder, bool pdf){
            book_presel_plots(folder);
            book_combo_plots(folder);
            const vector<double> mbinning = {
                250.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 
                480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 610.0, 640.0, 
                680.0, 730.0, 800.0, 920.0, 1200.0};

            //book<TH1F>(folder, "m_tt", "", 180, 200., 2000);			
            book<TH1F>(folder, "m_tt", "", 19, &mbinning[0]);			

            book<TH1F>(folder, "pt_thad" , "", 200, 0., 2000);			
            book<TH1F>(folder, "eta_thad", "", 200, -10., 10);			
            book<TH1F>(folder, "pt_tlep" , "", 200, 0., 2000);			
            book<TH1F>(folder, "eta_tlep", "", 200, -10., 10);			
            book<TH1F>(folder, "pt_tt" , "", 200, 0., 2000);			
            book<TH1F>(folder, "eta_tt", "", 200, -10., 10);			

            book<TH1F>(folder, "full_discriminant_3j", "" , 640, 0., 80.);
            book<TH1F>(folder, "full_discriminant_4j", "" , 640, 0., 80.);
            book<TH1F>(folder, "full_discriminant_5j", "" , 640, 0., 80.);
            book<TH1F>(folder, "full_discriminant_6Pj", "", 640, 0., 80.);

            //Angular variables
            //top angles
            book<TH1F>(folder, "tlep_ctstar", "", 200, -1.0001, 1.0001);
            book<TH1F>(folder, "thad_ctstar", "", 200, -1.0001, 1.0001);

            //delta
            book<TH1F>(folder, "cdelta_ld", "", 200, -1.0001, 1.0001);

            //helicity frame
            book<TH1F>(folder, "hframe_ctheta_d", "", 200, -1.0001, 1.0001);

            //2D plots
            book<TH2F>(folder, "mtt_tlep_ctstar", "", 19, &mbinning[0], 20, 0., 1.0001);
            //book<TH2F>(folder, "mtt_tlep_ctstar", "", 19, &mbinning[0], 5, 0., 1.0001);

            //PDF uncertainties
            //book<TH1D>(folder, "PDFYields", "", 249, 0, 249);		
            if(pdf) {
                auto dir = histos_.find(folder);
                if(dir == histos_.end()) {
                    throw std::runtime_error("Could not access folder");
                }
                pdf_uncs_.book_replicas(folder, "m_tt", dir->second["m_tt"]);
                pdf_uncs_.book_replicas(folder, "mtt_tlep_ctstar", dir->second["mtt_tlep_ctstar"]);
            }
        }

        void fill_selection_plots(string folder, URStreamer &event, Permutation &hyp, bool pdf){
            fill_presel_plots(folder, event);
            fill_combo_plots(folder, hyp);
            auto dir = histos_.find(folder);

            //PDF uncertainties
            //if(has_pdfs_) {
            //	const vector<Mcweight>& ws =  event.MCWeights();
            //	for(size_t h = 0 ; h < ws.size() ; ++h) {
            //		dir->second["PDFYields"].fill(h, evt_weight_*ws[h].weights()/ws[0].weights());
            //	}
            //}
//cout << "fill selection folder: " << folder << endl;

            dir->second["pt_thad" ].fill(hyp.THad().Pt() , evt_weight_);
            dir->second["eta_thad"].fill(hyp.THad().Eta(), evt_weight_);
            dir->second["pt_tlep" ].fill(hyp.TLep().Pt() , evt_weight_);
            dir->second["eta_tlep"].fill(hyp.TLep().Eta(), evt_weight_);
            dir->second["pt_tt"   ].fill(hyp.LVect().Pt() , evt_weight_);
            dir->second["eta_tt"  ].fill(hyp.LVect().Eta(), evt_weight_);
            dir->second["m_tt"].fill(hyp.LVect().M(), evt_weight_);
            if(pdf) pdf_uncs_.fill_replicas(folder, "m_tt", hyp.LVect().M(), evt_weight_, event);

            if(object_selector_.clean_jets().size() == 4)
                dir->second["full_discriminant_4j" ].fill(hyp.Prob(), evt_weight_);
            else if(object_selector_.clean_jets().size() == 5)
                dir->second["full_discriminant_5j" ].fill(hyp.Prob(), evt_weight_);
            else
                dir->second["full_discriminant_6Pj"].fill(hyp.Prob(), evt_weight_);

            //choose up/down
            TLorentzVector whad(hyp.WHad());
            TLorentzVector j1(*hyp.WJa());
            TLorentzVector j2(*hyp.WJb());
            j1.Boost(whad.BoostVector()*-1);
            j2.Boost(whad.BoostVector()*-1);
            if(j2.E() > j1.E()) hyp.SwapWJets();

            //Angular variables
            hyp::TTbar ttang(hyp);
            auto ttcm = ttang.to_CM(); 
            auto tlepcm = ttang.tlep().to_CM();
            auto thadcm = ttang.thad().to_CM();

            //top angles
            double tlep_ctstar = min(fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())), 0.99999);
            double thad_ctstar = min(fabs(ttang.unit3D().Dot(ttcm.thad().unit3D())), 0.99999);
            dir->second["tlep_ctstar"].fill(tlep_ctstar, evt_weight_);
            dir->second["thad_ctstar"].fill(thad_ctstar, evt_weight_);
            //dir->second["tlep_ctstar"].fill(ttang.unit3D().Dot(ttcm.tlep().unit3D()), evt_weight_);
            //dir->second["thad_ctstar"].fill(ttang.unit3D().Dot(ttcm.thad().unit3D()), evt_weight_);
            dir->second["mtt_tlep_ctstar"].fill(
                    hyp.LVect().M(), tlep_ctstar, evt_weight_ 
                    //fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())), 
                    //evt_weight_
                    );
            if(pdf) pdf_uncs_.fill_replicas2D(
                    folder, "mtt_tlep_ctstar", 
                    hyp.LVect().M(), tlep_ctstar,
                    //hyp.LVect().M(), fabs(ttang.unit3D().Dot(ttcm.tlep().unit3D())),
                    evt_weight_, event
                    );

            dir->second["cdelta_ld"].fill(tlepcm.W().l().unit3D().Dot(thadcm.W().down().unit3D()), evt_weight_);

            //helicity frame
            dir->second["hframe_ctheta_d"].fill(
                    thadcm.W().down().unit3D().Dot(ttcm.thad().unit3D()), 
                    evt_weight_
                    );
        }

        void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
        }

        virtual void begin(){
            Logger::log().debug() << "htt_simple::begin" << endl;
            outFile_.cd();
            if(sync_) {
                sync_tree_ = new TTree("sync", "sync");
                sync_tree_->Branch("Run", &sync_info_.Run, "Run/i");
                sync_tree_->Branch("LumiSection", &sync_info_.LumiSection, "LumiSection/l");
                sync_tree_->Branch("Event", &sync_info_.Event, "Event/l");
                sync_tree_->Branch("RecoSuccess", &sync_info_.RecoSuccess, "RecoSuccess/O");
                sync_tree_->Branch("MassTT", &sync_info_.MassTT, "MassTT/f");
                sync_tree_->Branch("hasMuon", &sync_info_.hasMuon, "hasMuon/O");
                return;
            }

            // added for expected event tables
            string exp = "Expected_Plots";
            book<TH1F>(exp, "Expected_Event_Categories_3J", "", 5, 0.5, 5.5);
            book<TH1F>(exp, "Expected_Event_Categories_4PJ", "", 5, 0.5, 5.5);
            
            if( isTTJetsM1000_ || isTTJetsM700_ ){
                book<TH1F>(exp, "Expected_Event_Categories_3J_massbins", "", 11, 0.5, 11.5); // split into 1k-1200 (700-800),12-15 (800-900),15-2k (900-1k) for M1000 (M700)
                book<TH1F>(exp, "Expected_Event_Categories_4PJ_massbins", "", 11, 0.5, 11.5);
            }

        ////////////// book 3 jet plots

            string evt_3J = "3J_Events";
            book_3J_reco_plots(evt_3J+"/Pre_Alpha/Reconstruction");
            book_3J_reso_plots(evt_3J+"/Pre_Alpha/Resolution");
            book_3J_reco_plots(evt_3J+"/Post_Alpha/Reconstruction");
            book_3J_reso_plots(evt_3J+"/Post_Alpha/Resolution");
            //book_presel_plots(evt_3J+"/Preselection");		// not needed, 3 jet events already included
            //book_selection_plots(evt_3J+"/selection", false);

            ////if( isTTbar_ || isSignal_ ){
            //book<TH1F>(evt_3J, "njets", "", 10, 2., 12.);
            ////}

        ////////////// book 3 jet plots

            vector<string> leptons = {"electrons", "muons"};		
            vector<string> subs;
            if(isSignal_) subs = {"/positive", "/negative"};
            else subs = {"/right", "/matchable", "/unmatchable", "/noslep"};		
            vector<string> lepIDs  = {"looseNOTTight", "tight"};		
            //vector<string> MTs     = {"MTHigh"};		
            vector<string> MTs     = {"MTHigh", "MTLow"};		
            //vector<string> tagging = {"ctagged", "notag"};
            for(auto& lepton : leptons) {
                for(auto& sys : systematics_) {
                    string sys_name = systematics::shift_to_name.at(sys);
                    string dname = lepton+"/"+sys_name+"/preselection";
                    //Logger::log().debug() << "Booking histos in: " << dname << endl;
                        // 3 and 4+ jet plots
                    book_presel_plots(dname);				
                    book_combo_plots(dname);
                        // 4+ jet plots
                    book_presel_plots("4PJ/selection/"+dname);				
                    book_combo_plots( "4PJ/selection/"+dname);
                        // 3 jet plots
                    book_presel_plots(evt_3J+"/selection/"+dname);				
                    book_combo_plots( evt_3J+"/selection/"+dname);
                    for(auto& subsample : subs) {
                        string sub = (isTTbar_ || isSignal_) ? subsample : "";
                        for(auto& lepid : lepIDs) {
                            for(auto& mt : MTs) {
                                if(sys != Sys::NOSYS && (lepid != "tight" || mt != "MTHigh")) {
                                    continue;
                                }
                                stringstream dstream;
                                dstream << lepton << sub << "/";
                                dstream << sys_name << "/" << lepid << "/" << mt;

                                //Logger::log().debug() << "Booking histos in: " << dstream.str() << endl;
                                bool runpdf = (runsys_ && isTTbar_ && sys == systematics::SysShifts::NOSYS && lepid == "tight" && mt == "MTHigh");
                                //bool runpdf = (runsys_ && isTTbar_ );//&& sys == systematics::SysShifts::NOSYS && lepid == "tight" && mt == "MTHigh");
                                    // 3 and 4+ jet plots
                                book_selection_plots(dstream.str(), runpdf);
                                    // 4+ jet plots
                                book_selection_plots("4PJ/selection/"+dstream.str(), runpdf);
                                    // 3 jet plots
                                book_selection_plots(evt_3J+"/selection/"+dstream.str(), runpdf);
                            }
                        }
                        if(!(isTTbar_ || isSignal_)) break;
                    }
                }
            }
        }


    ////// functions for events with 3 jets
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

            // apply alpha correction
        std::pair< TLorentzVector, double > alpha_thad( Permutation &perm ){
                //alphas found by fitting Alpha_THad_P/E hists
                    //values taken from 1degree vals http://home.fnal.gov/~jdulemba/Plots/ttbar_reco_3J/2018/JetpTcut30/Only_Alpha_Correction/Full/ttJetsM0/3J_Event_Plots/Lost_BP/Alpha_Correction/fit_parameters.json
            // y = mx +b -> m is first element in json, b is second

            double alpha_E = 0.4019*( 173.1/perm.THad().M() ) + 0.5834; // only alpha_E used because it's more consistent over mtt spectrum 
            //double alpha_P = 0.1544*( 173.1/perm.THad().M() ) + 0.8599;
            TLorentzVector Alpha_THad(alpha_E*perm.THad().Px(), alpha_E*perm.THad().Py(), alpha_E*perm.THad().Pz(), alpha_E*perm.THad().E());

                // corrected thad
            TLorentzVector reco_corr_ttang = Alpha_THad+perm.TLep();
            TLorentzVector reco_ttcm_corr_thad = Alpha_THad;
            reco_ttcm_corr_thad.Boost(-1*reco_corr_ttang.BoostVector());
            double reco_corr_thad_cth = reco_corr_ttang.Vect().Unit().Dot(reco_ttcm_corr_thad.Vect().Unit());

            return std::make_pair( Alpha_THad, reco_corr_thad_cth ); // corrected thad tvec, costh* val
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



        // find best perm assuming it's lost 
        Permutation lost_best_perm( Permutation &perm, string disc_type, string presel_dir, bool lep_is_tight ){

            Permutation lost_bp;
            double lost_lowest_discval = 1e10;

            for( auto test_perm : permutator_.permutations_3J(perm.WJa(), perm.WJb(), perm.BHad(), perm.BLep(), perm.L(), perm.MET(), perm.LepCharge()) ){
                solver_.Solve_3J_Lost(test_perm);

                if(!sync_ && lep_is_tight){
                    fill_3J_combo_plots("3J_Events/selection/"+presel_dir, test_perm);
                    fill_3J_combo_plots(presel_dir, test_perm);
                }
                double tp_discval = 3e10;

                if( disc_type == "Total" ){ tp_discval = test_perm.Prob(); }
                else if( disc_type == "NS" ){ tp_discval = test_perm.NuDiscr(); }
                else if( disc_type == "Mass" ){ tp_discval = test_perm.MassDiscr(); }
                else{
                    Logger::log().debug() << "Not a valid discriminant type" << endl;
                    throw 42;
                }

                if( tp_discval < lost_lowest_discval ){
                    lost_lowest_discval = tp_discval;
                    lost_bp = test_perm;
                }
            }

            return lost_bp;
        }



        //Permutation best_perm_3J(URStreamer &event){
        //    //initialize perm objects
        //    //jets
        //    IDJet* wj1 = 0;
        //    IDJet* wj2 = 0;
        //    IDJet* bj1 = 0;
        //    IDJet* bj2 = 0;
        //    Permutation empty_perm; // perm.WJa(), WJb(), BHad(), BLep()
        //
        //    vector<IDJet*> jets_vector;
        //    for( auto jet : object_selector_.clean_jets() ){
        //        jets_vector.push_back(jet);
        //    }
        //    sort(jets_vector.begin(), jets_vector.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });

        //    //tracker_.track("All Num BTag");

        //    //if( !jets_vector[1]->BTagId(IDJet::BTag::CSVMEDIUM) ) return std::make_pair(empty_perm, empty_perm); // require at least 2 jets that pass btag
        //    if( !jets_vector[1]->BTagId(IDJet::BTag::CSVMEDIUM) ) return empty_perm; // require at least 2 jets that pass btag

        //    //tracker_.track("Num BTag >= 2");

        //    bj1 = jets_vector[0];
        //    bj2 = jets_vector[1];
        //    wj1 = jets_vector[2];

        //    //cout << "bj1 pT: " << bj1->Pt() << endl;
        //    //cout << "bj2 pT: " << bj2->Pt() << endl;
        //    //cout << "wj1 pT: " << wj1->Pt() << endl;

        //    // create perm from the event and find best perms from merged/lost TTBarSolver
        //    Permutation event_perm(wj1, wj2, bj1, bj2, object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wjat, wjb, bhad, blep
        //    //Permutation merged_bp = merged_best_perm( event_perm );
        //    Permutation lost_bp = lost_best_perm( event_perm, "NS" );

        //    return lost_bp;
        //    //return std::make_pair(merged_bp, lost_bp);
        //} // end of best_perm_3J()


        //void lost_bp_cats( GenTTBar &ttbar, Permutation &lost_bp, string bp_solution, URStreamer &event ){

        //    //auto evt_dir = histos_.find("3J_Events");
        //    string lost_perm_status;

        //    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){ // skip to next event if perm is empty
        //        tracker_.track("Not semilep events");
        //        return;
        //    }
        //    tracker_.track("semilep");

        //    // get matched perm from event
        //    Permutation mp = dr_matcher_.dr_match(
        //            genp_selector_.ttbar_final_system(),
        //            object_selector_.clean_jets(),
        //            object_selector_.lepton(),
        //            object_selector_.met(),
        //            object_selector_.lepton_charge());


        //    //if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) || !mp.Lost_Event() ) lost_perm_status = "WRONG";
        //    //if( mp.Lost_Event() ){
        //    //    if( mp.AreBsSame(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "RIGHT";
        //    //    else if( mp.AreBsFlipped(lost_bp) && ( lost_bp.WJa() == mp.WJa() || lost_bp.WJa() == mp.WJb() ) ) lost_perm_status = "LOST_SWAP";
        //    //    else lost_perm_status = "LOST";
        //    //}
        //    //else lost_perm_status = "WRONG";

        //    if( mp.IsEmpty() || ( !(mp.BHad() && mp.BLep()) && !(mp.WJa() || mp.WJb()) ) ) lost_perm_status = "OTHER"; // check if mp exists, has bhad and blep and at least one wjet
        //        // check if mp and lost_bp misidentified b's
        //    else if( mp.AreBsSame(lost_bp)  ) lost_perm_status = "CORRECT_B";
        //    else if( mp.AreBsFlipped(lost_bp) ) lost_perm_status = "WRONG_B";
        //    else lost_perm_status = "OTHER";

        //    string gen_match_status;
        //    if( lost_bp.BHad()->match() == ttbar.had_b() && lost_bp.BLep()->match() != ttbar.lep_b() ) gen_match_status = "CORRECT_BHAD"; // only perm bhad match is same as gen bhad
        //    else if( lost_bp.BHad()->match() != ttbar.had_b() && lost_bp.BLep()->match() == ttbar.lep_b() ) gen_match_status = "CORRECT_BLEP"; // only perm blep match is same as gen blep
        //    else if( lost_bp.BHad()->match() == ttbar.had_b() && lost_bp.BLep()->match() == ttbar.lep_b() ) gen_match_status = "CORRECT_Bs"; // both perm b matches same as gen bs
        //    else if( lost_bp.BHad()->match() == ttbar.lep_b() && lost_bp.BLep()->match() == ttbar.had_b() ) gen_match_status = "SWAPPED_Bs"; // perm bhad match same as gen blep, perm blep match same as gen bhad
        //    else gen_match_status = "OTHER_MATCH";

        //    //evt_dir->second["njets"].fill(object_selector_.clean_jets().size());

        //    fill_3J_pre_alpha_reco_plots("3J_Events/Pre_Alpha/Reconstruction", lost_bp);
        //    fill_3J_pre_alpha_reso_plots("3J_Events/Pre_Alpha/Resolution", lost_bp, ttbar);

        //    fill_3J_post_alpha_reco_plots("3J_Events/Post_Alpha/Reconstruction", lost_bp);
        //    fill_3J_post_alpha_reso_plots("3J_Events/Post_Alpha/Resolution", lost_bp, ttbar);

        //    //fill_presel_plots("3J_Events/Preselection", event);
        //    fill_3J_selection_plots("3J_Events/selection", event, lost_bp, false);

        //} // end of lost_bp_cats

        //void process_3J_evt( systematics::SysShifts shift, URStreamer &event, vector<IDJet*> &clean_jets ){

        //    IDJet* wj2 = 0;
        //    Permutation event_perm(clean_jets[2], wj2, clean_jets[0], clean_jets[1], object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wjat, wjb, bhad, blep
        //    Permutation best_permutation = lost_best_perm( event_perm, "Total" );
        //   
        //    bool reco_success = !best_permutation.IsEmpty() && best_permutation.Prob() < 1e9;
        //    if( !reco_success ){
        //    //if(!sync_ && !reco_success){
        //        //tracker_.track("not sync and not reco success");
        //        return;
        //    }
        //    //if(lep_is_tight) tracker_.track("best perm");

        //    //if(sync_) {
        //    //    if(lep_is_tight) {
        //    //        sync_info_.Run  = event.run;
        //    //        sync_info_.LumiSection = event.lumi;
        //    //        sync_info_.Event= event.evt;
        //    //        sync_info_.hasMuon = (object_selector_.lepton_type() == -1) ? 0 : 1;
        //    //        sync_info_.RecoSuccess = reco_success;
        //    //        if(reco_success)
        //    //            sync_info_.MassTT = best_permutation.LVect().M();
        //    //        else
        //    //            sync_info_.MassTT = -1;
        //    //        sync_tree_->Fill();
        //    //        tracker_.track("SYNC END", leptype);
        //    //    }
        //    //    return;
        //    //}


        //    //if( isTTbar_ || isSignal_){
        //    //    permutator_.reset_3J();
        //    //
        //    //    //generator selection
        //    //    bool selection = genp_selector_.select(event);
        //    //    if( !selection ){
        //    //        Logger::log().debug() << "event has no gen selection " << endl;
        //    //        return;
        //    //    }
        //    //    GenTTBar &ttbar = genp_selector_.ttbar_system();
        //
        //    //    //Permutation bp = best_perm_3J(event);
        //
        //    //    if( !bp.IsEmpty() ){
        //    //        lost_bp_cats( ttbar, bp, "3J_Events", event );
        //
        ////  //          //GenTTBar &ttbar = genp_selector_.ttbar_system();
        ////  //  
        ////  //          exp_dir->second["Expected_Event_Categories_3J"].fill(1); // total expected events == 1
        ////  //          if( isTTJetsM1000_ ){
        ////  //              if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(1); // total expected events == 1
        ////  //              else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(5); // total expected events == 5
        ////  //              else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(9); // total expected events == 1
        ////  //          }
        ////  //          if( isTTJetsM700_ ){
        ////  //              if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(1); // total expected events == 1
        ////  //              else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(5); // total expected events == 5
        ////  //              else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(9); // total expected events == 1
        ////  //          }
        ////  //  
        ////  //          if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){
        ////  //              exp_dir->second["Expected_Event_Categories_3J"].fill(5); // expected other events == 5
        ////  //              if( isTTJetsM1000_ ){
        ////  //                  if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(3); // expected other events == 1
        ////  //                  else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(7); // expected other events == 5
        ////  //                  else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(11); // expected other events == 1
        ////  //              }
        ////  //              if( isTTJetsM700_ ){
        ////  //                  if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(3); // expected other events == 1
        ////  //                  else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(7); // expected other events == 5
        ////  //                  else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(11); // expected other events == 1
        ////  //              }
        ////  //              return;
        ////  //          }
        ////  //          else if( ttbar.type == GenTTBar::DecayType::SEMILEP ) {
        ////  //              if( ttbar.partial_hadronic_merged(0.4) ){ // gen partons merged
        ////  //                  exp_dir->second["Expected_Event_Categories_3J"].fill(3); // expected merged events == 3
        ////  //                  if( isTTJetsM1000_ ){
        ////  //                      if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(2); // expected merged events == 1
        ////  //                      else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(6); // expected merged events == 5
        ////  //                      else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(10); // expected merged events == 1
        ////  //                  }
        ////  //                  if( isTTJetsM700_ ){
        ////  //                      if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(2); // expected merged events == 1
        ////  //                      else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(6); // expected merged events == 5
        ////  //                      else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(10); // expected merged events == 1
        ////  //                  }
        ////  //              }
        ////  //              else{ // gen partons not merged
        ////  //                  exp_dir->second["Expected_Event_Categories_3J"].fill(5); // expected other events == 5
        ////  //                  if( isTTJetsM1000_ ){
        ////  //                      if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(3); // expected other events == 1
        ////  //                      else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(7); // expected other events == 5
        ////  //                      else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(11); // expected other events == 1
        ////  //                  }
        ////  //                  if( isTTJetsM700_ ){
        ////  //                      if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(3); // expected other events == 1
        ////  //                      else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(7); // expected other events == 5
        ////  //                      else exp_dir->second["Expected_Event_Categories_3J_massbins"].fill(11); // expected other events == 1
        ////  //                  }
        ////  //              }
        ////  //          }
        //    //    }
        //    //} // isTTbar or isSignal
        //} // process_3J_evt 




        void process_evt(systematics::SysShifts shift, URStreamer &event){
            tracker_.track("start", "electrons");
            tracker_.track("start", "muons");
            //float weight = 1.;

            auto exp_dir = histos_.find("Expected_Plots");
            //tracker_.track("before reco");

            //select reco objects
            if( !object_selector_.select(event, shift, sync_) ) return;
            //tracker_.track("b4 jet sel");

            int njets = object_selector_.clean_jets().size();		

            //tracker_.track("b4 leptype");

            string leptype = (object_selector_.lepton_type() == -1) ? "electrons" : "muons";
            //tracker_.track("b4 l_is_tight");

            bool lep_is_tight = (object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTMU || 
                    object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTEL);

            //tracker_.track("b4 obj sel", leptype);
            if(lep_is_tight) {
                tracker_.group(leptype);
                tracker_.track("object selection", leptype);
            }

            // check MT
            double mt = MT(object_selector_.lepton(), object_selector_.met());
            bool mt_high;
            if( mt < 50. ) mt_high = false;
            else mt_high = true;

            // cout << "Muon: " << *object_selector_.lepton() << endl;
            //MC Weight for lepton selection
            float lep_weight=1;
            if(!isData_) { 
                if(object_selector_.tight_muons().size() == 1) {
                    lep_weight = muon_sf_.get_sf(object_selector_.muon()->Pt(), object_selector_.muon()->Eta());
                    // if(-1.3 < object_selector_.muon()->Eta() && object_selector_.muon()->Eta() < -1)
                    // 	cout << "Mu " << *object_selector_.muon() << " weight: " << lep_weight << " prev: " << evt_weight_ << endl;
                }
                if(object_selector_.tight_electrons().size() == 1) lep_weight = electron_sf_.get_sf(object_selector_.electron()->Pt(), object_selector_.electron()->etaSC());
            }
            evt_weight_ *= lep_weight;
            string sys_name = systematics::shift_to_name.at(shift);
            stringstream presel_dir;
            presel_dir << leptype << "/";
            presel_dir << sys_name << "/preselection";
            //if(!sync_ && lep_is_tight ){
            if(!sync_ && lep_is_tight && mt_high){
                fill_presel_plots(presel_dir.str(), event);
                //tracker_.track("not sync and tight lep");
                if( njets == 3 ) fill_presel_plots("3J_Events/selection/"+presel_dir.str(), event);
                else fill_presel_plots("4PJ/selection/"+presel_dir.str(), event);
            }


            //cut on btag
            auto &clean_jets = object_selector_.clean_jets();
                // original b cuts
            //sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->CombinedMVA() > B->CombinedMVA());});
            sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){ return( A->csvIncl() > B->csvIncl() ); });
            if(!clean_jets[0]->BTagId(cut_tight_b_)) return;
            //tracker_.track("first b pass");
            if(!clean_jets[1]->BTagId(cut_loose_b_)) return;
            //tracker_.track("second b pass");
                //

            //if( !clean_jets[0]->BTagId(IDJet::BTag::CSVMEDIUM) ) return; // require at least 2 jets that pass btag
            //if( !clean_jets[1]->BTagId(IDJet::BTag::CSVMEDIUM) ) return; // require at least 2 jets that pass btag


            if(lep_is_tight) tracker_.track("b cuts", leptype);

            if( !permutator_.preselection(
                        object_selector_.clean_jets(), object_selector_.lepton(), 
                        object_selector_.met(), object_selector_.lepton_charge(),
                        //lep_is_tight) ) return;
                event.rho().value(), lep_is_tight) ) return;
            if(lep_is_tight) tracker_.track("perm preselection");

            //find mc weight for btag
            double bweight = 1;
            if(!isData_ && !sync_){
                bweight *= btag_sf_.scale_factor(clean_jets, shift, false);
                //tracker_.track("not data and not sync");
            }
            evt_weight_ *= bweight;
            //find mc weight for top pt for ttbar events
            if( isTTbar_ ){
                double top_pt_weight = top_pt_reweighting("nominal");
                evt_weight_ *= top_pt_weight;
            }

            if( njets == 3 ){

                /// Find best permutation
                permutator_.reset_3J();
                IDJet* wj2 = 0;
                Permutation event_perm(clean_jets[2], wj2, clean_jets[0], clean_jets[1], object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wjat, wjb, bhad, blep
                Permutation best_perm = lost_best_perm( event_perm, "Total", presel_dir.str(), lep_is_tight );
           
                bool reco_success_3J = !best_perm.IsEmpty() && best_perm.Prob() < 1e9;
                if(!sync_ && !reco_success_3J){
                    //tracker_.track("not sync and not reco success");
                    return;
                }
                if(lep_is_tight) tracker_.track("best perm");

                if(sync_) {
                    if(lep_is_tight) {
                        sync_info_.Run  = event.run;
                        sync_info_.LumiSection = event.lumi;
                        sync_info_.Event= event.evt;
                        sync_info_.hasMuon = (object_selector_.lepton_type() == -1) ? 0 : 1;
                        sync_info_.RecoSuccess = reco_success_3J;
                        if(reco_success_3J){
                                // corrected thad
                            std::pair< TLorentzVector, double > corrected_thad = alpha_thad(best_perm); // corrected thad tvec, costh* value
                            sync_info_.MassTT = (corrected_thad.first + best_perm.TLep()).M();
                        }
                        else{
                            sync_info_.MassTT = -1;
                        }
                        sync_tree_->Fill();
                        tracker_.track("SYNC END", leptype);
                    }
                    return;
                }


                //create event dir (contains the dir path to be filled)
                stringstream evtdir;		
                evtdir << leptype;

                //Gen matching (TT events only)
                Permutation mp_3J;
                if(isTTbar_) {
                    //generator selection
                    GenTTBar &ttbar = genp_selector_.ttbar_system();

                    if( ttbar.type == GenTTBar::DecayType::SEMILEP ) {
                        mp_3J = dr_matcher_.dr_match(
                            genp_selector_.ttbar_final_system(),
                            object_selector_.clean_jets(),
                            object_selector_.lepton(),
                            object_selector_.met(),
                            object_selector_.lepton_charge());
                        mp_3J.SetMET(object_selector_.met());

                        if( best_perm.AreBsSame(mp_3J) && ( best_perm.WJa() == mp_3J.WJa() || best_perm.WJa() == mp_3J.WJb() ) ) evtdir << "/right";
                        else if( mp_3J.IsTLepComplete() && mp_3J.BHad() && ( mp_3J.WJa() || mp_3J.WJb() ) ) evtdir << "/matchable";
                        else evtdir << "/unmatchable";
                    }
                    else evtdir << "/noslep";
                } 
                else if (isSignal_) {
                    if(evt_weight_ > 0.) evtdir << "/positive";
                    else evtdir << "/negative";
                }

                if(lep_is_tight) tracker_.track("matched perm");

                //cout << "category: " << evtdir.str() << endl;

                //check isolation type
                bool runpdf = (runsys_ && isTTbar_ && shift == systematics::SysShifts::NOSYS);
                evtdir << "/" << sys_name << "/";
                bool tight;
                switch(object_selector_.event_type()) {
                    case TTObjectSelector::EvtType::TIGHTMU: 
                    case TTObjectSelector::EvtType::TIGHTEL: evtdir << "tight"; runpdf &= true; tight=true; break;
                    case TTObjectSelector::EvtType::LOOSEMU: 
                    case TTObjectSelector::EvtType::LOOSEEL: evtdir << "looseNOTTight"; runpdf &= false; tight=false; break;
                    default: throw 42; break;
                }
                if(!tight && shift != Sys::NOSYS) return;

                //MT category
                if( shift == Sys::NOSYS ){
                    if( !mt_high ){
                        evtdir << "/MTLow";
                        runpdf=false;
                        //cout << "MTLow: " << mt_high << endl;
                    }
                    else evtdir << "/MTHigh";
                }
                else{
                    if( !mt_high ) return;
                    evtdir << "/MTHigh";
                }
                //evtdir << "/MTHigh";

                //cout << "category: " << evtdir.str() << endl;
                //if( evtdir.str().find("right") ) cout << "category: " << evtdir.str() << endl;

                //fill right category
                fill_3J_selection_plots("3J_Events/selection/"+evtdir.str(), event, best_perm, runpdf);
                fill_3J_selection_plots(evtdir.str(), event, best_perm, runpdf);
                if(lep_is_tight) tracker_.track("END", leptype);


                if( shift == systematics::SysShifts::NOSYS ){
                    if(isTTbar_) {
                        //generator selection
                        GenTTBar &ttbar = genp_selector_.ttbar_system();

                        if( ttbar.type == GenTTBar::DecayType::SEMILEP ) {
                            fill_3J_pre_alpha_reso_plots("3J_Events/Pre_Alpha/Resolution", best_perm, ttbar);
                            fill_3J_post_alpha_reso_plots("3J_Events/Post_Alpha/Resolution", best_perm, ttbar);
                        }
                    } 
                    fill_3J_pre_alpha_reco_plots("3J_Events/Pre_Alpha/Reconstruction", best_perm);
                    fill_3J_post_alpha_reco_plots("3J_Events/Post_Alpha/Reconstruction", best_perm);
                }


                return;
            } // njets == 3 requirement



    //// events with 4+ jets

            // cout << object_selector_.clean_jets().size() << " --> " << permutator_.capped_jets().size() << endl;
            // for(auto i : object_selector_.clean_jets()) {
            // 	cout << *i << endl;
            // }
            // cout << endl << endl;

            //Find best permutation
            Permutation best_permutation;
            //int idx = 0;
            for(auto& test_perm : permutator_.permutations()) {
                solver_.Solve(test_perm);

                // cout << "Permutation #" << idx << "  " << test_perm.Prob() << " -- " << best_permutation.Prob() << endl;
                // idx++;
                // if(test_perm.IsComplete()) {
                // 	cout << test_perm << endl;
                // }
                // else {
                // 	cout << "Permutation incomplete!" << endl
                // 			 << "l: " << test_perm.L() << ", b_l: " << test_perm.BLep() << ", b_h: " 
                // 			 << test_perm.BHad() << ", j1: " << test_perm.WJa() << ", j2: " 
                // 			 << test_perm.WJb() << endl;
                // }
                if(!sync_ && lep_is_tight){
                    fill_combo_plots(presel_dir.str(), test_perm);
                    fill_combo_plots("4PJ/selection/"+presel_dir.str(), test_perm);
                }
                if(test_perm.Prob()  < best_permutation.Prob()) {
                    best_permutation = test_perm;
                    // cout << "This is the best so far!" << endl; 
                }
                // cout << endl;
            }

            //            float wja_sf = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, best_permutation.WJa()->Eta()}});
            //            float wja_sf_down = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, best_permutation.WJa()->Eta()}}, Variation::DOWN);
            //            float wja_sf_up = jet_res_sf_.getScaleFactor({{JME::Binning::JetEta, best_permutation.WJa()->Eta()}}, Variation::UP);

            bool reco_success = (best_permutation.IsComplete() && best_permutation.Prob() <= 1E9);
            if(!sync_ && !reco_success){
                //tracker_.track("not sync and not reco success");
                return;
            }
            if(lep_is_tight) tracker_.track("best perm");

            //	    	cout << "BEST PERM" << endl;
            //	    	if(best_permutation.IsComplete()){
            //	    		cout << best_permutation << endl;
            //	    	} 
            //	    	else cout << "Permutation incomplete!" <<endl;
            //            cout << "wja pt res: " << wja_pt_res << ", sf: " << wja_sf << ", sf down: " << wja_sf_down << ", sf up: " << wja_sf_up << endl;


            if(sync_) {
                if(lep_is_tight) {
                    sync_info_.Run  = event.run;
                    sync_info_.LumiSection = event.lumi;
                    sync_info_.Event= event.evt;
                    sync_info_.hasMuon = (object_selector_.lepton_type() == -1) ? 0 : 1;
                    sync_info_.RecoSuccess = reco_success;
                    if(reco_success)
                        sync_info_.MassTT = best_permutation.LVect().M();
                    else
                        sync_info_.MassTT = -1;
                    sync_tree_->Fill();
                    tracker_.track("SYNC END", leptype);
                }
                return;
            }

            //create event dir (contains the dir path to be filled)
            stringstream evtdir;		
            evtdir << leptype;

            //Gen matching (TT events only)
            Permutation matched_perm;
            if(isTTbar_) {
                //generator selection
                GenTTBar &ttbar = genp_selector_.ttbar_system();

                if( ttbar.type == GenTTBar::DecayType::SEMILEP ) {
                    matched_perm = matcher_.match(
                            genp_selector_.ttbar_final_system(),
                            object_selector_.clean_jets(), 
                            object_selector_.veto_electrons(),
                            object_selector_.veto_muons()
                            );
                    matched_perm.SetMET(object_selector_.met());

                    if(best_permutation.IsCorrect(matched_perm)) evtdir << "/right";
                    else if(matched_perm.IsComplete()) evtdir << "/matchable";
                    else evtdir << "/unmatchable";
                }
                else evtdir << "/noslep";
            } 
            else if (isSignal_) {
                if(evt_weight_ > 0.) evtdir << "/positive";
                else evtdir << "/negative";
            }

            if( isTTbar_ || isSignal_ ){
                GenTTBar &ttbar = genp_selector_.ttbar_system();
                bool selection = genp_selector_.select(event);
                if( selection ){
                    exp_dir->second["Expected_Event_Categories_4PJ"].fill(1); // total expected events == 1
                    if( isTTJetsM1000_ ){
                        if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(1); // total expected events == 1
                        else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(5); // total expected events == 5
                        else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(9); // total expected events == 1
                    }
                    if( isTTJetsM700_ ){
                        if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(1); // total expected events == 1
                        else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(5); // total expected events == 5
                        else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(9); // total expected events == 1
                    }

                    if( !(ttbar.type == GenTTBar::DecayType::SEMILEP) ){
                        exp_dir->second["Expected_Event_Categories_4PJ"].fill(5); // expected other events == 5
                        if( isTTJetsM1000_ ){
                            if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(3); // expected other events == 1
                            else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(7); // expected other events == 5
                            else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(11); // expected other events == 1
                        }
                        if( isTTJetsM700_ ){
                            if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(3); // expected other events == 1
                            else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(7); // expected other events == 5
                            else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(11); // expected other events == 1
                        }
                        //return;
                    }
                    else if( ttbar.type == GenTTBar::DecayType::SEMILEP ) {
                        if( ttbar.partial_hadronic_merged(0.4) ){ // gen partons merged
                            exp_dir->second["Expected_Event_Categories_4PJ"].fill(3); // expected merged events == 3
                            if( isTTJetsM1000_ ){
                                if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(2); // expected merged events == 1
                                else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(6); // expected merged events == 5
                                else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(10); // expected merged events == 1
                            }
                            if( isTTJetsM700_ ){
                                if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(2); // expected merged events == 1
                                else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(6); // expected merged events == 5
                                else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(10); // expected merged events == 1
                            }
                        }
                        else{ // gen partons not merged
                            exp_dir->second["Expected_Event_Categories_4PJ"].fill(5); // expected other events == 5
                            if( isTTJetsM1000_ ){
                                if( ttbar.M() < 1200 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(3); // expected other events == 1
                                else if( ttbar.M() >= 1200 && ttbar.M() < 1500 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(7); // expected other events == 5
                                else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(11); // expected other events == 1
                            }
                            if( isTTJetsM700_ ){
                                if( ttbar.M() < 800 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(3); // expected other events == 1
                                else if( ttbar.M() >= 800 && ttbar.M() < 900 ) exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(7); // expected other events == 5
                                else exp_dir->second["Expected_Event_Categories_4PJ_massbins"].fill(11); // expected other events == 1
                            }
                        }
                    }
                }
            }

            if(lep_is_tight) tracker_.track("matched perm");

            //check isolation type
            bool runpdf = (runsys_ && isTTbar_ && shift == systematics::SysShifts::NOSYS);
            evtdir << "/" << sys_name << "/";
            bool tight;
            switch(object_selector_.event_type()) {
                case TTObjectSelector::EvtType::TIGHTMU: 
                case TTObjectSelector::EvtType::TIGHTEL: evtdir << "tight"; runpdf &= true; tight=true; break;
                case TTObjectSelector::EvtType::LOOSEMU: 
                case TTObjectSelector::EvtType::LOOSEEL: evtdir << "looseNOTTight"; runpdf &= false; tight=false; break;
                default: throw 42; break;
            }
            if(!tight && shift != Sys::NOSYS) return;

            //MT category
            if( shift == Sys::NOSYS ){
                if( !mt_high ){
                    evtdir << "/MTLow";
                    runpdf=false;
                }
                else evtdir << "/MTHigh";
            }
            else{
                if( !mt_high ) return;
                evtdir << "/MTHigh";
            }
            //evtdir << "/MTHigh";

            //fill right category
            fill_selection_plots(evtdir.str(), event, best_permutation, runpdf);
            fill_selection_plots("4PJ/selection/"+evtdir.str(), event, best_permutation, runpdf);
            if(lep_is_tight) tracker_.track("END", leptype);
            //           } // end of 4+ jet loop

        } // end of process_evt

        //This method is called once every file, contains the event loop
        //run your proper analysis here
        virtual void analyze(){
            opts::variables_map &values = URParser::instance().values();
            //int limit = 10;
            int limit = values["limit"].as<int>();
            int report = values["report"].as<int>();
            int skip  = values["skip"].as<int>();
            string pick = values["pick"].as<string>();
            if(values.count("bystep")) tracker_.verbose(true);
            EventList picker;
            if(pick.size() != 0) {
                EventList nn(pick);
                picker = nn;
            }

            if(evt_idx_ >= limit) return;
            Logger::log().debug() << "htt_simple::analyze" << endl;
            Logger::log().debug() << tree_ << " " << tree_->GetEntries() << endl;
            URStreamer event(tree_);

            while( event.next() )
            {

                if(picker.active()){
                    if(picker.contains(event.run, event.lumi, event.evt)) {
                        Logger::log().fatal() << "Picking event " << " run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
                    }
                    else continue;
                }
                evt_idx_++;
                if(skip > 0 && evt_idx_ < skip) {
                    continue;
                }
                evt_idx_--;
                if(evt_idx_ % report == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << " run: " << event.run << " lumisection: " << event.lumi << " eventnumber: " << event.evt << endl;
                if(limit > 0 && evt_idx_ > limit) {
                    return;
                }
                evt_idx_++;

                if(isTTbar_){
                    genp_selector_.select(event);           
                }

                for(auto shift : systematics_){
                    evt_weight_ = (isData_) ? 1. : mc_weights_.evt_weight(event, shift);

                    if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
                    //Logger::log().debug() << "processing: " << shift << endl;
                    process_evt(shift, event);
                    tracker_.group("");
                    if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
                }
            } //while(event.next())

            Logger::log().debug() << "End of analyze() " << evt_idx_ << endl;
        } // end of analyze()

        //this method is called at the end of the job, by default saves
        //every histogram/tree produced, override it if you need something more
        virtual void end(){
            outFile_.cd();
            tracker_.writeTo(outFile_);
            if(sync_) {
                sync_tree_->Write();
                return;
            }
            outFile_.Write();
            for(string ltype : {"electrons", "muons"}) {
                TDirectory* tdir = (TDirectory*) outFile_.Get(ltype.c_str());
                tracker_.writeTo(*tdir, ltype);
            }
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
                ("nopdfs", opts::value<int>()->default_value(0), "do not run pdf uncertainties")
                ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                ("report,r", opts::value<int>()->default_value(10000), "limit the number of events processed per file")
                ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                ("optimization", "run in optimization mode: reduced number of histograms")
                ("sync", "dump sync ntuple")
                ("bystep", "print every step")
                ("nosys", "do not run systematics")
                ("noweights", "do not run systematics")
                ("pick", opts::value<string>()->default_value(""), "pick from evtlist");

            parser.addCfgParameter<std::string>("general", "csv_sffile", "");
            parser.addCfgParameter<std::string>("general", "wjets_efficiencies", "");
        }
        };

        //make it executable
        int main(int argc, char *argv[])
        {
            URParser &parser = URParser::instance(argc, argv);
            URDriver<htt_simple> test;
            int excode = test.run();
            //Logger::log().debug() << "RUNNING DONE " << std::endl;
            auto files = gROOT->GetListOfFiles(); //make ROOT aware that some files do not exist, because REASONS
            Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl; //need to print out this otherwise ROOT loses its shit in 7.4.X (such I/O, much features)
            opts::variables_map &values = parser.values();
            string output_file = values["output"].as<std::string>();
            Logger::log().debug() << "  Output File: " << output_file << endl;
            return excode;
        }
