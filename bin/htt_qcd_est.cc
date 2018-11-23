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
#include "Analyses/URTTbar/interface/Alpha_Corrections.h"
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

using namespace TMath;
using namespace systematics;
typedef SysShifts Sys;

class htt_qcd_est : public AnalyzerBase
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

        //selectors and helpers
        TTGenParticleSelector genp_selector_;
        TTGenMatcher matcher_;
        DR_TTGenMatcher dr_matcher_;
        Alpha_Corrections alpha_corr_;
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
        string btag_str_id_;

        //sync
        bool sync_;
        TTree *sync_tree_;
        SyncInfo sync_info_;

        float MTCut_;

    public:
        inline double MT(TLorentzVector *l, TLorentzVector *met) {
            return sqrt(pow(l->Pt() + met->Pt(), 2) - pow(l->Px() + met->Px(), 2) - pow(l->Py() + met->Py(), 2));
        }

        htt_qcd_est(const std::string output_filename):
            AnalyzerBase("htt_qcd_est", output_filename), 
            tracker_(),
            //genp_selector_(TTGenParticleSelector::SelMode::LHE),
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
            pdf_uncs_(254),
            systematics_(),
            sync_(false),
            sync_tree_(0),
            sync_info_(){
                object_selector_.allow_loose();//allow loose but not tight lepton events
                Logger::log().debug() << "htt_qcd_est ctor" << endl;
                cut_tight_b_ = btag_sf_.tight_cut();
                cut_loose_b_ = btag_sf_.loose_cut();
                btag_str_id_ = IDJet::id_string(cut_tight_b_);
                Logger::log().debug() << "btag tight: " << IDJet::tag2string(cut_tight_b_) << ", loose: " << IDJet::tag2string(cut_loose_b_) << " WPs" << endl;
                //  " Perm: " << permutator_.tight_bID_cut() << " " << permutator_.loose_bID_cut() << endl;

                URParser &parser = URParser::instance();

                parser.addCfgParameter<string>("event", "MTCut","");
                parser.parseArguments();

                MTCut_ = parser.getCfgPar<float>("event", "MTCut" );

                if( MTCut_ != 50. ){
                    Logger::log().error() << "MTCut is " << MTCut_ << " but should be 50 in cfg file!" << endl;
                    throw 42;
                }


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

                if( isTTbar_ ) genp_selector_ = TTGenParticleSelector(TTGenParticleSelector::SelMode::LHE);
                else genp_selector_ = TTGenParticleSelector();

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


        //// selection plots


        void fill_3J_selection_plots(string folder, URStreamer &event, Permutation &hyp, bool pdf){
            fill_presel_plots(folder, event);
        } // fill_3J_selection_plots
    ////////////////////// book and fill 3 jet plots



        void book_presel_plots(string folder) {
            if(optim_) return;
            book<TH1F>(folder, "min_lep_jet_dr" , ";min #Delta R(lep, jets);", 50, 0., 5.);			
            book<TH1F>(folder, "lep_pt_B" , ";p_{T}(l matched to b-quark) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "lep_pt_Prompt" , ";p_{T}(l matched to prompt-quark) (GeV)", 500, 0., 500.);
            book<TH1F>(folder, "lep_eta_B" , ";#eta(l matched to b-quark)", 300, -3., 3.);
            book<TH1F>(folder, "lep_eta_Prompt" , ";#eta(l matched to prompt-quark)", 300, -3., 3.);
            book<TH2F>(folder, "lep_eta_pt_B" , ";p_{T}(l matched to b-quark) (GeV);#eta(l matched to b-quark)", 500, 0., 500., 300, -3., 3.);
            book<TH2F>(folder, "lep_eta_pt_Prompt" , ";p_{T}(l matched to prompt-quark) (GeV);#eta(l matched to prompt-quark)", 500, 0., 500., 300, -3., 3.);

            book<TH1F>(folder, "MT" , ";M_{T} (GeV)",  500, 0, 500);

        }

        void fill_presel_plots(string folder, URStreamer &event){
            if(optim_) return;
            auto dir = histos_.find(folder);
            if(dir == histos_.end()) {
                Logger::log().error() << "could not find: " << folder << endl;
                throw 40;
            }

            double mt = MT(object_selector_.lepton(), object_selector_.met());
            dir->second["MT"].fill(mt, evt_weight_);

            double lep_jet_DR = 1e10;
            Jet lep_parton;
            for( auto jet : event.jets() ){
                if( jet.DeltaR( *object_selector_.lepton() ) < lep_jet_DR ){
                    lep_jet_DR = jet.DeltaR( *object_selector_.lepton() );
                    lep_parton = jet;
                }
            }

            int hflav = fabs( lep_parton.hadronFlavour() );
            if( hflav == 5 ){
                dir->second["lep_pt_B"].fill(object_selector_.lepton()->Pt(), evt_weight_);
                dir->second["lep_eta_B"].fill(object_selector_.lepton()->Eta(), evt_weight_);
                dir->second["lep_eta_pt_B"].fill(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta(), evt_weight_);
            }
            else{
                dir->second["lep_pt_Prompt"].fill(object_selector_.lepton()->Pt(), evt_weight_);
                dir->second["lep_eta_Prompt"].fill(object_selector_.lepton()->Eta(), evt_weight_);
                dir->second["lep_eta_pt_Prompt"].fill(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta(), evt_weight_);
            }

            dir->second["min_lep_jet_dr"].fill(lep_jet_DR, evt_weight_);

        }

        void book_selection_plots(string folder, bool pdf){
            book_presel_plots(folder);
        }

        void fill_selection_plots(string folder, URStreamer &event, Permutation &hyp, bool pdf){
            fill_presel_plots(folder, event);
        }

        //void fill(string folder, Permutation &hyp){//, TTbarHypothesis *genHyp=0) {
        //}

        virtual void begin(){
            Logger::log().debug() << "htt_qcd_est::begin" << endl;
            outFile_.cd();
            //if(sync_) {
            //    sync_tree_ = new TTree("sync", "sync");
            //    sync_tree_->Branch("Run", &sync_info_.Run, "Run/i");
            //    sync_tree_->Branch("LumiSection", &sync_info_.LumiSection, "LumiSection/l");
            //    sync_tree_->Branch("Event", &sync_info_.Event, "Event/l");
            //    sync_tree_->Branch("RecoSuccess", &sync_info_.RecoSuccess, "RecoSuccess/O");
            //    sync_tree_->Branch("MassTT", &sync_info_.MassTT, "MassTT/f");
            //    sync_tree_->Branch("hasMuon", &sync_info_.hasMuon, "hasMuon/O");
            //    return;
            //}

            vector<string> njets = {"3Jets", "4PJets"};
            vector<string> leptons = {"electrons", "muons"};		
            vector<string> subs;
            if(isSignal_) subs = {"/positive", "/negative"};
            else subs = {"/right", "/matchable", "/unmatchable", "/noslep"};		
            vector<string> lepIDs  = {"looseNOTTight", "tight"};		
            vector<string> btags   = {"csvPass", "csvFail"};		
            for(auto& lepton : leptons) {
                for(auto& njet: njets){
                    for(auto& sys : systematics_) {
                        string sys_name = systematics::shift_to_name.at(sys);
                        string dname = lepton+"/"+njet+"/"+sys_name+"/preselection";
                        //Logger::log().debug() << "Booking histos in: " << dname << endl;
                            // 3 and 4+ jet plots
                        book_presel_plots(dname);				
                        for(auto& subsample : subs) {
                            string sub = (isTTbar_ || isSignal_) ? subsample : "";
                            for(auto& lepid : lepIDs) {
                                for(auto& btag : btags) {
                                    if(sys != Sys::NOSYS && lepid != "tight"){// || mt != "MTHigh")) {
                                        continue;
                                    }
                                    stringstream dstream;
                                    dstream << lepton << "/" << njet << sub << "/";
                                    dstream << sys_name << "/" << lepid << "/" << btag;

                                    //Logger::log().debug() << "Booking histos in: " << dstream.str() << endl;
                                    bool runpdf = (runsys_ && isTTbar_ && sys == systematics::SysShifts::NOSYS && lepid == "tight");
                                    //bool runpdf = (runsys_ && isTTbar_ );//&& sys == systematics::SysShifts::NOSYS && lepid == "tight" && mt == "MTHigh");
                                        // 3 and 4+ jet plots
                                    book_selection_plots(dstream.str(), runpdf);
                                }
                                
                            }
                            if(!(isTTbar_ || isSignal_)) break;
                        }
                    }
                }
            }
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

        void process_evt(systematics::SysShifts shift, URStreamer &event){
            tracker_.track("start", "electrons");
            tracker_.track("start", "muons");

            //select reco objects
            if( !object_selector_.select(event, shift, sync_) ) return;
            //tracker_.track("b4 jet sel");

            int njets = object_selector_.clean_jets().size();		

            //tracker_.track("b4 leptype");
            string njets_type = ( njets == 3 ) ? "3Jets" : "4PJets";
            string leptype = (object_selector_.lepton_type() == -1) ? "electrons" : "muons";
            //tracker_.track("b4 l_is_tight");

            bool lep_is_tight = (object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTMU || 
                    object_selector_.event_type() == TTObjectSelector::EvtType::TIGHTEL);

            //tracker_.track("b4 obj sel", leptype);
            if(lep_is_tight) {
                tracker_.group(leptype);
                tracker_.track("object selection", leptype);
            }

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
            presel_dir << leptype << "/" << njets_type << "/";
            presel_dir << sys_name << "/preselection";
            if(!sync_ && lep_is_tight ){
                fill_presel_plots(presel_dir.str(), event);
                //tracker_.track("not sync and tight lep");
            }

            // btag regions
            auto &clean_jets = object_selector_.clean_jets();
            if( IDJet::id_type(cut_tight_b_) == IDJet::IDType::CSV){
                sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->csvIncl() > B->csvIncl());});
            }
            else if(IDJet::id_type(cut_tight_b_) == IDJet::IDType::MVA){
                sort(clean_jets.begin(), clean_jets.end(), [](IDJet* A, IDJet* B){return(A->CombinedMVA() > B->CombinedMVA());});
            }
            else{
                Logger::log().error() << "Don't know how to sort bjets in Permutations!" << endl;
                throw 42;
            }

                // both highest-tagged jets pass cuts
            bool btag_pass = false;
            if( clean_jets[0]->BTagId(cut_tight_b_) && clean_jets[1]->BTagId(cut_loose_b_) ){
                btag_pass = true;
            }

                // neither highest-tagged jets pass cuts
            bool btag_fail = false;
            if( !clean_jets[0]->BTagId(cut_tight_b_) && !clean_jets[1]->BTagId(cut_loose_b_) ){
                btag_fail = true;
            }

            if( !(btag_pass || btag_fail) ) return;
            //if( btag_fail ) cout << "btag low: " << btag_fail << ", njets = " << njets << endl;
            //if( btag_pass ) cout << "btag high: " << btag_pass << ", njets = " << njets << endl;

            //if(!clean_jets[0]->BTagId(cut_tight_b_)) return;
            ////tracker_.track("tight b pass");
            //if(!clean_jets[1]->BTagId(cut_loose_b_)) return;
            ////tracker_.track("loose b pass");
            
            //if(lep_is_tight) tracker_.track("b cuts", leptype);

            permutator_.preselection(
                        object_selector_.clean_jets(), object_selector_.lepton(), 
                        object_selector_.met(), object_selector_.lepton_charge(),
                        //lep_is_tight) ) return;
                        event.rho().value(), lep_is_tight
            );
            //if( !permutator_.preselection(
            //            object_selector_.clean_jets(), object_selector_.lepton(), 
            //            object_selector_.met(), object_selector_.lepton_charge(),
            //            //lep_is_tight) ) return;
            //            event.rho().value(), lep_is_tight
            //            )
            //) return;
            //if(lep_is_tight) tracker_.track("perm preselection");

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
                //return;

                /// Find best permutation
                permutator_.reset_3J();
                IDJet* wj2 = 0;
                Permutation event_perm(clean_jets[2], wj2, clean_jets[0], clean_jets[1], object_selector_.lepton(), object_selector_.met(), object_selector_.lepton_charge());//wja, wjb, bhad, blep
                Permutation best_perm = lost_best_perm( event_perm, "Total", presel_dir.str(), lep_is_tight );
           
                bool reco_success_3J = !best_perm.IsEmpty() && best_perm.Prob() < 1e9;
                if(!sync_ && !reco_success_3J){
                    //tracker_.track("not sync and not reco success");
                    return;
                }
                if(lep_is_tight) tracker_.track("best perm");

                //if(sync_) {
                //    if(lep_is_tight) {
                //        sync_info_.Run  = event.run;
                //        sync_info_.LumiSection = event.lumi;
                //        sync_info_.Event= event.evt;
                //        sync_info_.hasMuon = (object_selector_.lepton_type() == -1) ? 0 : 1;
                //        sync_info_.RecoSuccess = reco_success_3J;
                //        if(reco_success_3J){
                //                // corrected thad
                //            TLorentzVector corrected_thad = alpha_corr_.Alpha_THad(best_perm, "1d", "None", "Mtt"); // 3rd argument set to "None" to not apply correction
                //            //TLorentzVector corrected_thad = alpha_corr_.Alpha_THad(best_perm, "1d", "E", "Mtt"); // 3rd argument set to "None" to not apply correction
                //            sync_info_.MassTT = (corrected_thad + best_perm.TLep()).M();
                //        }
                //        else{
                //            sync_info_.MassTT = -1;
                //        }
                //        sync_tree_->Fill();
                //        tracker_.track("SYNC END", leptype);
                //    }
                //    return;
                //}


                //create event dir (contains the dir path to be filled)
                stringstream evtdir;		
                evtdir << leptype << "/" << njets_type;

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

                //// btag and MT categories
                if( shift == Sys::NOSYS ){
                    if( btag_pass ) evtdir << "/csvPass";
                    if( btag_fail )  evtdir << "/csvFail";
                }
                else{
                    if( !btag_pass ) return;
                    evtdir << "/csvPass";
                }

                //fill right category
                fill_3J_selection_plots(evtdir.str(), event, best_perm, runpdf);
                if(lep_is_tight) tracker_.track("END", leptype);

                return;
            } // njets == 3 requirement


    //// events with 4+ jets

            //if( lep_is_tight && btag_fail ) cout << "1: Njets=" << njets << ", tight iso: " << lep_is_tight << ", btag fail: " << btag_fail << endl;

            //Find best permutation
            Permutation best_permutation;
            //int idx = 0;
            for(auto& test_perm : permutator_.permutations()) {
                solver_.Solve(test_perm);

                if(test_perm.Prob()  < best_permutation.Prob()) {
                    best_permutation = test_perm;
                    //cout << "This is the best so far!" << endl; 
                }
                // cout << endl;
            }

            bool reco_success = (best_permutation.IsComplete() && best_permutation.Prob() <= 1E9);
            if(!sync_ && !reco_success){
                //tracker_.track("not sync and not reco success");
                return;
            }
            if(lep_is_tight) tracker_.track("best perm");


            //if( lep_is_tight && btag_fail ) cout << "2: Njets=" << njets << ", tight iso: " << lep_is_tight << ", btag fail: " << btag_fail << endl;


            //if(sync_) {
            //    if(lep_is_tight) {
            //        sync_info_.Run  = event.run;
            //        sync_info_.LumiSection = event.lumi;
            //        sync_info_.Event= event.evt;
            //        sync_info_.hasMuon = (object_selector_.lepton_type() == -1) ? 0 : 1;
            //        sync_info_.RecoSuccess = reco_success;
            //        if(reco_success)
            //            sync_info_.MassTT = best_permutation.LVect().M();
            //        else
            //            sync_info_.MassTT = -1;
            //        sync_tree_->Fill();
            //        tracker_.track("SYNC END", leptype);
            //    }
            //    return;
            //create event dir (contains the dir path to be filled)
            stringstream evtdir;		
            evtdir << leptype << "/" << njets_type;

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


            //if( lep_is_tight && btag_fail ) cout << "3: Njets=" << njets << ", tight iso: " << lep_is_tight << ", btag fail: " << btag_fail << endl;

            //// btag and MT categories
            if( shift == Sys::NOSYS ){
                if( btag_pass ) evtdir << "/csvPass";
                if( btag_fail )  evtdir << "/csvFail";
                //if( btag_fail ){
                //    evtdir << "/csvFail";
                //    cout << "Njets=" << njets << ", tight iso: " << tight << ", btag fail: " << btag_fail << endl;
                //}
            }
            else{
                if( !btag_pass ) return;
                evtdir << "/csvPass";
            }

            //fill right category
            fill_selection_plots(evtdir.str(), event, best_permutation, runpdf);
            if(lep_is_tight) tracker_.track("END", leptype);

        } // end of process_evt

        //This method is called once every file, contains the event loop
        //run your proper analysis here
        virtual void analyze(){
            opts::variables_map &values = URParser::instance().values();
            //int limit = 10000;
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
            Logger::log().debug() << "htt_qcd_est::analyze" << endl;
            Logger::log().debug() << "-- DONE -- reporting every -- " << report << " out of " << tree_->GetEntries() << endl;
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
                ("nosys", opts::value<int>()->default_value(1), "do not run systematics")
                //("nosys", "do not run systematics")
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
            URDriver<htt_qcd_est> test;
            int excode = test.run();
            //Logger::log().debug() << "RUNNING DONE " << std::endl;
            auto files = gROOT->GetListOfFiles(); //make ROOT aware that some files do not exist, because REASONS
            Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl; //need to print out this otherwise ROOT loses its shit in 7.4.X (such I/O, much features)
            opts::variables_map &values = parser.values();
            string output_file = values["output"].as<std::string>();
            Logger::log().debug() << "  Output File: " << output_file << endl;
            return excode;
        }
