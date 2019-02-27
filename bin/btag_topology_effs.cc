#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/TTBarResponse.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"

#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "Analyses/URTTbar/interface/TTGenMatcher.h"
#include "TRandom3.h"
#include "Analyses/URTTbar/interface/helper.h"
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include <unordered_map>

using namespace std;

class btag_topology_effs : public AnalyzerBase
{
    private:
        //histograms and helpers
        CutFlowTracker tracker_;
        //histos
        map<systematics::SysShifts, unordered_map<string, RObject> > histos_;

        //selectors and helpers
        TTObjectSelector object_selector_;
        TTPermutator permutator_;
        TTBarSolver solver_;
        float evt_weight_;
        TRandom3 randomizer_;// = TRandom3(98765);
        MCWeightProducer mc_weights_;

        vector<systematics::SysShifts> systematics_;

        //Scale factors
        LeptonSF electron_sf_, muon_sf_;
        string loose_blabel_, tight_blabel_;

        IDJet::BTag cut_tight_b_=IDJet::BTag::NONE;
        IDJet::BTag cut_loose_b_=IDJet::BTag::NONE;

    public:
        btag_topology_effs(const std::string output_filename):
            AnalyzerBase("btag_topology_effs", output_filename),
            tracker_(),
            loose_blabel_(),
            tight_blabel_(),
            object_selector_(),
            permutator_(),
            solver_(true),
            evt_weight_(1.),
            mc_weights_(),
            electron_sf_("electron_sf", false),
            muon_sf_("muon_sf"){

                //tracker_.verbose(true);    
                //set tracker
                tracker_.use_weight(&evt_weight_);
                object_selector_.set_tracker(&tracker_);
                permutator_.set_tracker(&tracker_);

                opts::variables_map &values = URParser::instance().values();

                //find out which sample are we running on
                string output_file = values["output"].as<std::string>();
                //DataFile solver_input(values["general.ttsolver_input"].as<std::string>());
                string sample = systematics::get_sample(output_file);
                bool isTTbar = boost::starts_with(sample, "ttJets");

                //choose systematics to run based on sample
                systematics_ = {systematics::SysShifts::NOSYS};
                if(!isTTbar) {
                    Logger::log().error() << "This analyzer is only supposed to run on ttbar samples!" << endl;
                    throw 49;
                }

                mc_weights_.init(sample);

                cut_tight_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.tightb"));
                cut_loose_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.looseb"));

                //get cut names
                tight_blabel_ = IDJet::tag2string(cut_tight_b_);
                loose_blabel_ = IDJet::tag2string(cut_loose_b_);
            };

        ~btag_topology_effs() {
        }

        //This method is called once per job at the beginning of the analysis
        //book here your histograms/tree and run every initialization needed
        virtual void begin() {
            outFile_.cd();
            for(auto shift : systematics_) {
                TDirectory* dir_sys = outFile_.mkdir(systematics::shift_to_name.at(shift).c_str());
                string dir = "alljets";
                TDirectory* tdir = dir_sys->mkdir(dir.c_str());
                tdir->cd();

                Logger::log().debug() << "booking: " << dir_sys->GetName() << "/" << dir << std::endl;
                for(string cut : {loose_blabel_, tight_blabel_}) {
                    for(string flav : {"_bjet_", "_cjet_", "_ljet_"}) {
                        for(string cat : {"all", "pass"}) {
                            string name = "btag_" + cut + flav + cat;
                            string key = dir+"/"+name;
                            Logger::log().debug() << "booking: " << name << std::endl;
                            histos_[shift][key]= RObject::book<TH2D>(name.c_str(), "btag SF input histograms;p_{T};#eta", 100, 0, 1000, 60, -3, 3);
                        }
                    }
                }

                dir = "Wjets";
                tdir = dir_sys->mkdir(dir.c_str());
                tdir->cd();
                Logger::log().debug() << "booking: " << dir_sys->GetName() << "/" << dir << std::endl;
                for(auto entry : IDJet::tag_names) {
                    Logger::log().debug() << "Saving flavour probabilities for: " << entry.first << std::endl;
                    for(string flav : {"_bjet_", "_cjet_", "_ljet_"}) {
                        for(string cat : {"all", "pass"}) {
                            string name = "WjetTag_" + entry.first + flav + cat;
                            string key = dir+"/"+name;
                            histos_[shift][key] = RObject::book<TH2D>(name.c_str(), "SF Input;p_{T};#eta", 100, 0, 1000, 60, -3, 3);
                        }
                    }
                }
            }
        }


        void process_evt(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS) { 
            //select reco objects
            if( !object_selector_.select(event, shift) ) return;
            tracker_.track("obj selection");

            //find mc weight
            if(object_selector_.tight_muons().size() == 1)
                evt_weight_ *= muon_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
            if(object_selector_.tight_electrons().size() == 1)
                evt_weight_ *= electron_sf_.get_sf(object_selector_.lepton()->Pt(), object_selector_.lepton()->Eta());
            tracker_.track("MC weights");

            auto &plots = histos_.find(shift)->second;

            //keeping only the n leading jets. 
            std::vector<IDJet*> leading_jets = object_selector_.clean_jets();
            sort(leading_jets.begin(), leading_jets.end(), [](IDJet* A, IDJet* B){return(A->Pt() > B->Pt());});
            int reducedsize = Min(leading_jets.size(), permutator_.njets_max());
            leading_jets.resize(reducedsize);

            bool preselection_pass = permutator_.preselection(
                    object_selector_.clean_jets(), object_selector_.lepton(), object_selector_.met()
                    );
            tracker_.track("permutation pre-selection done (not applied)");

            if( !preselection_pass ) return;
            tracker_.track("perm preselection");

            //Find best permutation
            Permutation best_permutation;
            for(auto test_perm : permutator_.pemutations()) {
                solver_.Solve(test_perm);
                if(test_perm.MassDiscr() < best_permutation.MassDiscr()){
                    best_permutation = test_perm;
                }
            }
            if(!best_permutation.IsComplete() || best_permutation.Prob() > 1E9) return; //FIXME, is right??? best_permutation.Prob() > 1E9

            //fill BTagging SF BEFORE preselection (it already cuts on the BTag value of the jets
            for(auto jet : {best_permutation.BHad(), best_permutation.BLep()}) {
                int jet_flav = Abs(jet->hadronFlavour());
                string flav;
                if(jet_flav == ura::PDGID::b) flav = "_bjet_";
                else if(jet_flav == ura::PDGID::c) flav = "_cjet_";
                else flav = "_ljet_";

                if(jet->BTagId(cut_loose_b_)) 
                    plots["alljets/btag_"+loose_blabel_+flav+"pass"].fill(jet->Pt(), jet->Eta(), evt_weight_);
                if(jet->BTagId(cut_tight_b_))         
                    plots["alljets/btag_"+tight_blabel_+flav+"pass"].fill(jet->Pt(), jet->Eta(), evt_weight_);
                plots["alljets/btag_"+loose_blabel_+flav+"all"].fill(jet->Pt(), jet->Eta(), evt_weight_);
                plots["alljets/btag_"+tight_blabel_+flav+"all"].fill(jet->Pt(), jet->Eta(), evt_weight_);
            } //for(auto jet : leading_jets)
            tracker_.track("BTag plots");

            //NOW, cut on btagging
            if(!(best_permutation.BHad()->BTagId(cut_tight_b_) || best_permutation.BLep()->BTagId(cut_tight_b_))) return;
            if(!(best_permutation.BHad()->BTagId(cut_loose_b_) && best_permutation.BLep()->BTagId(cut_loose_b_))) return;
            tracker_.track("perm btag");

            for(auto entry : IDJet::tag_names) {
                for(auto jet : {best_permutation.WJa(), best_permutation.WJb()}) {
                    int jet_flav = Abs(jet->hadronFlavour());
                    string flav;
                    if(jet_flav == ura::PDGID::b) flav = "_bjet_";
                    else if(jet_flav == ura::PDGID::c) flav = "_cjet_";
                    else flav = "_ljet_";
                    plots["Wjets/WjetTag_"+entry.first+flav+"all"].fill(jet->Pt(), jet->Eta(), evt_weight_);
                    bool passes = boost::contains(entry.first, "CTAG") ? jet->CTagId(entry.second) : jet->BTagId(entry.second);
                    if(passes) {plots["Wjets/WjetTag_"+entry.first+flav+"pass"].fill(jet->Pt(), jet->Eta(), evt_weight_);}
                }
            }
            tracker_.track("end");
        }

        //This method is called once every file, contains the event loop
        //run your proper analysis here
        virtual void analyze()
        {
            unsigned long evt_idx = 0;
            URStreamer event(tree_);

            Logger::log().debug() << "--retrieving running conditions--" << endl;
            opts::variables_map &values = URParser::instance().values();
            int limit = values["limit"].as<int>();
            int skip  = values["skip"].as<int>();
            int report = values["report"].as<int>();
            Logger::log().debug() << "--DONE--" << endl;
            Logger::log().debug() << "-- DONE -- reporting every -- " << report << " out of " << tree_->GetEntries() << endl;
            while(event.next()) {
                if(limit > 0 && evt_idx > limit) {
                    return;
                }
                evt_idx++;
                if(skip > 0 && evt_idx < skip) {
                    continue;
                }
                if(evt_idx % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx << endl;
                tracker_.track("start");

                tracker_.deactivate();    
                for(auto shift : systematics_){
                    evt_weight_ = mc_weights_.evt_weight(event, shift);
                    if(shift == systematics::SysShifts::NOSYS) tracker_.activate();
                    //Logger::log().debug() << "processing: " << shift << endl;
                    process_evt(event, shift);
                    if(shift == systematics::SysShifts::NOSYS) tracker_.deactivate();
                }

            }  //while(event.next())
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
            parser.addCfgParameter<int>("general", "skew_ttpt_distribution", "Should I skew the pt distribution? (0/1)");
            parser.addCfgParameter<int>("general", "pseudotop", "should I use pseudo-top? (0/1)");

            opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
            opts.add_options()
                ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
                ("report", opts::value<int>()->default_value(10000), "report every in debug mode")
                ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file");
            parser.addCfgParameter<std::string>("best_permutation", "tightb", "");
            parser.addCfgParameter<std::string>("best_permutation", "looseb", "");
        }
};

//make it executable
int main(int argc, char *argv[])
{
    URParser &parser = URParser::instance(argc, argv);
    URDriver<btag_topology_effs> test;
    return test.run();
}
