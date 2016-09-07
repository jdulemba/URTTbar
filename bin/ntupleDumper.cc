#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include <boost/algorithm/string/predicate.hpp>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/TTBarResponse.h"
#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "Analyses/URTTbar/interface/TTBarPlots.h"
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
#include "URAnalysis/AnalysisFW/interface/EventList.h"
#include "TROOT.h"
//#include <map>

using namespace std;

class ntupleDumper : public AnalyzerBase
{
public:
	enum TTNaming {RIGHT, RIGHT_THAD, RIGHT_TLEP, WRONG, OTHER, NOTSET};
private:
  float evt_weight_;
	bool mc_;

public:
  ntupleDumper(const std::string output_filename):
    AnalyzerBase("ntupleDumper", output_filename),
    evt_weight_(1.) {
	};
  
  ~ntupleDumper() {
  }

  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin() {
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
		bool mc = values["mc"].as<bool>();
		string pick = values["pick"].as<string>();
		EventList picker;
		if(pick.size() != 0) {
			EventList nn(pick);
			picker = nn;
		}

    Logger::log().debug() << "-- DONE -- reporting every -- " << report << endl;
    while(event.next()) {
			if(picker.active()) {
				if(picker.contains(event.run, event.lumi, event.evt)) {
					Logger::log().debug() << "Picking event " << " run: " << event.run << " lumisection: " << 
						event.lumi << " eventnumber: " << event.evt << endl;
				}
				else continue;
			}
			if(limit > 0 && evt_idx > limit) {
				return;
			}
			evt_idx++;
			if(skip > 0 && evt_idx < skip) {
				continue;
			}
			cout << "***** Event "<<event.run<<":"<<event.lumi<<":"<<event.evt<<" *****"<< endl << endl;
			cout << "Trigger: " << endl;
			cout << "  HLT_Ele32_eta2p1_WPTight_Gsf: " << event.trigger().HLT_Ele32_eta2p1_WPTight_Gsf() << endl; 
			cout << "  HLT_IsoMu24: " << event.trigger().HLT_IsoMu24() << endl; 
			cout << "  HLT_IsoTkMu24: " <<  event.trigger().HLT_IsoTkMu24() << endl;
			cout << endl;

			cout << "Event Filters: " << endl;
			cout << "  Flag_HBHENoiseIsoFilter: " <<                  event.filter().Flag_HBHENoiseIsoFilter() << endl; 								 
			cout << "  Flag_EcalDeadCellTriggerPrimitiveFilter: " << 	event.filter().Flag_EcalDeadCellTriggerPrimitiveFilter() << endl;	 
			cout << "  Flag_goodVertices: " << 												event.filter().Flag_goodVertices() << endl;												 
			cout << "  Flag_eeBadScFilter: " << 											event.filter().Flag_eeBadScFilter() << endl;											 
			cout << "  Flag_globalTightHalo2016Filter: " << 					event.filter().Flag_globalTightHalo2016Filter() << endl;					 
			cout << "  Flag_BadPFMuon: " << 													event.filter().Flag_BadPFMuon() << endl;													 
			cout << "  Flag_BadChargedCandidate: " << 	              event.filter().Flag_BadChargedCandidate() << endl;	                
			cout << endl;
			
			cout << "Electrons:" << endl;
			auto& els = event.electrons();
			double rho = event.rho().value();
			for(size_t i=0; i<els.size(); i++) {
				IDElectron el(els[i], rho);
				cout <<" #" << i+1 << endl;
				cout <<"  Momentum (pt, eta, phi): " << el.Pt() << ", " << el.Eta() << ", " << el.Phi() << endl;
				cout <<"  Isolation: "<< el.PFIsolationRho2015()/el.Pt() <<", eta_SC: "<< el.etaSC() << endl;
				cout <<"  Pass cut-based ID (veto to tight): " << el.eidCutVeto() <<", " << el.LooseID25ns()<<", "<< el.MediumID25ns() <<", " <<el.TightID25ns()<<endl;
				cout <<"  Triggering MVA ID (80%, 90%): " << "--"<< endl;
				cout <<"  Pass trigger-emulating preselection: " << endl;
			}
			cout << endl;
			cout << "Muons:" << endl;
			auto& mus = event.muons();
			for(size_t i=0; i<mus.size(); i++) {
				IDMuon mu(mus[i]);
				cout <<" #" << i+1 << endl;
				cout <<"  Momentum (pt, eta, phi): " << mu.Pt() << ", " << mu.Eta() << ", " << mu.Phi() << endl;
				cout <<"  Isolation: "<< mu.RelPFIsoDb() << endl;
				cout <<"  Pass ID (loose to tight): " << mu.isLoose() << ", "<< "-" <<", " << mu.isTight() << endl;
			}
			cout << endl;
			cout << "Jets:" << endl;
			auto& jets = event.jets();
			for(size_t i=0; i<jets.size(); i++) {
				IDJet jet(jets[i]);
				cout <<" #" << i+1 << endl;
				cout <<"  Raw momentum (pt, eta, phi, m): " << jet.uncorrPt() << ", " << jet.uncorrEta() << ", " << jet.uncorrPhi() <<", "<<jet.uncorrM() << endl;
				double pt = (mc) ? jet.Pt()*jet.JER() : jet.Pt();
				cout <<"  Fully corrected pt: " << pt << " JER: " << jet.JER() << " JES: " << jet.Pt()/jet.uncorrPt() << endl;
				if(mc) {
					double jerup = Abs(jet.JERUp()-jet.JER());
					double jerdw = Abs(jet.JERDown()-jet.JER());
					double max_unc = (jerup > jerdw) ? jerup : jerdw;
					cout <<"  JEC uncertainty: "<< jet.JESUnc() << ", JER uncertainty: " << max_unc << endl;
				}
				cout <<"  Pass jet ID: " << jet.ID() << endl;
				cout <<"  CSV: "<< jet.csvIncl()<<", cMVA: "<< jet.CombinedMVA() << endl;
				cout <<"  cCvsL: "<< jet.CvsLtag()<<", cCvsB: "<< jet.CvsBtag() << endl;
			}
			cout << endl;
			cout << "MET:" << endl;
			auto& mets = event.METs();
			IDMet met(mets[0]);
			cout << " Raw MET (pt, phi): --, --" << endl;
			cout << " Corrected MET (pt, phi): " << met.Pt() <<", "<< met.Phi() << endl;
			cout << " Smeared MET pt: " << sqrt( pow(met.pxsmear(), 2) + pow(met.pysmear(), 2) ) << endl;
			cout << endl << endl;
    }  //while(event.next())
  }

  //this method is called at the end of the job, by default saves
  //every histogram/tree produced, override it if you need something more
  virtual void end(){
		outFile_.Write();
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
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("report", opts::value<int>()->default_value(10000), "report every in debug mode")
			("mc", opts::value<bool>()->default_value(false), "report every in debug mode")
      ("pick", opts::value<string>()->default_value(""), "pick from evtlist");
  }
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<ntupleDumper> test;
	int excode = test.run();
	//Logger::log().debug() << "RUNNING DONE " << std::endl;
	auto files = gROOT->GetListOfFiles(); //make ROOT aware that some files do not exist, because REASONS
	Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl; //need to print out this otherwise ROOT loses its shit in 7.4.X (such I/O, much features)
  return excode;
}
