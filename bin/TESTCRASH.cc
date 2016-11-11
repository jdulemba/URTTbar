#include <iostream>
#include "URAnalysis/AnalysisFW/interface/AnalyzerBase.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/URDriver.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "TROOT.h"

#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/TTGenParticleSelector.h"
#include <unordered_map>
#include "URAnalysis/AnalysisFW/interface/RObject.h"

class TESTCRASH : public AnalyzerBase
{
private:
  unsigned long evt_idx_ = 0;
	MCWeightProducer mc_weights_;
  TTGenParticleSelector genp_selector_;
  unordered_map<string, RObject > histos_;
  // Add your private variables/methods here
public:
  TESTCRASH(const std::string output_filename):
    AnalyzerBase("TESTCRASH", output_filename),
		mc_weights_(),
		genp_selector_(TTGenParticleSelector::LHE) {
		opts::variables_map &values = URParser::instance().values();
		string output_file = values["output"].as<std::string>();
		string sample = systematics::get_sample(output_file);
		mc_weights_.init(sample);
	};
  
  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin()
  {
    outFile_.cd();
		histos_["pippo"] = RObject::book<TH1D>("m_tt", "", 400, 0., 2000);			
		histos_["mtop"] = RObject::book<TH1D>("m_top", "", 400, 0., 400);			
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

    URStreamer event(tree_);
    while(event.next())
    {
			if(limit > 0 && evt_idx_ > limit) return;
			evt_idx_++;
			if(skip > 0 && evt_idx_ < skip) continue;
			if(evt_idx_ % report == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << " -- " << event.run<<":"<<event.lumi<<":"<<event.evt <<endl;
			/*

				DO YOUR ANALYSIS HERE!

			 */
			genp_selector_.select(event);
			histos_["pippo"].fill(genp_selector_.ttbar_final_system().M());
			histos_["mtop"].fill(genp_selector_.ttbar_final_system().top.M());
			histos_["mtop"].fill(genp_selector_.ttbar_final_system().tbar.M());
    }
  }

  //this method is called at the end of the job, by default saves
  //every histogram/tree produced, override it if you need something more
  //virtual void end();

  //do you need command-line or cfg options? If so implement this 
  //method to book the options you need. CLI parsing is provided
  //by AnalysisFW/interface/URParser.h and uses boost::program_options
  //look here for a quickstart tutorial: 
  //http://www.boost.org/doc/libs/1_51_0/doc/html/program_options/tutorial.html
  static void setOptions() {
		URParser &parser = URParser::instance();
		opts::options_description &opts = parser.optionGroup("analyzer", "CLI and CFG options that modify the analysis");
		opts.add_options()
      ("limit,l", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("skip,s", opts::value<int>()->default_value(-1), "limit the number of events processed per file")
      ("report,s", opts::value<int>()->default_value(10000), "report every");
	}
};

//make it executable
int main(int argc, char *argv[])
{
  URParser &parser = URParser::instance(argc, argv);
  URDriver<TESTCRASH> test;
  int excode = test.run();
	Logger::log().debug() << "RUNNING DONE " << std::endl;
	auto files = gROOT->GetListOfFiles() ;
	Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl;
	return excode;
}
