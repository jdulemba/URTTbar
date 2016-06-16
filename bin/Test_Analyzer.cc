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


using namespace TMath;

class Test_Analyzer : public AnalyzerBase
{
private:
  //Counters
  unsigned long evt_idx_ = 0; //event index
  unsigned long up_type_ = 0; //number of up-type jets
  unsigned long down_type_ = 0; //number of donw-type jets
  unsigned long gen_match_ = 0; //number of gen matched events
  unsigned long delta_helframe_costheta_ = 0;
  unsigned long delta_labframe_cosdeltaphi_ = 0;
  unsigned long delta_labframe_deltaphi_ = 0;
  unsigned long delta_helframe_cosdelta_ = 0;
  unsigned long delta_helframe_prodcosth_ = 0;
//  unsigned long delta_E_ = 0; //energy difference between up and down type jets
  unsigned long frac_delta_E_ = 0; //fractional delta_E
 
  //histograms 
  TH1F *njets=0;
  TH1F *lep_pt=0;
  TH1F *lep_charge=0;
  TH1F *up_type_hist=0; //up-type dist for WJa and WJb
  TH1F *down_type_hist=0; //down-type dist for WJa and WJb
  CutFlowTracker tracker_; //tracker_ tracks how many events pass a certain cut
  TH1F *delta_E_hist=0; //energy difference between up and down type jets
  TH1F *frac_delta_E_hist=0; //fractional energy difference between up and down type jets

	//Angular variables
  TH1F *helframe_costheta_lep=0; //angular variable costheta in helicity frame for leptons
  TH1F *helframe_costheta_nu=0; //angular variable costheta in helicity frame for neutrinos
  
  TH1F *helframe_costheta_dtype=0; //angular variable costheta in helicity frame for down type jets
  TH1F *helframe_costheta_utype=0; //angular variable costheta in helicity frame for up type jets
  TH1F *delta_helframe_costheta_hist=0; //difference in the two
  TH1F *frac_delta_helframe_costheta_hist=0; //fractional difference

  TH1F *labframe_cosdeltaphi_lep_dtype=0; //ang var cosdeltaphi in lab frame for down type jets
  TH1F *labframe_cosdeltaphi_lep_utype=0; //ang var cosdeltaphi in lab frame for up type jets
  TH1F *delta_labframe_cosdeltaphi_hist=0; //difference in the two
  TH1F *frac_delta_labframe_cosdeltaphi_hist=0; //fractional difference

  TH1F *labframe_deltaphi_lep_dtype=0; //ang var deltaphi in lab frame for leptonic down type jets
  TH1F *labframe_deltaphi_lep_utype=0; //ang var deltaphi in lab frame for up type jets
  TH1F *delta_labframe_deltaphi_hist=0; //difference in the two
  TH1F * frac_delta_labframe_deltaphi_hist=0;

  TH1F *helframe_cosdelta_lep_dtype=0;//ang var cosdelta in helicity frame for leptonic down type jets
  TH1F *helframe_cosdelta_lep_utype=0; //ang var cosdeltaphi in helicity frame for leptonic up type jets
  TH1F *delta_helframe_cosdelta_hist=0;
  TH1F *frac_delta_helframe_cosdelta_hist=0;

  TH1F *helframe_prodcosth_lep_dtype=0; //ang var prodcosth in helicity frame for leptons
  TH1F *helframe_prodcosth_lep_utype=0; //ang var prodcosth in helicity frame for leptonic up type jets
  TH1F *delta_helframe_prodcosth_hist=0;
  TH1F *frac_delta_helframe_prodcosth_hist=0;

  //selectors and helpers
  TTObjectSelector object_selector_; //selects ttbar objects
  TTGenParticleSelector genp_selector_; //selects generator level objects
  TTGenMatcher matcher_; //matches particles on generator level
  TTPermutator permutator_;

public:
  Test_Analyzer(const std::string output_filename):
    AnalyzerBase("Test_Analyzer", output_filename),
    tracker_(),
    object_selector_(),
    permutator_(),
    matcher_()
    {
      };
  
  //This method is called once per job at the beginning of the analysis
  //book here your histograms/tree and run every initialization needed
  virtual void begin()
  {
    Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;

    outFile_.cd();
   
    //histograms
    njets = new TH1F("nJets", "nJets", 50, 0., 50.);
    lep_pt = new TH1F("lep_pt", "lep_pt", 500, 0., 500.);
    lep_charge = new TH1F("lep_charge", "lep_charge", 5, -2., 2.);
    up_type_hist = new TH1F("up_type_hist", "up_type_hist", 12, -6., 6.);
    down_type_hist = new TH1F("down_type_hist", "down_type_hist", 12, -6., 6.);
    delta_E_hist = new TH1F("delta_E_hist", "delta_E_hist", 500, -800., 800.);
    frac_delta_E_hist = new TH1F("frac_delta_E_hist", "frac_del_E_hist", 500, -2., 2.);

	//angluar variables
    helframe_costheta_lep = new TH1F("helframe_costheta_lep", "", 200, -1., 1.);
    helframe_costheta_nu = new TH1F("helframe_costheta_nu" , "", 200, -1., 1.);

    helframe_costheta_dtype = new TH1F("helframe_costheta_dtype", "", 200, -1., 1.);
    helframe_costheta_utype = new TH1F("helframe_costheta_utype", "", 200, -1., 1.);
    delta_helframe_costheta_hist = new TH1F("delta_helframe_costheta_hist", "", 200, -2.5, 2.5);
    frac_delta_helframe_costheta_hist = new TH1F("frac_delta_helframe_costheta_hist", "", 200, -100., 100.);

    labframe_cosdeltaphi_lep_dtype = new TH1F("labframe_cosdeltaphi_lep_dtype", "", 200, -1., 1.);
    labframe_cosdeltaphi_lep_utype = new TH1F("labframe_cosdeltaphi_lep_utype", "", 200, -1., 1.);
    delta_labframe_cosdeltaphi_hist = new TH1F("delta_labframe_cosdeltaphi_hist", "", 200, -2.5, 2.5);
    frac_delta_labframe_cosdeltaphi_hist = new TH1F("frac_delta_labframe_cosdeltaphi_hist","", 200, -100., 100.);

    labframe_deltaphi_lep_dtype = new TH1F("labframe_deltaphi_lep_dtype", "", 200, 0., Pi());
    labframe_deltaphi_lep_utype = new TH1F("labframe_deltaphi_lep_utype", "", 200, 0., Pi());
    delta_labframe_deltaphi_hist = new TH1F("delta_labframe_deltaphi_hist", "", 200, -8.0, 8.0);
    frac_delta_labframe_deltaphi_hist = new TH1F("frac_delta_labframe_deltaphi_hist", "", 200, -100., 100.);

    helframe_cosdelta_lep_dtype = new TH1F("helframe_cosdelta_lep_dtype", "", 200, -1., 1.);
    helframe_cosdelta_lep_utype = new TH1F("helframe_cosdelta_lep_utype", "", 200, -1., 1.);
    delta_helframe_cosdelta_hist = new TH1F("delta_helframe_cosdelta_hist", "", 200, -2.5, 2.5);
    frac_delta_helframe_cosdelta_hist = new TH1F("frac_delta_helframe_cosdelta_hist", "", 200, -100., 100.);

    helframe_prodcosth_lep_dtype = new TH1F("helframe_prodcosth_lep_dtype", "", 200, -1., 1.);
    helframe_prodcosth_lep_utype = new TH1F("helframe_prodcosth_lep_utype", "", 200, -1., 1.);
    delta_helframe_prodcosth_hist = new TH1F("delta_helframe_prodcosth_hist", "", 200, -2.5, 2.5);
    frac_delta_helframe_prodcosth_hist = new TH1F("frac_delta_helframe_prodcosth_hist", "", 200, -100., 100.);

    Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
  }

  //This method is called once every file, contains the event loop
  //run your proper analysis here
  virtual void analyze()
  {
    Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

    opts::variables_map &values = URParser::instance().values();
    int limit = values["limit"].as<int>();
    int skip  = values["skip"].as<int>();
    int report = values["report"].as<int>();
    if(evt_idx_ >= limit) return;


    //Angular variables (helicity and lab frames from TOP 14 023)
    double cth_l = 0; //costheta_lep
    double cth_n = 0; //costheta_nu
    double cth_d = 0; //costheta_dtype
    double cth_u = 0; //costheta_utype

//    double delta_helframe_costheta = 0;
//    double frac_delta_helframe_costheta = 0;

//    double delta_labframe_cosdeltaphi = 0;
//    double frac_delta_labframe_cosdeltaphi = 0;

//    double delta_labframe_deltaphi = 0;
//    double frac_delta_labframe_deltaphi = 0;

//    double delta_helframe_cosdelta = 0;
//    double frac_delta_helframe_cosdelta = 0;

//    double delta_helframe_prodcosth = 0;
//    double frac_delta_helframe_prodcosth = 0;

    URStreamer event(tree_);
    
    while(event.next()/* && evt_idx_ < 100000*/)
      {
	if(limit > 0 && evt_idx_ > limit) return;
	evt_idx_++;
	if(skip > 0 && evt_idx_ < skip) continue;
	if(evt_idx_ % report == 1) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

	//generator selection
	bool selection = genp_selector_.select(event);
        tracker_.track("gen selection");
        if( !selection ) continue;
	
	//full ttbar event selection trial (reco objects)
	if( !object_selector_.select(event) ) continue;
	tracker_.track("obj selection");
		
	//Generator level matching (Gen matching)
	Permutation matched_perm;
	GenTTBar &ttbar = genp_selector_.ttbar_system();
	if( ttbar.type == GenTTBar::DecayType::SEMILEP){ //gets correct kind of events
		matched_perm = matcher_.match( //needs 4 arguments
				genp_selector_.ttbar_final_system(),
				object_selector_.clean_jets(),
				object_selector_.loose_electrons(),
				object_selector_.loose_muons()
				);
	}
	matched_perm.SetMET(object_selector_.met());
	if( !matched_perm.IsComplete() ) continue;	
	if( !matched_perm.WJa()->match() || !matched_perm.WJb()->match() ){
		Logger::log().debug() << "matched_perm has no match " << endl;
		continue;
	}
	tracker_.track("gen matching");
	++gen_match_;

	//make sure one event is even and the other is odd
	if( (Abs(matched_perm.WJa()->match()->pdgId()) + Abs(matched_perm.WJb()->match()->pdgId()) % 2 == 0) ){
		Logger::log().debug() << "one jet isn't even and the other isn't odd " << endl;
		continue;
	}

	//select up and down type jets
	const IDJet* up_Wjet = ( Abs(matched_perm.WJa()->match()->pdgId()) % 2 == 0 ) ? matched_perm.WJa() : matched_perm.WJb();// up type
	tracker_.track("Up-type WJets");
	up_type_hist->Fill(up_Wjet->match()->pdgId());
	++up_type_;
	
	const IDJet* down_Wjet = ( Abs(matched_perm.WJa()->match()->pdgId()) % 2 == 1 ) ? matched_perm.WJa() : matched_perm.WJb();// down type
	tracker_.track("Down-type WJets");
	down_type_hist->Fill(down_Wjet->match()->pdgId());
	++down_type_;

	//energy (and fractional energy) difference between up and down type jets
	if( (up_Wjet->match()->E() - down_Wjet->match()->E())/((up_Wjet->match()->E() + down_Wjet->match()->E())/2) > 0 ){
		++frac_delta_E_;
	}

	delta_E_hist->Fill(up_Wjet->match()->E() - down_Wjet->match()->E());
	frac_delta_E_hist->Fill((up_Wjet->match()->E() - down_Wjet->match()->E())/((up_Wjet->match()->E() + down_Wjet->match()->E())/2));

	//Angular variables
	hyp::TTbar ttbar_ang(matched_perm); //accesses hypothesis namespace and uses info from matched_perm for var. ttbar_ang

	auto leptopcm = ttbar_ang.tlep().to_CM();
	auto hadtopcm = ttbar_ang.thad().to_CM();
	auto ttcm = ttbar_ang.to_CM();
	auto lep = ttbar_ang.tlep().W().l();
	auto nu = ttbar_ang.tlep().W().nu();
	auto up = ttbar_ang.thad().W().up();
	auto dw = ttbar_ang.thad().W().down();

	cth_l = leptopcm.W().l().unit3D().Dot(ttcm.tlep().unit3D());
	cth_n = leptopcm.W().nu().unit3D().Dot(ttcm.tlep().unit3D());
	cth_d = hadtopcm.W().down().unit3D().Dot(ttcm.thad().unit3D());
	cth_u = hadtopcm.W().up().unit3D().Dot(ttcm.thad().unit3D());

	if(cth_l != cth_l || cth_n != cth_n || cth_d != cth_d || cth_u != cth_u) {
		Logger::log().error() << " NANERROR! "<< cth_l << " " << cth_n << " " << cth_d << " " << cth_u << endl;
	}

		//Fill histograms
	helframe_costheta_lep->Fill(cth_l);
	helframe_costheta_nu->Fill(cth_n);

	helframe_costheta_dtype->Fill(cth_d);
	helframe_costheta_utype->Fill(cth_u);
	if( cth_u-cth_d > 0 ){
		++delta_helframe_costheta_;
	}
	delta_helframe_costheta_hist->Fill(cth_u-cth_d); //up-down
	frac_delta_helframe_costheta_hist->Fill((cth_u-cth_d)/((cth_u+cth_d)/2));// (up-down)/((up+down)/2)

	labframe_cosdeltaphi_lep_dtype->Fill(Cos(lep.DeltaPhi(dw)));
	labframe_cosdeltaphi_lep_utype->Fill(Cos(lep.DeltaPhi(up)));
	if( Cos(lep.DeltaPhi(up))-Cos(lep.DeltaPhi(dw)) > 0 ){
		++delta_labframe_cosdeltaphi_;
	}
	delta_labframe_cosdeltaphi_hist->Fill(Cos(lep.DeltaPhi(up))-Cos(lep.DeltaPhi(dw)));
	frac_delta_labframe_cosdeltaphi_hist->Fill((Cos(lep.DeltaPhi(up))-Cos(lep.DeltaPhi(dw)))/((Cos(lep.DeltaPhi(up))+Cos(lep.DeltaPhi(dw)))/2));

	labframe_deltaphi_lep_dtype->Fill(lep.DeltaPhi(dw));
	labframe_deltaphi_lep_utype->Fill(lep.DeltaPhi(up));
	if(lep.DeltaPhi(up)-lep.DeltaPhi(dw) > 0 ){
		++delta_labframe_deltaphi_;
	}
	delta_labframe_deltaphi_hist->Fill(lep.DeltaPhi(up)-lep.DeltaPhi(dw));
	frac_delta_labframe_deltaphi_hist->Fill((lep.DeltaPhi(up)-lep.DeltaPhi(dw))/((lep.DeltaPhi(up)+lep.DeltaPhi(dw))/2));

	helframe_cosdelta_lep_dtype->Fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()));
	helframe_cosdelta_lep_utype->Fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D()));
	if( leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D())-leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()) > 0 ){
		++delta_helframe_cosdelta_;
	}
	delta_helframe_cosdelta_hist->Fill(leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D())-leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()));
	frac_delta_helframe_cosdelta_hist->Fill((leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D())-leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()))/((leptopcm.W().l().unit3D().Dot(hadtopcm.W().up().unit3D())+leptopcm.W().l().unit3D().Dot(hadtopcm.W().down().unit3D()))/2));

	helframe_prodcosth_lep_dtype->Fill(cth_l*cth_d);
	helframe_prodcosth_lep_utype->Fill(cth_l*cth_u);
	if( cth_l*cth_u-cth_l*cth_d > 0 ){
		++delta_helframe_prodcosth_;
	}
	delta_helframe_prodcosth_hist->Fill(cth_l*cth_u-cth_l*cth_d);
	frac_delta_helframe_prodcosth_hist->Fill((cth_l*cth_u-cth_l*cth_d)/((cth_l*cth_u+cth_l*cth_d)/2));

	//Other histograms filled
	njets->Fill(object_selector_.clean_jets().size()); //number of jets
	lep_pt->Fill(object_selector_.lepton()->Pt()); //lepton pt
	lep_charge->Fill(object_selector_.lepton_charge()); //lepton charge
	
      }
    cout << "delta_helframe_costheta > 0: " << delta_helframe_costheta_ << endl;
    cout << "delta_labframe_cosdeltaphi > 0: " << delta_labframe_cosdeltaphi_ << endl;
    cout << "delta_labframe_deltaphi > 0: " << delta_labframe_deltaphi_ << endl;
    cout << "delta_helframe_cosdelta > 0: " << delta_helframe_cosdelta_ << endl;
    cout << "delta_helframe_prodcosth > 0: " << delta_helframe_prodcosth_ << endl;
    cout << "frac_delta_E > 0: " << frac_delta_E_ << endl;
//    cout << "Number of up-type WJets: " << up_type_ << endl;
//    cout << "Number of down-type WJets: " << down_type_ << endl;
    cout << "Number of gen matched events: " << gen_match_ << endl;
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
  static void setOptions() {
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
  URDriver<Test_Analyzer> test;
  int thing = test.run();
  // auto files = gROOT->GetListOfFiles() ;
  // for(int i=0; i<files->GetSize(); i++) {
  //   TNamed *obj = (TNamed*) files->At(i);
  //   Logger::log().debug() << "file " << i << " " << obj->GetName() << std::endl;
  //   TFile *ff = (TFile*) obj;
  //   ff->Close();
  // }
  // gROOT->CloseFiles();
  // files = gROOT->GetListOfFiles();
  // Logger::log().debug() << "Nfiles " << files->GetSize() << std::endl;
  return thing;
}

