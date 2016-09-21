#include "Analyses/URTTbar/interface/Permutation.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "Analyses/URTTbar/interface/NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/IDMet.h"
#include "URAnalysis/AnalysisFW/interface/Wards.h"
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include <limits>

TTBarSolver::TTBarSolver() : 
	WTmass_right_(),
	N_right_()
{
}

TTBarSolver::TTBarSolver(bool dummy) {
  URParser &parser = URParser::instance();
  //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
  parser.addCfgParameter<string>("ttsolver", "filename", "");
  parser.addCfgParameter<string>("ttsolver", "dirname" , "");
  parser.addCfgParameter<bool>("ttsolver", "btag" , "");
  parser.addCfgParameter<bool>("ttsolver", "nusolver" , "");
  parser.addCfgParameter<bool>("ttsolver", "invmass" , "");
  
  parser.parseArguments();
  
  string fname = parser.getCfgPar<string>("ttsolver", "filename");
  string dname = parser.getCfgPar<string>("ttsolver", "dirname" );
  bool btag    = parser.getCfgPar<bool>("ttsolver", "btag");
  bool nusl    = parser.getCfgPar<bool>("ttsolver", "nusolver");
  bool mass    = parser.getCfgPar<bool>("ttsolver", "invmass" );

  TFile probfile(DataFile(fname).path().c_str());
  TDirectory *td = (TDirectory*) probfile.Get(dname.c_str());
  Init(td, btag, nusl, mass);
}


TTBarSolver::~TTBarSolver()
{}

void TTBarSolver::Init(string filename, bool usebtag, bool usens, bool usemass)
{
	USEBTAG_ = usebtag;
	USENS_ = usens;
	USEMASS_ = usemass;
  if(!filename.empty()) {
		DirectoryWard dward;
    TFile* probfile = TFile::Open(filename.c_str(), "READ");
		Init(probfile, usebtag, usens, usemass);
  }  
}

void TTBarSolver::Init(TDirectory* dir, bool usebtag, bool usens, bool usemass, string wtname, string bname, string nuname) {
	USEBTAG_ = usebtag;
	USENS_ = usens;
	USEMASS_ = usemass;
  if(dir) {
		HistoOwnershipWard hward;
		if(USEMASS_) {
			WTmass_right_ = preproccess_histo<TH2D>(dir, (wtname+"_right"), "wt_right");
		}
		// if(USEBTAG_) {
		// 	BTag_right_ = preproccess_histo<TH1D>(dir, (bname+"_right"), "b_right");
		// }
		if(USENS_) {
			N_right_ = preproccess_histo<TH1D>(dir, (nuname+"_right"), "nu_right");
		}
  }  
}

void TTBarSolver::Solve(Permutation &hyp, bool lazy)
{
  if(!lazy && !hyp.IsComplete()) {                          
    Logger::log().fatal() << "The permutation you are trying to solve is not complete!" << std::endl;
    throw 42;
  }

	double nschi    = numeric_limits<double>::max();
	double res      = numeric_limits<double>::max();
	double nstest   = numeric_limits<double>::max();
	double masstest = numeric_limits<double>::max();
	double btagtest = numeric_limits<double>::max();

	NeutrinoSolver NS(hyp.L(), hyp.BLep(), mw_, mtop_);
	hyp.Nu(NS.GetBest(hyp.MET()->Px(), hyp.MET()->Py(), 1, 1, 0., nschi)); //ignore MET covariance matrix, take bare distance

	//Fill chi discriminant
	if(nschi > 0. && nschi < 10000. && N_right_)
		nstest = -1.*Log(N_right_->Interpolate(Sqrt(nschi)));

	double mwhad = hyp.WHad().M();
	double mthad = hyp.THad().M();
	if(mthad < 500. && mwhad < 500. && WTmass_right_) {
		double massdisval = WTmass_right_->Interpolate(mwhad, mthad);
		if(massdisval > 1.0E-10) {masstest = -1.*Log(massdisval);}
		//masstest = -1.*Log(WTmass_right_->Interpolate(mwhad, mthad)/Max(1., WTmass_wrong_->Interpolate(mwhad, mthad)));
	}

  // if(USEBTAG_ && BTag_right_) {
  //   //std::cout << bhad->csvIncl() << std::endl;
  //   btagtest_  = -1.*Log(BTag_right_->Interpolate(bhad->csvIncl())/BTag_wrong_->Interpolate(bhad->csvIncl()));
  //   btagtest_ -= Log(BTag_right_->Interpolate(blep->csvIncl())/BTag_wrong_->Interpolate(blep->csvIncl()));
  //   btagtest_ -= Log(BTag_wrong_->Interpolate(j1had->csvIncl())/BTag_right_->Interpolate(j1had->csvIncl()));
  //   btagtest_ -= Log(BTag_wrong_->Interpolate(j2had->csvIncl())/BTag_right_->Interpolate(j2had->csvIncl()));
	// }

	res = 0.;
	if(USEMASS_) {res += masstest;}
	if(USENS_  ) {res += nstest;}
	if(USEBTAG_) {res += btagtest;}

	//fix values in permutation
	hyp.Prob(res);
	hyp.NuChisq(nschi);
	hyp.NuDiscr(nstest);
	hyp.BDiscr(btagtest);
	hyp.MassDiscr(masstest);
}

