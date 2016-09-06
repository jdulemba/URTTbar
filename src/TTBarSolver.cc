#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "Analyses/URTTbar/interface/NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/IDMet.h"
#include "URAnalysis/AnalysisFW/interface/Wards.h"
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"

TTBarSolver* TTBarSolver::TTBS = 0; 

void myfuncln(Int_t& npar, Double_t* deriv, Double_t& val, Double_t* par, Int_t flag)
{
	val = TTBarSolver::TTBS->Test(par);
}

TTBarSolver::TTBarSolver() : 
	minuit_(9), 
	WTmass_right_(),
	BTag_right_(),
	N_right_(), 
	WTmass_wrong_(),
	BTag_wrong_(), 
	N_wrong_()
{
}

TTBarSolver::TTBarSolver(bool dummy): 
  minuit_(9) {
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
    WTmass_right_ = preproccess_histo<TH2D>(dir, (wtname+"_right"), "wt_right");
    WTmass_wrong_ = preproccess_histo<TH2D>(dir, (wtname+"_wrong"), "wt_wrong");
    BTag_right_ = preproccess_histo<TH1D>(dir, (bname+"_right"), "b_right");
    BTag_wrong_ = preproccess_histo<TH1D>(dir, (bname+"_wrong"), "b_wrong");
    N_right_ = preproccess_histo<TH1D>(dir, (nuname+"_right"), "nu_right");
    N_wrong_ = preproccess_histo<TH1D>(dir, (nuname+"_wrong"), "nu_wrong");
  }  
}

void TTBarSolver::Solve(Jet* bhad, Jet* j1had, Jet* j2had, Jet* blep, TLorentzVector* llep, IDMet* met)
{
	bhad_ = bhad;
	j1had_ = j1had;
	j2had_ = j2had;
	blep_ = blep;
	llep_ = llep;
	met_ = met;
	bool SMEAR = false;
	bool LEPMASS = false;
	ubhad_ = 0.05;
	uj1had_ = 0.05;
	uj2had_ = 0.05;
	ublep_ = 0.05;
	ullep_ = 0.01;
	umetx_ = 1.;//met->pxunctot(); 
	umety_ = 1.;//met->pyunctot();
	//umety_ = 0.05*met->Py();//Sqrt(met_->pyUnc());
	//umetx_ = 25;//Sqrt(met_->pxUnc());
	//umety_ = 25;//Sqrt(met_->pyUnc());
	rhomet_ = 0.;//met_->pxpyUnc()/(umetx_*umety_);

	nschi_ = -1;
	res_ = 1.E10;
	nstest_ = 1.E10;
	masstest_ = 1.E10;
	//cout << bhad->csvIncl() << " " << blep->csvIncl() << " " << j1had->csvIncl() << " " << j2had->csvIncl() << endl;
  if(USEBTAG_ && BTag_right_ && BTag_wrong_) {
    //std::cout << bhad->csvIncl() << std::endl;
    btagtest_ = -1.*Log(BTag_right_->Interpolate(bhad->csvIncl())/BTag_wrong_->Interpolate(bhad->csvIncl()));
    btagtest_ -= Log(BTag_right_->Interpolate(blep->csvIncl())/BTag_wrong_->Interpolate(blep->csvIncl()));
    btagtest_ -= Log(BTag_wrong_->Interpolate(j1had->csvIncl())/BTag_right_->Interpolate(j1had->csvIncl()));
    btagtest_ -= Log(BTag_wrong_->Interpolate(j2had->csvIncl())/BTag_right_->Interpolate(j2had->csvIncl()));
	}

	TTBS = this;
	minuit_.SetFCN(myfuncln);
	minuit_.SetPrintLevel(-1);
	Int_t flags = 0;
	minuit_.mnparm(0, "top mass", 173., 1., 10, 1000, flags);
	minuit_.mnparm(1, "w mass", 80., 1., 10, 1000, flags);
	minuit_.mnparm(2, "ubhad", 1, 0.01, 0.5, 1.5, flags);
	minuit_.mnparm(3, "uj1had", 1, 0.01, 0.5, 1.5, flags);
	minuit_.mnparm(4, "uj2had", 1, 0.01, 0.5, 1.5, flags);
	minuit_.mnparm(5, "ublep", 1, 0.01, 0.5, 1.5, flags);
	minuit_.mnparm(6, "ullep", 1, 0.01, 0.5, 1.5, flags);
	minuit_.mnparm(7, "umetx", 1, 0.01, 0.5, 1.5, flags);
	minuit_.mnparm(8, "umety", 1, 0.01, 0.5, 1.5, flags);
	//minuit_.mnparm(9, "metz", neutrino.Pz(), 1., 0., 10000., flags);
	//minuit_.mnparm(10, "wn mass", 80., 1., 10, 1000, flags);
	if(!LEPMASS)
	{
		minuit_.FixParameter(0);
		minuit_.FixParameter(1);
	}
	if(!SMEAR)
	{
		minuit_.FixParameter(2);
		minuit_.FixParameter(3);
		minuit_.FixParameter(4);
		minuit_.FixParameter(5);
		minuit_.FixParameter(6);
		minuit_.FixParameter(7);
		minuit_.FixParameter(8);
	}
	if(!LEPMASS && !SMEAR)
	{
		double par[9];
		par[0] = 173.;
		par[1] = 80.;
		par[2] = 1.;
		par[3] = 1.;
		par[4] = 1.;
		par[5] = 1.;
		par[6] = 1.;
		par[7] = 1.;
		par[8] = 1.;
		Test(par);
	}
	else
	{
		minuit_.SetMaxIterations(500);
		minuit_.Migrad();
	}
}

double TTBarSolver::Test(double* par)
{
	bhadT_ = TLorentzVector(bhad_->Px()*par[2], bhad_->Py()*par[2], bhad_->Pz()*par[2], bhad_->E()*par[2]);
	j1hadT_ = TLorentzVector(j1had_->Px()*par[3], j1had_->Py()*par[3], j1had_->Pz()*par[3], j1had_->E()*par[3]);
	j2hadT_ = TLorentzVector(j2had_->Px()*par[4], j2had_->Py()*par[4], j2had_->Pz()*par[4], j2had_->E()*par[4]);
	blepT_ = TLorentzVector(blep_->Px()*par[5], blep_->Py()*par[5], blep_->Pz()*par[5], blep_->E()*par[5]);
	llepT_ = TLorentzVector(llep_->Px()*par[6], llep_->Py()*par[6], llep_->Pz()*par[6], llep_->E()*par[6]);
	NeutrinoSolver NS(&llepT_, &blepT_, par[1], par[0]);
	metT_ = TLorentzVector(NS.GetBest(met_->Px()*par[7], met_->Py()*par[8], umetx_, umety_, rhomet_, nschi_));
	//cout << nschi_ << " NS " << (metT_ + *llep_ + *blep_).M() << " " << (metT_ + *llep_).M() << endl;
	if(nschi_ > 0. && nschi_ < 10000. && N_right_)
	{
		nstest_ = -1.*Log(N_right_->Interpolate(Sqrt(nschi_)));
	}

	double mwhad = (j1hadT_ + j2hadT_).M();
	double mthad = (j1hadT_ + j2hadT_ + bhadT_).M();
	if(mthad < 500. && mwhad < 500. && WTmass_right_)
	{
		double massdisval = WTmass_right_->Interpolate(mwhad, mthad);
		if(massdisval > 1.0E-10) {masstest_ = -1.*Log(massdisval);}
		//masstest_ = -1.*Log(WTmass_right_->Interpolate(mwhad, mthad)/Max(1., WTmass_wrong_->Interpolate(mwhad, mthad)));
	}

	res_ = 0.;
	double sqrt2 = Sqrt(2);
	res_ += Power((par[0]-173)/(20), 2);
	res_ += Power((par[1]-80)/(20), 2);
	res_ += Power((par[2]-1.)/ubhad_/sqrt2 , 2);
	res_ += Power((par[3]-1.)/uj1had_/sqrt2 , 2);
	res_ += Power((par[4]-1.)/uj2had_/sqrt2 , 2);
	res_ += Power((par[5]-1.)/ublep_/sqrt2 , 2);
	res_ += Power((par[6]-1.)/ullep_/sqrt2 , 2);
	res_ += Power((par[7]-1.)/umetx_/sqrt2 , 2);
	res_ += Power((par[8]-1.)/umety_/sqrt2 , 2);
	if(USEMASS_) {res_ += masstest_;}
	if(USENS_) {res_ += nstest_;}
	if(USEBTAG_) {res_ += btagtest_;}

	return(res_);
}

