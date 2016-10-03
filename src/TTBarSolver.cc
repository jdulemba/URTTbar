#include "Analyses/URTTbar/interface/Permutation.h"
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "Analyses/URTTbar/interface/NeutrinoSolver.h"
#include "Analyses/URTTbar/interface/IDMet.h"
#include "URAnalysis/AnalysisFW/interface/Wards.h"
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include <limits>
#include "URAnalysis/AnalysisFW/interface/Logger.h"

TTBarSolver::TTBarSolver(bool active) {
	if(!active) {
		Logger::log().warning() << "WARNING! TTBarSolver is being disabled!" << endl;
		return;
	}
  URParser &parser = URParser::instance();
  //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
  parser.addCfgParameter<string>("ttsolver", "filename", "");
  parser.addCfgParameter<string>("ttsolver", "dirname" , "");
  parser.addCfgParameter<bool>("ttsolver", "nusolver" , "");
  parser.addCfgParameter<bool>("ttsolver", "invmass" , "");
  parser.addCfgParameter<bool>("ttsolver", "btag" , "", false);
  parser.addCfgParameter<bool>("ttsolver", "ptratio" , "", false);
  parser.addCfgParameter<bool>("ttsolver", "qarkgluon" , "", false);
  
  parser.parseArguments();
  
  string fname = parser.getCfgPar<string>("ttsolver", "filename");
  string dname = parser.getCfgPar<string>("ttsolver", "dirname" );
  USEBTAG_ = parser.getCfgPar<bool>("ttsolver", "btag");
  USENS_   = parser.getCfgPar<bool>("ttsolver", "nusolver");
  USEMASS_ = parser.getCfgPar<bool>("ttsolver", "invmass" );
	useptratios_  = parser.getCfgPar<bool>("ttsolver", "ptratio" );
	usewjetqgtag_ = parser.getCfgPar<bool>("ttsolver", "qarkgluon" );

  TFile probfile(DataFile(fname).path().c_str());
  TDirectory *dir = (TDirectory*) probfile.Get(dname.c_str());

  if(dir) {
		HistoOwnershipWard hward;
		if(USEMASS_) 
			WTmass_right_ = preproccess_histo<TH2D>(dir, "mWhad_vs_mtophad_right", "wt_right");
		if(USENS_)
			N_right_ = preproccess_histo<TH1D>(dir, "nusolver_chi2_right", "nu_right");
		if(USEBTAG_) {
			wj1_btag_right_ = preproccess_histo<TH1D>(dir, "wjets_bcMVA_p11_right", "wjets_bcMVA_p11_wrong", "j1btag_right"); //best wjet cMVA^11
			wj2_btag_right_ = preproccess_histo<TH1D>(dir, "wjets_wcMVA_p11_right", "wjets_wcMVA_p11_wrong", "j2btag_right"); //worst wjet cMVA^11
		}
		if(useptratios_) {
			lep_b_ratio_right_ = preproccess_histo<TH1D>(dir, "lb_ratio_right" ,  "lb_ratio_wrong", "lbr_right");
			wj2_b_ratio_right_ = preproccess_histo<TH1D>(dir, "w2b_ratio_right", "w2b_ratio_wrong", "jbr_right");
		}
		if(usewjetqgtag_) {
			wj1_qgtag_right_ = preproccess_histo<TH1D>(dir, "wjets_bqgt_right", "wjets_bqgt_wrong", "j1qgtag_right");
			wj2_qgtag_right_ = preproccess_histo<TH1D>(dir, "wjets_wqgt_right", "wjets_wqgt_wrong", "j2qgtag_right");			
		}
  } else {
		Logger::log().error() << "could not find directory: "<< dname << " in " << fname << endl;
		throw 42;
	}
}


TTBarSolver::~TTBarSolver()
{}

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
	double ptratios = numeric_limits<double>::max();
	double qgtest   = numeric_limits<double>::max();

	NeutrinoSolver NS(hyp.L(), hyp.BLep(), mw_, mtop_);
	hyp.Nu(NS.GetBest(hyp.MET()->Px(), hyp.MET()->Py(), 1, 1, 0., nschi)); //ignore MET covariance matrix, take bare distance

	//Fill chi discriminant
	if(nschi > 0. && nschi < 10000. && N_right_)
		nstest = -1.*Log(N_right_->Interpolate(Sqrt(nschi)));

	double mwhad = hyp.WHad().M();
	double mthad = hyp.THad().M();
	auto min_max = [] (const double& a, const double& b) -> pair<const double,const double> { 
		return (b<a) ? std::make_pair(b, a) : std::make_pair(a, b); };

	if(mthad < 500. && mwhad < 500. && WTmass_right_) {
		double massdisval = WTmass_right_->Interpolate(mwhad, mthad);
		if(massdisval > 1.0E-10) {masstest = -1.*Log(massdisval);}
		//masstest = -1.*Log(WTmass_right_->Interpolate(mwhad, mthad)/Max(1., WTmass_wrong_->Interpolate(mwhad, mthad)));
	}
	if(USEBTAG_) {
		auto btags = min_max(hyp.WJa()->CombinedMVA(), hyp.WJb()->CombinedMVA());
		btagtest = -1.*Log(wj1_btag_right_->Interpolate(pow(btags.second, 11)));
		btagtest -= 1.*Log(wj2_btag_right_->Interpolate(pow(btags.first , 11)));
	}
	if(useptratios_) {
		auto jpts = min_max(hyp.WJa()->Pt(), hyp.WJb()->Pt());
		ptratios = -1.*Log(lep_b_ratio_right_->Interpolate(hyp.L()->Pt()/hyp.BLep()->Pt()));
		ptratios -= 1.*Log(wj2_b_ratio_right_->Interpolate(jpts.first/hyp.BHad()->Pt()));
	}
	if(usewjetqgtag_) {
		auto qgtag = min_max(hyp.WJa()->qgTag(), hyp.WJb()->qgTag());
		qgtest = -1.*Log(wj1_qgtag_right_->Interpolate(qgtag.second));
		qgtest -= 1.*Log(wj2_qgtag_right_->Interpolate(qgtag.first ));
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
	if(useptratios_) {res += ptratios;}
	if(usewjetqgtag_) {res += qgtest;}

	//fix values in permutation
	hyp.Prob(res);
	hyp.NuChisq(nschi);
	hyp.NuDiscr(nstest);
	hyp.BDiscr(btagtest);
	hyp.MassDiscr(masstest);
	hyp.QGDiscr(qgtest);
	hyp.JRatioDiscr(ptratios);
}

