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
#include <list>
//#include "JetMETCorrections/Modules/interface/JetResolution.h"


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

        // Joseph added for perm discriminant
  parser.addCfgParameter<bool>("ttsolver", "permdisc", "", false);
        //
//  parser.addCfgParameter<string>("JERC", "JER_SF","");
//  parser.addCfgParameter<string>("JERC", "PT_JER","");

  
  parser.parseArguments();
  
  string fname = parser.getCfgPar<string>("ttsolver", "filename");
  string dname = parser.getCfgPar<string>("ttsolver", "dirname" );
  USEBTAG_ = parser.getCfgPar<bool>("ttsolver", "btag");
  USENS_   = parser.getCfgPar<bool>("ttsolver", "nusolver");
  USEMASS_ = parser.getCfgPar<bool>("ttsolver", "invmass" );
  useptratios_  = parser.getCfgPar<bool>("ttsolver", "ptratio" );
  usewjetqgtag_ = parser.getCfgPar<bool>("ttsolver", "qarkgluon" );

        // Joseph added for perm discriminant
  USEPERM_ = parser.getCfgPar<bool>("ttsolver", "permdisc");
        //

//  DataFile jer_sf_fname_(parser.getCfgPar<string>("JERC","JER_SF"));
//  DataFile resolution_fname_(parser.getCfgPar<string>("JERC","PT_JER"));
//  resolution_ = JME::JetResolution(resolution_fname_.path());
//  jer_sf_ = JME::JetResolutionScaleFactor(jer_sf_fname_.path());


  TFile probfile(DataFile(fname).path().c_str());
  TDirectory *dir = (TDirectory*) probfile.Get(dname.c_str());

  if(dir) {
		HistoOwnershipWard hward;
		if(USEMASS_) 
			WTmass_right_ = preproccess_histo<TH2D>(dir, "mWhad_vs_mtophad_right", "wt_right");
			//Rel_Delt_WTmass_right_ = preproccess_histo<TH2D>(dir, "Rel_Delta_mWhad_vs_Rel_Delta_mtophad_right", "rel_delt_wt_right");
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

            // Joseph added for discriminants
        if( USEPERM_ ){
            Max_Mjet_3J_right_ = preproccess_histo<TH2D>(dir, "3J_mbpjet_vs_maxmjet_correct", "3J_mbpjet_vs_maxmjet_wrong", "");
            Max_Mjet_4J_right_ = preproccess_histo<TH2D>(dir, "4J_mbpjet_vs_maxmjet_correct", "4J_mbpjet_vs_maxmjet_wrong", "");
        }
//		if(USENS_)
//			N_3J_right_ = preproccess_histo<TH1F>(dir, "3J_nusolver_chi2_correct", "3J_nu_correct");
            //


  }
  else {
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

	auto min_max = [] (const double& a, const double& b) -> pair<const double,const double> { 
		return (b<a) ? std::make_pair(b, a) : std::make_pair(a, b); };

	NeutrinoSolver NS(hyp.L(), hyp.BLep(), mw_, mtop_);
	hyp.Nu(NS.GetBest(hyp.MET()->Px(), hyp.MET()->Py(), 1, 1, 0., nschi)); //ignore MET covariance matrix, take bare distance

	//Fill chi discriminant
	if(N_right_) {
		int binx = N_right_->GetXaxis()->FindFixBin(Sqrt(nschi));
		if(binx <= N_right_->GetNbinsX()) 
			nstest = -1.*Log(N_right_->GetBinContent(binx));
	}


    // Fill mass discriminant
//    if(Rel_Delt_WTmass_right_){
//        // get nominal scale factors
//        float WJa_sf = jer_sf_.getScaleFactor({{JME::Binning::JetEta, hyp.WJa()->Eta()}});
//        float WJb_sf = jer_sf_.getScaleFactor({{JME::Binning::JetEta, hyp.WJb()->Eta()}});
//        float BHad_sf = jer_sf_.getScaleFactor({{JME::Binning::JetEta, hyp.BHad()->Eta()}});
//
//        // get jet resolutions
//        float WJa_res = resolution_.getResolution({{JME::Binning::JetPt, hyp.WJa()->Pt()}, {JME::Binning::JetEta, hyp.WJa()->Eta()}, {JME::Binning::Rho, hyp.RhoVal()}});
//        float WJb_res = resolution_.getResolution({{JME::Binning::JetPt, hyp.WJb()->Pt()}, {JME::Binning::JetEta, hyp.WJb()->Eta()}, {JME::Binning::Rho, hyp.RhoVal()}});
//        float BHad_res =resolution_.getResolution({{JME::Binning::JetPt, hyp.BHad()->Pt()}, {JME::Binning::JetEta, hyp.BHad()->Eta()}, {JME::Binning::Rho, hyp.RhoVal()}});
//
//        // get st devs from true masses
//        float M_j1upj2 = (*hyp.WJa()*(WJa_res*WJa_sf+1.)+*hyp.WJb()).M();
//        float M_j2upj1 = (*hyp.WJb()*(WJb_res*WJb_sf+1.)+*hyp.WJa()).M();
//        float sigma2j = sqrt(pow(M_j1upj2-hyp.WHad().M(),2)+pow(M_j2upj1-hyp.WHad().M(),2));
//
//        float M_j1upj2j3 = (*hyp.WJa()*(WJa_res*WJa_sf+1.)+*hyp.WJb()+*hyp.BHad()).M();
//        float M_j2upj1j3 = (*hyp.WJb()*(WJb_res*WJb_sf+1.)+*hyp.WJa()+*hyp.BHad()).M();
//        float M_j3upj1j2 = (*hyp.BHad()*(BHad_res*BHad_sf+1.)+*hyp.WJa()+*hyp.WJb()).M();
//        float sigma3j = sqrt(pow(M_j1upj2j3-hyp.THad().M(),2)+pow(M_j2upj1j3-hyp.THad().M(),2)+pow(M_j3upj1j2-hyp.THad().M(),2));
//
//        float rel_delt_mw = (hyp.WHad().M()-80.385)/sigma2j; //value from http://pdglive.lbl.gov/Particle.action?node=S043&init=0
//        float rel_delt_mt = (hyp.THad().M()-173.1)/sigma3j; //value from http://pdglive.lbl.gov/Particle.action?node=Q007&init=0
//
//        
//        int binx = TMath::Max(1,
//                    TMath::Min(Rel_Delt_WTmass_right_->GetXaxis()->FindFixBin(rel_delt_mw), Rel_Delt_WTmass_right_->GetNbinsX())
//            );
//        int biny = TMath::Max(1,
//                    TMath::Min(Rel_Delt_WTmass_right_->GetYaxis()->FindFixBin(rel_delt_mt), Rel_Delt_WTmass_right_->GetNbinsY())
//            );
//   
//        double massdisval = Rel_Delt_WTmass_right_->GetBinContent(binx,biny);
//
//       if( massdisval > 1.0E-10 ){
//            masstest = -1.*Log(massdisval);
//        }
//
//    }

// old mass disc
	double mwhad = hyp.WHad().M();
	double mthad = hyp.THad().M();
	if(mthad < 500. && mwhad < 500. && WTmass_right_) {
		int binx = TMath::Max(
			1,
			TMath::Min(
				WTmass_right_->GetXaxis()->FindFixBin(mwhad),
				WTmass_right_->GetNbinsX()
				)
			);
		int biny = TMath::Max(
			1,
			TMath::Min(
				WTmass_right_->GetYaxis()->FindFixBin(mthad),
				WTmass_right_->GetNbinsY()
				)
			);
		double massdisval = WTmass_right_->GetBinContent(binx, biny);
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
    Logger::log().fatal() << "USE of QGTag is not supported any longer!" << std::endl;
    throw 42;
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




////// Joseph added for permutation discriminant
void TTBarSolver::Solve_3J(Permutation &hyp, bool lazy)
{
  if( !lazy ) {                          
    Logger::log().fatal() << "Lazy requirement not met!" << std::endl;
    throw 42;
  }

  if( !hyp.BLep() || !hyp.BHad() ){
    Logger::log().fatal() << "Hadronic and leptonic b's not present!" << std::endl;
    throw 42;
  }

    double nschi_3J = numeric_limits<double>::max();
    double nstest_3J = numeric_limits<double>::max();
	double res = numeric_limits<double>::max();
    double permtest_3J = numeric_limits<double>::max();

    double MaxMjet = numeric_limits<double>::max();
    double Mbpjet = numeric_limits<double>::max();


    auto jet_pair = [] (IDJet* a, IDJet* b) -> pair<IDJet*, IDJet*>{ // template for pairs of jets
        return std::make_pair(a,b); };

	auto max_min = [] (const double& a, const double& b) -> pair<const double,const double> { // template for max,min values
		return (b>a) ? std::make_pair(b, a) : std::make_pair(a, b); };


/// calculated for best and matched permutations

    /// Neutrino Solver
    NeutrinoSolver NS(hyp.L(), hyp.BLep(), mw_, mtop_);
    hyp.Nu(NS.GetBest(hyp.MET()->Px(), hyp.MET()->Py(), 1, 1, 0., nschi_3J)); // ignore MET covariance matrix, take bare distance

//    if( nschi_3J == -1 ) cout << "chi_3J: " << nschi_3J << endl;

    if( N_3J_right_ ){
        int binx = N_3J_right_->GetXaxis()->FindFixBin(Sqrt(nschi_3J));
        if( binx <= N_3J_right_->GetNbinsX() ){
            nstest_3J = -1.*Log(N_3J_right_->GetBinContent(binx));

        }
    }

/// only calculated for best permutation, not matched

    if( hyp.WJa() && !hyp.WJb() ){ // only WJa (wj1), BHad (bj1), and BLep (bj2) exist

    /// Permutation Discriminant Solver
            // creates pairs of jets
        auto BHadWJa_pair = jet_pair(hyp.BHad(),hyp.WJa());
        
            //for each valid pair of jets, find max jet mass within a pair and invariant mass of pair to use for likelihood
         MaxMjet = max_min( BHadWJa_pair.first->M(), BHadWJa_pair.second->M() ).first;
         Mbpjet = ( *BHadWJa_pair.first+*BHadWJa_pair.second ).M();
        
         if( Max_Mjet_3J_right_){
             int binx = TMath::Max(1,
                         TMath::Min(Max_Mjet_3J_right_->GetXaxis()->FindFixBin(MaxMjet), Max_Mjet_3J_right_->GetNbinsX())
                 );
             int biny = TMath::Max(1,
                         TMath::Min(Max_Mjet_3J_right_->GetYaxis()->FindFixBin(Mbpjet), Max_Mjet_3J_right_->GetNbinsY())
                 );
        
             double permdisval_3J = Max_Mjet_3J_right_->GetBinContent(binx,biny);

            if( permdisval_3J > 1.0E-10 ){
                 permtest_3J = -1.*Log(permdisval_3J);
             }

         }

    }

    res = 0.;
    if(USEPERM_ ) {res += permtest_3J;}
    if(USENS_ ) {res += nstest_3J;}

    hyp.Prob(res);
    hyp.NuChisq(nschi_3J);
    hyp.NuDiscr(nstest_3J);
    hyp.PermDiscr(permtest_3J);

}



////// Joseph added for permutation discriminant
void TTBarSolver::Solve_4J(Permutation &hyp, bool lazy)
{
  if( !lazy ) {                          
    Logger::log().fatal() << "Lazy requirement not met!" << std::endl;
    throw 42;
  }

  if( !hyp.BLep() || !hyp.BHad() ){
    Logger::log().fatal() << "Hadronic and leptonic b's not present!" << std::endl;
    throw 42;
  }

  if( !(hyp.WJa() && hyp.WJb() && hyp.BHad() && hyp.BLep()) ){
    Logger::log().fatal() << "4 jets not present!" << endl;
    throw 42;
  }



	double res      = numeric_limits<double>::max();
    double permtest_4J = numeric_limits<double>::max();

    double MaxMjet = numeric_limits<double>::max();
    double Mbpjet = numeric_limits<double>::max();

    double lowest_like_4J = 1e10;

    
    auto jet_pair = [] (IDJet* a, IDJet* b) -> pair<IDJet*, IDJet*>{ // template for pairs of jets
        return std::make_pair(a,b); };

	auto max_min = [] (const double& a, const double& b) -> pair<const double,const double> { // template for max,min values
		return (b>a) ? std::make_pair(b, a) : std::make_pair(a, b); };



//    if( hyp.WJa() && !hyp.WJb() ){ // only WJa exists
//            // creates pairs of jets
//        auto BLepWJa_pair = jet_pair(hyp.BLep(),hyp.WJa());
//        auto BHadWJa_pair = jet_pair(hyp.BHad(),hyp.WJa());
//
//            // creates list of valid pairs of jets
//        list< pair<IDJet*, IDJet*> > pair_list;
//        if( BLepWJa_pair.first != BLepWJa_pair.second ) pair_list.push_back(BLepWJa_pair);
//        if( BHadWJa_pair.first != BHadWJa_pair.second ) pair_list.push_back(BHadWJa_pair);
//
//            //for each valid pair of jets, find max jet mass within a pair and invariant mass of pair to use for likelihood
//        for( list< pair<IDJet*, IDJet*> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it ){
//            MaxMjet = max_min( (*it).first->M(), (*it).second->M() ).first;
//            Mbpjet = ( *(*it).first+*(*it).second ).M();
//
//            if( Max_Mjet_4J_right_){
//                int binx = TMath::Max(1,
//                            TMath::Min(Max_Mjet_4J_right_->GetXaxis()->FindFixBin(MaxMjet), Max_Mjet_4J_right_->GetNbinsX())
//                    );
//                int biny = TMath::Max(1,
//                            TMath::Min(Max_Mjet_4J_right_->GetYaxis()->FindFixBin(Mbpjet), Max_Mjet_4J_right_->GetNbinsY())
//                    );
//        
//                double permdisval_4J = Max_Mjet_4J_right_->GetBinContent(binx,biny);
//                if( permdisval_4J > 1.0E-10 ){ permtest_4J = -1.*Log(permdisval_4J); }
//            }
//        
//        	res = 0.;
//        	if(USEPERM_) {res += permtest_4J;}
//        
//            hyp.PermDiscr(permtest_4J);
//
//        }
//    }
//
//    if( !hyp.WJa() && hyp.WJb() ){ // only WJb exists
//            // creates pairs of jets
//        auto BLepWJb_pair = jet_pair(hyp.BLep(),hyp.WJb());
//        auto BHadWJb_pair = jet_pair(hyp.BHad(),hyp.WJb());
//
//            // creates list of valid pairs of jets
//        list< pair<IDJet*, IDJet*> > pair_list;
//        if( BLepWJb_pair.first != BLepWJb_pair.second ) pair_list.push_back(BLepWJb_pair);
//        if( BHadWJb_pair.first != BHadWJb_pair.second ) pair_list.push_back(BHadWJb_pair);
//
//            //for each valid pair of jets, find max jet mass within a pair and invariant mass of pair to use for likelihood
//        for( list< pair<IDJet*, IDJet*> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it ){
//            MaxMjet = max_min( (*it).first->M(), (*it).second->M() ).first;
//            Mbpjet = ( *(*it).first+*(*it).second ).M();
//
//            if( Max_Mjet_4J_right_){
//                int binx = TMath::Max(1,
//                            TMath::Min(Max_Mjet_4J_right_->GetXaxis()->FindFixBin(MaxMjet), Max_Mjet_4J_right_->GetNbinsX())
//                    );
//                int biny = TMath::Max(1,
//                            TMath::Min(Max_Mjet_4J_right_->GetYaxis()->FindFixBin(Mbpjet), Max_Mjet_4J_right_->GetNbinsY())
//                    );
//        
//                double permdisval_4J = Max_Mjet_4J_right_->GetBinContent(binx,biny);
//                if( permdisval_4J > 1.0E-10 ){ permtest_4J = -1.*Log(permdisval_4J); }
//            }
//        
//        	res = 0.;
//        	if(USEPERM_) {res += permtest_4J;}
//        
//            hyp.PermDiscr(permtest_4J);
//
//        }
//    }
    if( hyp.WJa() && hyp.WJb() ){ // both W jets exist
            // creates pairs of jets
        auto BLepWJa_pair = jet_pair(hyp.BLep(),hyp.WJa());
        auto BHadWJa_pair = jet_pair(hyp.BHad(),hyp.WJa());
        auto BLepWJb_pair = jet_pair(hyp.BLep(),hyp.WJb());
        auto BHadWJb_pair = jet_pair(hyp.BHad(),hyp.WJb());
        auto BHadBLep_pair = jet_pair(hyp.BHad(),hyp.BLep());
        auto WJaWJb_pair = jet_pair(hyp.WJa(),hyp.WJb());

            // creates list of valid pairs of jets
        list< pair<IDJet*, IDJet*> > pair_list;
        if( BLepWJa_pair.first != BLepWJa_pair.second ) pair_list.push_back(BLepWJa_pair);
        if( BHadWJa_pair.first != BHadWJa_pair.second ) pair_list.push_back(BHadWJa_pair);
        if( BLepWJb_pair.first != BLepWJb_pair.second ) pair_list.push_back(BLepWJb_pair);
        if( BHadWJb_pair.first != BHadWJb_pair.second ) pair_list.push_back(BHadWJb_pair);
        if( BHadBLep_pair.first != BHadBLep_pair.second ) pair_list.push_back(BHadBLep_pair);
        if( WJaWJb_pair.first != WJaWJb_pair.second ) pair_list.push_back(WJaWJb_pair);


        vector<double> perm_vec;
            //for each valid pair of jets, find max jet mass within a pair and invariant mass of pair to use for likelihood
        for( list< pair<IDJet*, IDJet*> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it ){
            MaxMjet = max_min( (*it).first->M(), (*it).second->M() ).first;
            Mbpjet = ( *(*it).first+*(*it).second ).M();

            if( Max_Mjet_4J_right_){
                int binx = TMath::Max(1,
                            TMath::Min(Max_Mjet_4J_right_->GetXaxis()->FindFixBin(MaxMjet), Max_Mjet_4J_right_->GetNbinsX())
                    );
                int biny = TMath::Max(1,
                            TMath::Min(Max_Mjet_4J_right_->GetYaxis()->FindFixBin(Mbpjet), Max_Mjet_4J_right_->GetNbinsY())
                    );
        
                double permdisval_4J = Max_Mjet_4J_right_->GetBinContent(binx,biny);
                if( permdisval_4J > 1.0E-10 ){
                    permtest_4J = -1.*Log(permdisval_4J);
                    perm_vec.push_back(permtest_4J);
                }
            }
        
        	res = 0.;
        	if(USEPERM_) {res += permtest_4J;}
        }

        if( perm_vec.size() > 0 ){
            for( vector<double>::iterator low_it = perm_vec.begin(); low_it != perm_vec.end(); ++low_it ){
                if( *low_it < lowest_like_4J ){
                    lowest_like_4J = *low_it;
                }
            }
            hyp.PermDiscr(lowest_like_4J);
        }

        hyp.PermDiscr_Vec(perm_vec);
    }
}
