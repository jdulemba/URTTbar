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
#include "Analyses/URTTbar/interface/TTBarSolver.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"
#include "Analyses/URTTbar/interface/GenObject.h"
#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include "URAnalysis/AnalysisFW/interface/RObject.h"
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include <algorithm>
#include <list>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include <boost/algorithm/string/predicate.hpp>
#include "TMath.h"
#include "TRandom3.h"
#include <set>
#include "Analyses/URTTbar/interface/IDMet.h"
#include "TUUID.h"
#include "Analyses/URTTbar/interface/systematics.h"
#include "Analyses/URTTbar/interface/TTBarPlots.h"
#include "Analyses/URTTbar/interface/MCWeightProducer.h"
#include "Analyses/URTTbar/interface/BTagSFProducer.h"
#include "Analyses/URTTbar/interface/LeptonSF.h"
#include "Analyses/URTTbar/interface/PDFuncertainty.h"

#include "TROOT.h"
#include "TTree.h"

using namespace TMath;

class jet_effs : public AnalyzerBase
{
	private:
        //counters
	    unsigned long evt_idx_ = 0; //event index
	    double njets_ = 0;
	    double jets3_ = 0;
	    double matched_BHadWJa_ = 0;
	    double matched_BHadWJb_ = 0;
	    double matched_WJaWJb_ = 0;
	    double jet_pt_min_ = 25;
	    double jet_eta_cut_ = 2.4;

        //histograms
	    CutFlowTracker tracker_; //tracks how many events pass cuts

	    int mass_bins_, pt_bins_, DR_bins_, eta_bins_, nbins;
//	    double m_bins = {0};
	    double mass_max_, mass_min_;
	    double pt_max_, pt_min_;
	    double DR_max_, DR_min_;
	    double eta_max_, eta_min_;

	    //Gen plots
	    	//DR
	    TH2D *DR_LepBHad_vs_Mtt_hist = 0;
            TH2D *DR_LepBLep_vs_Mtt_hist = 0;
            TH2D *DR_LepWJa_vs_Mtt_hist = 0;
            TH2D *DR_LepWJb_vs_Mtt_hist = 0;
            TH2D *DR_BHadBLep_vs_Mtt_hist = 0;
            TH2D *DR_BHadWJa_vs_Mtt_hist = 0;
            TH2D *DR_BHadWJb_vs_Mtt_hist = 0;
            TH2D *DR_BLepWJa_vs_Mtt_hist = 0;
            TH2D *DR_BLepWJb_vs_Mtt_hist = 0;
            TH2D *DR_WJaWJb_vs_Mtt_hist = 0;
    
            TH1F *DR_LepBHad_hist = 0;
            TH1F *DR_LepBLep_hist = 0;
            TH1F *DR_LepWJa_hist = 0;
            TH1F *DR_LepWJb_hist = 0;
            TH1F *DR_BHadBLep_hist = 0;
            TH1F *DR_BHadWJa_hist = 0;
            TH1F *DR_BHadWJb_hist = 0;
            TH1F *DR_BLepWJa_hist = 0;
            TH1F *DR_BLepWJb_hist = 0;
            TH1F *DR_WJaWJb_hist = 0;

	    TH2D *DRmin_thad_vs_mttbar_hist = 0;
	    TH2D *DRmin_thad_vs_ptthad_hist = 0;
	    TH2D *DRmin_tlep_vs_mttbar_hist = 0;
	    TH2D *DRmin_tlep_vs_pttlep_hist = 0;

	    TH1F *DRmin_thad_hist = 0;
	    TH1F *DRmin_tlep_hist = 0;

	    TH2D *DRmax_thad_vs_mttbar_hist = 0;
	    TH2D *DRmax_thad_vs_ptthad_hist = 0;

	    TH1F *DRmax_thad_hist = 0;
		//Ratio plots
	    TH1F *DRmin_thad_lp4_vs_mttbar_hist = 0; //DRmin thad values < 0.4
	    TH1F *DRmin_thad_lp4_vs_ptthad_hist = 0; //DRmin thad values < 0.4
	    TH1F *DRmin_tlep_lp4_vs_mttbar_hist = 0; //DRmin tlep values < 0.4
	    TH1F *DRmin_tlep_lp4_vs_pttlep_hist = 0; //DRmin tlep values < 0.4

	    TH1F *DRmin_thad_gp4_vs_mttbar_hist = 0; //DRmin thad values > 0.4
	    TH1F *DRmin_thad_gp4_vs_ptthad_hist = 0; //DRmin thad values > 0.4
	    TH1F *DRmin_tlep_gp4_vs_mttbar_hist = 0; //DRmin tlep values > 0.4
	    TH1F *DRmin_tlep_gp4_vs_pttlep_hist = 0; //DRmin tlep values > 0.4

	    TH1F *DRmax_thad_lp4_vs_mttbar_hist = 0; //DRmin thad values < 0.4
	    TH1F *DRmax_thad_lp4_vs_ptthad_hist = 0; //DRmin thad values < 0.4

	    TH1F *DRmax_thad_gp4_vs_mttbar_hist = 0; //DRmin thad values > 0.4
	    TH1F *DRmax_thad_gp4_vs_ptthad_hist = 0; //DRmin thad values > 0.4

	    TH1F *DR_LepBHad_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_LepBLep_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_LepWJa_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_LepWJb_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_BHadBLep_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_BHadWJa_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_BHadWJb_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_BLepWJa_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_BLepWJb_lp4_vs_mttbar_hist = 0;
	    TH1F *DR_WJaWJb_lp4_vs_mttbar_hist = 0;

	    TH1F *DR_LepBHad_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_LepBLep_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_LepWJa_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_LepWJb_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_BHadBLep_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_BHadWJa_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_BHadWJb_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_BLepWJa_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_BLepWJb_gp4_vs_mttbar_hist = 0;
	    TH1F *DR_WJaWJb_gp4_vs_mttbar_hist = 0;
		//Pt
	    TH2D *Pt_Lep_vs_Mtt_hist = 0;
            TH2D *Pt_BLep_vs_Mtt_hist = 0;
            TH2D *Pt_BHad_vs_Mtt_hist = 0;
            TH2D *Pt_WJa_vs_Mtt_hist = 0;
            TH2D *Pt_WJb_vs_Mtt_hist = 0;
    
            TH1F *Pt_Lep_hist = 0;
            TH1F *Pt_BLep_hist = 0;
            TH1F *Pt_BHad_hist = 0;
            TH1F *Pt_WJa_hist = 0;
            TH1F *Pt_WJb_hist = 0;

	    TH1F *Pt_ttbar_hist = 0;
	    TH1F *Pt_thad_hist = 0;
	    TH1F *Pt_tlep_hist = 0;
		//Eta
	    TH2D *Eta_Lep_vs_Mtt_hist = 0;
            TH2D *Eta_BLep_vs_Mtt_hist = 0;
            TH2D *Eta_BHad_vs_Mtt_hist = 0;
            TH2D *Eta_WJa_vs_Mtt_hist = 0;
            TH2D *Eta_WJb_vs_Mtt_hist = 0;

            TH1F *Eta_Lep_hist = 0;
            TH1F *Eta_BLep_hist = 0;
            TH1F *Eta_BHad_hist = 0;
            TH1F *Eta_WJa_hist = 0;
            TH1F *Eta_WJb_hist = 0;

		//system plots
            TH1F *Mass_ttbar_hist = 0;
	    TH1F *Mass_thad_hist = 0;
	    TH1F *Mass_tlep_hist = 0;
	    TH1F *nJets_hist = 0;
	    TH1F *nMatched_objects_3J_hist = 0;//events with 3 jets
	    TH1F *nMatched_objects_4J_hist = 0;//events with 4 jets
	    TH1F *nMatched_objects_5PJ_hist = 0;//events with 5+ jets

		//reco plots
			//DR=0.4
	    TH1F *Matched_perm_BHad_pt_DRP4_hist = 0;//Matched perm pt
	    TH1F *Matched_perm_BLep_pt_DRP4_hist = 0;
	    TH1F *Matched_perm_WJa_pt_DRP4_hist = 0;
	    TH1F *Matched_perm_WJb_pt_DRP4_hist = 0;

	    TH1F *Matched_perm_BHad_eta_DRP4_hist = 0;//Matched perm eta 
	    TH1F *Matched_perm_BLep_eta_DRP4_hist = 0;
	    TH1F *Matched_perm_WJa_eta_DRP4_hist = 0;
	    TH1F *Matched_perm_WJb_eta_DRP4_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_pt_DRP4_hist = 0; //break ttJetsM0 into 700-1000 and > 1000
	    TH1F *Matched_perm_WJb_ttbarM700_pt_DRP4_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_pt_DRP4_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_pt_DRP4_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_frac_p_DRP4_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM700_frac_p_DRP4_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_frac_p_DRP4_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_frac_p_DRP4_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_3J_DRP4_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_3J_DRP4_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_3J_DRP4_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_3J_DRP4_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_3J_DRP4_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_3J_DRP4_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_4J_DRP4_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_4J_DRP4_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_4J_DRP4_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_4J_DRP4_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_4J_DRP4_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_4J_DRP4_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_5PJ_DRP4_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_5PJ_DRP4_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_5PJ_DRP4_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_5PJ_DRP4_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_5PJ_DRP4_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_5PJ_DRP4_hist = 0;

	    TH1F *BTag_lp4_BHad_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp4_BHad_loose_fail_hist = 0;
	    TH1F *BTag_lp4_BHad_medium_pass_hist = 0;
	    TH1F *BTag_lp4_BHad_medium_fail_hist = 0;
	    TH1F *BTag_lp4_BHad_tight_pass_hist = 0;
	    TH1F *BTag_lp4_BHad_tight_fail_hist = 0;

	    TH1F *BTag_lp4_BLep_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp4_BLep_loose_fail_hist = 0;
	    TH1F *BTag_lp4_BLep_medium_pass_hist = 0;
	    TH1F *BTag_lp4_BLep_medium_fail_hist = 0;
	    TH1F *BTag_lp4_BLep_tight_pass_hist = 0;
	    TH1F *BTag_lp4_BLep_tight_fail_hist = 0;

	    TH1F *BTag_gp4_BHad_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp4_BHad_loose_fail_hist = 0;
	    TH1F *BTag_gp4_BHad_medium_pass_hist = 0;
	    TH1F *BTag_gp4_BHad_medium_fail_hist = 0;
	    TH1F *BTag_gp4_BHad_tight_pass_hist = 0;
	    TH1F *BTag_gp4_BHad_tight_fail_hist = 0;

	    TH1F *BTag_gp4_BLep_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp4_BLep_loose_fail_hist = 0;
	    TH1F *BTag_gp4_BLep_medium_pass_hist = 0;
	    TH1F *BTag_gp4_BLep_medium_fail_hist = 0;
	    TH1F *BTag_gp4_BLep_tight_pass_hist = 0;
	    TH1F *BTag_gp4_BLep_tight_fail_hist = 0;

		//Merged jets
	    TH1F *Merged_BHadWJa_perm_and_WJb_mass_DRP4_hist = 0;
	    TH1F *Merged_BHadWJa_perm_mass_DRP4_hist = 0;
	    TH1F *Merged_WJb_mass_DRP4_hist = 0;
	    TH1F *Merged_BHadWJa_perm_DRP4_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJa_massDivpt_DRP4_hist = 0;

	    TH1F *Merged_BHadWJb_perm_and_WJa_mass_DRP4_hist = 0;
	    TH1F *Merged_BHadWJb_perm_mass_DRP4_hist = 0;
	    TH1F *Merged_WJa_mass_DRP4_hist = 0;
	    TH1F *Merged_BHadWJb_perm_DRP4_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJb_massDivpt_DRP4_hist = 0;

	    TH2D *Merged_BHadWJet_mass_vs_BHad_mass_DRP4_hist = 0;

	    TH1F *Merged_BLepWJa_perm_and_WJb_mass_DRP4_hist = 0;
	    TH1F *Merged_BLepWJa_perm_mass_DRP4_hist = 0;
	    TH1F *Merged_BLepWJa_perm_DRP4_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJa_massDivpt_DRP4_hist = 0;

	    TH1F *Merged_BLepWJb_perm_and_WJa_mass_DRP4_hist = 0;
	    TH1F *Merged_BLepWJb_perm_mass_DRP4_hist = 0;
	    TH1F *Merged_BLepWJb_perm_DRP4_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJb_massDivpt_DRP4_hist = 0;

	    TH2D *Merged_BLepWJet_mass_vs_BLep_mass_DRP4_hist = 0;


		//Unmerged jets
	    TH1F *Unmerged_BHad_mass_DRP4_hist = 0;
	    TH1F *Unmerged_WJa_mass_DRP4_hist = 0;
	    TH1F *Unmerged_WJb_mass_DRP4_hist = 0;
	    TH1F *Unmerged_BHad_LepBLep_DRP4_hist = 0;
	    TH1F *Unmerged_WJa_LepBLep_DRP4_hist = 0;
	    TH1F *Unmerged_WJb_LepBLep_DRP4_hist = 0;
	    TH1F *Unmerged_BHad_massDivpt_DRP4_hist = 0;
	    TH1F *Unmerged_WJa_massDivpt_DRP4_hist = 0;
	    TH1F *Unmerged_WJb_massDivpt_DRP4_hist = 0;

	    TH2D *Unmerged_BHadWJet_mass_vs_BHad_mass_DRP4_hist = 0;
	    TH2D *Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP4_hist = 0;

	    TH1F *Unmerged_BLep_mass_DRP4_hist = 0;
	    TH1F *Unmerged_BLep_LepBLep_DRP4_hist = 0;
	    TH1F *Unmerged_BLep_massDivpt_DRP4_hist = 0;

	    TH2D *Unmerged_BLepWJet_mass_vs_BLep_mass_DRP4_hist = 0;
	    TH2D *Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP4_hist = 0;


			//DR=0.5
	    TH1F *Matched_perm_BHad_pt_DRP5_hist = 0;//Matched perm pt
	    TH1F *Matched_perm_BLep_pt_DRP5_hist = 0;
	    TH1F *Matched_perm_WJa_pt_DRP5_hist = 0;
	    TH1F *Matched_perm_WJb_pt_DRP5_hist = 0;

	    TH1F *Matched_perm_BHad_eta_DRP5_hist = 0;//Matched perm eta 
	    TH1F *Matched_perm_BLep_eta_DRP5_hist = 0;
	    TH1F *Matched_perm_WJa_eta_DRP5_hist = 0;
	    TH1F *Matched_perm_WJb_eta_DRP5_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_pt_DRP5_hist = 0; //break ttJetsM0 into 700-1000 and > 1000
	    TH1F *Matched_perm_WJb_ttbarM700_pt_DRP5_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_pt_DRP5_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_pt_DRP5_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_frac_p_DRP5_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM700_frac_p_DRP5_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_frac_p_DRP5_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_frac_p_DRP5_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_3J_DRP5_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_3J_DRP5_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_3J_DRP5_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_3J_DRP5_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_3J_DRP5_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_3J_DRP5_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_4J_DRP5_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_4J_DRP5_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_4J_DRP5_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_4J_DRP5_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_4J_DRP5_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_4J_DRP5_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_5PJ_DRP5_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_5PJ_DRP5_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_5PJ_DRP5_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_5PJ_DRP5_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_5PJ_DRP5_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_5PJ_DRP5_hist = 0;

	    TH1F *BTag_lp5_BHad_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp5_BHad_loose_fail_hist = 0;
	    TH1F *BTag_lp5_BHad_medium_pass_hist = 0;
	    TH1F *BTag_lp5_BHad_medium_fail_hist = 0;
	    TH1F *BTag_lp5_BHad_tight_pass_hist = 0;
	    TH1F *BTag_lp5_BHad_tight_fail_hist = 0;

	    TH1F *BTag_lp5_BLep_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp5_BLep_loose_fail_hist = 0;
	    TH1F *BTag_lp5_BLep_medium_pass_hist = 0;
	    TH1F *BTag_lp5_BLep_medium_fail_hist = 0;
	    TH1F *BTag_lp5_BLep_tight_pass_hist = 0;
	    TH1F *BTag_lp5_BLep_tight_fail_hist = 0;

	    TH1F *BTag_gp5_BHad_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp5_BHad_loose_fail_hist = 0;
	    TH1F *BTag_gp5_BHad_medium_pass_hist = 0;
	    TH1F *BTag_gp5_BHad_medium_fail_hist = 0;
	    TH1F *BTag_gp5_BHad_tight_pass_hist = 0;
	    TH1F *BTag_gp5_BHad_tight_fail_hist = 0;

	    TH1F *BTag_gp5_BLep_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp5_BLep_loose_fail_hist = 0;
	    TH1F *BTag_gp5_BLep_medium_pass_hist = 0;
	    TH1F *BTag_gp5_BLep_medium_fail_hist = 0;
	    TH1F *BTag_gp5_BLep_tight_pass_hist = 0;
	    TH1F *BTag_gp5_BLep_tight_fail_hist = 0;

		//Merged jets
	    TH1F *Merged_BHadWJa_perm_and_WJb_mass_DRP5_hist = 0;
	    TH1F *Merged_BHadWJa_perm_mass_DRP5_hist = 0;
	    TH1F *Merged_WJb_mass_DRP5_hist = 0;
	    TH1F *Merged_BHadWJa_perm_DRP5_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJa_massDivpt_DRP5_hist = 0;

	    TH1F *Merged_BHadWJb_perm_and_WJa_mass_DRP5_hist = 0;
	    TH1F *Merged_BHadWJb_perm_mass_DRP5_hist = 0;
	    TH1F *Merged_WJa_mass_DRP5_hist = 0;
	    TH1F *Merged_BHadWJb_perm_DRP5_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJb_massDivpt_DRP5_hist = 0;

	    TH2D *Merged_BHadWJet_mass_vs_BHad_mass_DRP5_hist = 0;

	    TH1F *Merged_BLepWJa_perm_and_WJb_mass_DRP5_hist = 0;
	    TH1F *Merged_BLepWJa_perm_mass_DRP5_hist = 0;
	    TH1F *Merged_BLepWJa_perm_DRP5_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJa_massDivpt_DRP5_hist = 0;

	    TH1F *Merged_BLepWJb_perm_and_WJa_mass_DRP5_hist = 0;
	    TH1F *Merged_BLepWJb_perm_mass_DRP5_hist = 0;
	    TH1F *Merged_BLepWJb_perm_DRP5_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJb_massDivpt_DRP5_hist = 0;

	    TH2D *Merged_BLepWJet_mass_vs_BLep_mass_DRP5_hist = 0;

		//Unmerged jets
	    TH1F *Unmerged_BHad_mass_DRP5_hist = 0;
	    TH1F *Unmerged_WJa_mass_DRP5_hist = 0;
	    TH1F *Unmerged_WJb_mass_DRP5_hist = 0;
	    TH1F *Unmerged_BHad_LepBLep_DRP5_hist = 0;
	    TH1F *Unmerged_WJa_LepBLep_DRP5_hist = 0;
	    TH1F *Unmerged_WJb_LepBLep_DRP5_hist = 0;
	    TH1F *Unmerged_BHad_massDivpt_DRP5_hist = 0;
	    TH1F *Unmerged_WJa_massDivpt_DRP5_hist = 0;
	    TH1F *Unmerged_WJb_massDivpt_DRP5_hist = 0;

	    TH2D *Unmerged_BHadWJet_mass_vs_BHad_mass_DRP5_hist = 0;
	    TH2D *Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP5_hist = 0;

	    TH1F *Unmerged_BLep_mass_DRP5_hist = 0;
	    TH1F *Unmerged_BLep_LepBLep_DRP5_hist = 0;
	    TH1F *Unmerged_BLep_massDivpt_DRP5_hist = 0;

	    TH2D *Unmerged_BLepWJet_mass_vs_BLep_mass_DRP5_hist = 0;
	    TH2D *Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP5_hist = 0;


			//DR=0.6
	    TH1F *Matched_perm_BHad_pt_DRP6_hist = 0;//Matched perm pt
	    TH1F *Matched_perm_BLep_pt_DRP6_hist = 0;
	    TH1F *Matched_perm_WJa_pt_DRP6_hist = 0;
	    TH1F *Matched_perm_WJb_pt_DRP6_hist = 0;

	    TH1F *Matched_perm_BHad_eta_DRP6_hist = 0;//Matched perm eta 
	    TH1F *Matched_perm_BLep_eta_DRP6_hist = 0;
	    TH1F *Matched_perm_WJa_eta_DRP6_hist = 0;
	    TH1F *Matched_perm_WJb_eta_DRP6_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_pt_DRP6_hist = 0; //break ttJetsM0 into 700-1000 and > 1000
	    TH1F *Matched_perm_WJb_ttbarM700_pt_DRP6_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_pt_DRP6_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_pt_DRP6_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_frac_p_DRP6_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM700_frac_p_DRP6_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_frac_p_DRP6_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_frac_p_DRP6_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_3J_DRP6_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_3J_DRP6_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_3J_DRP6_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_3J_DRP6_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_3J_DRP6_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_3J_DRP6_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_4J_DRP6_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_4J_DRP6_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_4J_DRP6_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_4J_DRP6_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_4J_DRP6_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_4J_DRP6_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_5PJ_DRP6_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_5PJ_DRP6_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_5PJ_DRP6_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_5PJ_DRP6_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_5PJ_DRP6_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_5PJ_DRP6_hist = 0;

	    TH1F *BTag_lp6_BHad_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp6_BHad_loose_fail_hist = 0;
	    TH1F *BTag_lp6_BHad_medium_pass_hist = 0;
	    TH1F *BTag_lp6_BHad_medium_fail_hist = 0;
	    TH1F *BTag_lp6_BHad_tight_pass_hist = 0;
	    TH1F *BTag_lp6_BHad_tight_fail_hist = 0;

	    TH1F *BTag_lp6_BLep_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp6_BLep_loose_fail_hist = 0;
	    TH1F *BTag_lp6_BLep_medium_pass_hist = 0;
	    TH1F *BTag_lp6_BLep_medium_fail_hist = 0;
	    TH1F *BTag_lp6_BLep_tight_pass_hist = 0;
	    TH1F *BTag_lp6_BLep_tight_fail_hist = 0;

	    TH1F *BTag_gp6_BHad_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp6_BHad_loose_fail_hist = 0;
	    TH1F *BTag_gp6_BHad_medium_pass_hist = 0;
	    TH1F *BTag_gp6_BHad_medium_fail_hist = 0;
	    TH1F *BTag_gp6_BHad_tight_pass_hist = 0;
	    TH1F *BTag_gp6_BHad_tight_fail_hist = 0;

	    TH1F *BTag_gp6_BLep_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp6_BLep_loose_fail_hist = 0;
	    TH1F *BTag_gp6_BLep_medium_pass_hist = 0;
	    TH1F *BTag_gp6_BLep_medium_fail_hist = 0;
	    TH1F *BTag_gp6_BLep_tight_pass_hist = 0;
	    TH1F *BTag_gp6_BLep_tight_fail_hist = 0;

		//Merged jets
	    TH1F *Merged_BHadWJa_perm_and_WJb_mass_DRP6_hist = 0;
	    TH1F *Merged_BHadWJa_perm_mass_DRP6_hist = 0;
	    TH1F *Merged_WJb_mass_DRP6_hist = 0;
	    TH1F *Merged_BHadWJa_perm_DRP6_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJa_massDivpt_DRP6_hist = 0;

	    TH1F *Merged_BHadWJb_perm_and_WJa_mass_DRP6_hist = 0;
	    TH1F *Merged_BHadWJb_perm_mass_DRP6_hist = 0;
	    TH1F *Merged_WJa_mass_DRP6_hist = 0;
	    TH1F *Merged_BHadWJb_perm_DRP6_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJb_massDivpt_DRP6_hist = 0;

	    TH2D *Merged_BHadWJet_mass_vs_BHad_mass_DRP6_hist = 0;

	    TH1F *Merged_BLepWJa_perm_and_WJb_mass_DRP6_hist = 0;
	    TH1F *Merged_BLepWJa_perm_mass_DRP6_hist = 0;
	    TH1F *Merged_BLepWJa_perm_DRP6_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJa_massDivpt_DRP6_hist = 0;

	    TH1F *Merged_BLepWJb_perm_and_WJa_mass_DRP6_hist = 0;
	    TH1F *Merged_BLepWJb_perm_mass_DRP6_hist = 0;
	    TH1F *Merged_BLepWJb_perm_DRP6_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJb_massDivpt_DRP6_hist = 0;

	    TH2D *Merged_BLepWJet_mass_vs_BLep_mass_DRP6_hist = 0;

		//Unmerged jets
	    TH1F *Unmerged_BHad_mass_DRP6_hist = 0;
	    TH1F *Unmerged_WJa_mass_DRP6_hist = 0;
	    TH1F *Unmerged_WJb_mass_DRP6_hist = 0;
	    TH1F *Unmerged_BHad_LepBLep_DRP6_hist = 0;
	    TH1F *Unmerged_WJa_LepBLep_DRP6_hist = 0;
	    TH1F *Unmerged_WJb_LepBLep_DRP6_hist = 0;
	    TH1F *Unmerged_BHad_massDivpt_DRP6_hist = 0;
	    TH1F *Unmerged_WJa_massDivpt_DRP6_hist = 0;
	    TH1F *Unmerged_WJb_massDivpt_DRP6_hist = 0;

	    TH2D *Unmerged_BHadWJet_mass_vs_BHad_mass_DRP6_hist = 0;
	    TH2D *Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP6_hist = 0;

	    TH1F *Unmerged_BLep_mass_DRP6_hist = 0;
	    TH1F *Unmerged_BLep_LepBLep_DRP6_hist = 0;
	    TH1F *Unmerged_BLep_massDivpt_DRP6_hist = 0;

	    TH2D *Unmerged_BLepWJet_mass_vs_BLep_mass_DRP6_hist = 0;
	    TH2D *Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP6_hist = 0;


			//DR=0.8
	    TH1F *Matched_perm_BHad_pt_DRP8_hist = 0;//Matched perm pt
	    TH1F *Matched_perm_BLep_pt_DRP8_hist = 0;
	    TH1F *Matched_perm_WJa_pt_DRP8_hist = 0;
	    TH1F *Matched_perm_WJb_pt_DRP8_hist = 0;

	    TH1F *Matched_perm_BHad_eta_DRP8_hist = 0;//Matched perm eta 
	    TH1F *Matched_perm_BLep_eta_DRP8_hist = 0;
	    TH1F *Matched_perm_WJa_eta_DRP8_hist = 0;
	    TH1F *Matched_perm_WJb_eta_DRP8_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_pt_DRP8_hist = 0; //break ttJetsM0 into 700-1000 and > 1000
	    TH1F *Matched_perm_WJb_ttbarM700_pt_DRP8_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_pt_DRP8_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_pt_DRP8_hist = 0;

	    TH1F *Matched_perm_WJa_ttbarM700_frac_p_DRP8_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM700_frac_p_DRP8_hist = 0;
	    TH1F *Matched_perm_WJa_ttbarM1000_frac_p_DRP8_hist = 0;
	    TH1F *Matched_perm_WJb_ttbarM1000_frac_p_DRP8_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_3J_DRP8_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_3J_DRP8_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_3J_DRP8_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_3J_DRP8_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_3J_DRP8_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_3J_DRP8_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_4J_DRP8_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_4J_DRP8_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_4J_DRP8_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_4J_DRP8_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_4J_DRP8_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_4J_DRP8_hist = 0;

	    TH1F *Matched_BHadWJa_ptthad_5PJ_DRP8_hist = 0;
	    TH1F *Matched_BHadWJb_ptthad_5PJ_DRP8_hist = 0;
	    TH1F *Matched_WJaWJb_ptthad_5PJ_DRP8_hist = 0;
	    TH1F *All_Matched_BHadWJa_ptthad_5PJ_DRP8_hist = 0;
	    TH1F *All_Matched_BHadWJb_ptthad_5PJ_DRP8_hist = 0;
	    TH1F *All_Matched_WJaWJb_ptthad_5PJ_DRP8_hist = 0;

	    TH1F *BTag_lp8_BHad_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp8_BHad_loose_fail_hist = 0;
	    TH1F *BTag_lp8_BHad_medium_pass_hist = 0;
	    TH1F *BTag_lp8_BHad_medium_fail_hist = 0;
	    TH1F *BTag_lp8_BHad_tight_pass_hist = 0;
	    TH1F *BTag_lp8_BHad_tight_fail_hist = 0;

	    TH1F *BTag_lp8_BLep_loose_pass_hist = 0;//btagging for merged events
	    TH1F *BTag_lp8_BLep_loose_fail_hist = 0;
	    TH1F *BTag_lp8_BLep_medium_pass_hist = 0;
	    TH1F *BTag_lp8_BLep_medium_fail_hist = 0;
	    TH1F *BTag_lp8_BLep_tight_pass_hist = 0;
	    TH1F *BTag_lp8_BLep_tight_fail_hist = 0;

	    TH1F *BTag_gp8_BHad_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp8_BHad_loose_fail_hist = 0;
	    TH1F *BTag_gp8_BHad_medium_pass_hist = 0;
	    TH1F *BTag_gp8_BHad_medium_fail_hist = 0;
	    TH1F *BTag_gp8_BHad_tight_pass_hist = 0;
	    TH1F *BTag_gp8_BHad_tight_fail_hist = 0;

	    TH1F *BTag_gp8_BLep_loose_pass_hist = 0;//btagging for unmerged events
	    TH1F *BTag_gp8_BLep_loose_fail_hist = 0;
	    TH1F *BTag_gp8_BLep_medium_pass_hist = 0;
	    TH1F *BTag_gp8_BLep_medium_fail_hist = 0;
	    TH1F *BTag_gp8_BLep_tight_pass_hist = 0;
	    TH1F *BTag_gp8_BLep_tight_fail_hist = 0;

		//Merged jets
	    TH1F *Merged_BHadWJa_perm_and_WJb_mass_DRP8_hist = 0;
	    TH1F *Merged_BHadWJa_perm_mass_DRP8_hist = 0;
	    TH1F *Merged_WJb_mass_DRP8_hist = 0;
	    TH1F *Merged_BHadWJa_perm_DRP8_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJa_massDivpt_DRP8_hist = 0;

	    TH1F *Merged_BHadWJb_perm_and_WJa_mass_DRP8_hist = 0;
	    TH1F *Merged_BHadWJb_perm_mass_DRP8_hist = 0;
	    TH1F *Merged_WJa_mass_DRP8_hist = 0;
	    TH1F *Merged_BHadWJb_perm_DRP8_LepBLep_hist = 0;
	    TH1F *Merged_BHadWJb_massDivpt_DRP8_hist = 0;

	    TH2D *Merged_BHadWJet_mass_vs_BHad_mass_DRP8_hist = 0;

	    TH1F *Merged_BLepWJa_perm_and_WJb_mass_DRP8_hist = 0;
	    TH1F *Merged_BLepWJa_perm_mass_DRP8_hist = 0;
	    TH1F *Merged_BLepWJa_perm_DRP8_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJa_massDivpt_DRP8_hist = 0;

	    TH1F *Merged_BLepWJb_perm_and_WJa_mass_DRP8_hist = 0;
	    TH1F *Merged_BLepWJb_perm_mass_DRP8_hist = 0;
	    TH1F *Merged_BLepWJb_perm_DRP8_LepBLep_hist = 0;
	    TH1F *Merged_BLepWJb_massDivpt_DRP8_hist = 0;

	    TH2D *Merged_BLepWJet_mass_vs_BLep_mass_DRP8_hist = 0;

		//Unmerged jets
	    TH1F *Unmerged_BHad_mass_DRP8_hist = 0;
	    TH1F *Unmerged_WJa_mass_DRP8_hist = 0;
	    TH1F *Unmerged_WJb_mass_DRP8_hist = 0;
	    TH1F *Unmerged_BHad_LepBLep_DRP8_hist = 0;
	    TH1F *Unmerged_WJa_LepBLep_DRP8_hist = 0;
	    TH1F *Unmerged_WJb_LepBLep_DRP8_hist = 0;
	    TH1F *Unmerged_BHad_massDivpt_DRP8_hist = 0;
	    TH1F *Unmerged_WJa_massDivpt_DRP8_hist = 0;
	    TH1F *Unmerged_WJb_massDivpt_DRP8_hist = 0;

	    TH2D *Unmerged_BHadWJet_mass_vs_BHad_mass_DRP8_hist = 0;
	    TH2D *Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP8_hist = 0;

	    TH1F *Unmerged_BLep_mass_DRP8_hist = 0;
	    TH1F *Unmerged_BLep_LepBLep_DRP8_hist = 0;
	    TH1F *Unmerged_BLep_massDivpt_DRP8_hist = 0;

	    TH2D *Unmerged_BLepWJet_mass_vs_BLep_mass_DRP8_hist = 0;
	    TH2D *Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP8_hist = 0;

	//switches
	    bool isData_, isTTbar_;

	//selectors and helpers
            TTObjectSelector object_selector_; //selects ttbar objects
            TTGenParticleSelector genp_selector_; //selects generator level objects
	    TTBarSolver solver_; //solves ttbar events
	    TTGenMatcher matcher_; //matches particles on generator level
            TTPermutator permutator_;

	    float evt_weight_;
	    MCWeightProducer mc_weights_;
	    IDJet::BTag cut_tight_b_ = IDJet::BTag::CSVTIGHT;
	    IDJet::BTag cut_medium_b_ = IDJet::BTag::CSVMEDIUM;
	    IDJet::BTag cut_loose_b_ = IDJet::BTag::CSVLOOSE;

	public:
	    jet_effs(const std::string output_filename):
		AnalyzerBase("jet_effs", output_filename),
                tracker_(),
                object_selector_(),
                genp_selector_(TTGenParticleSelector::SelMode::LHE),
		solver_(),
		evt_weight_(1.),
		mc_weights_()
		{
		//set tracker
			tracker_.use_weight(&evt_weight_);
			//find out which sample we're running on

			opts::variables_map &values = URParser::instance().values();
			string output_file = values["output"].as<std::string>();
                        string sample = systematics::get_sample(output_file);
                        bool isSignal = boost::starts_with(sample, "AtoTT") || boost::starts_with(sample, "HtoTT");
                        isTTbar_ = boost::starts_with(sample, "ttJets") || isSignal;

                        isData_  = boost::starts_with(sample, "data");
			Logger::log().debug() << "isData_: " << isData_ << ", isTTbar_: " << isTTbar_ << endl;
			if(isData_) {
				if(sample.find("SingleElectron") != std::string::npos) object_selector_.lepton_type(-1);
				else object_selector_.lepton_type(1);
			}
			if(!isData_) mc_weights_.init(sample);

			//Init Solver
			string filename = "prob_ttJets.permProbComputer.test.root";
                        Logger::log().debug() << "solver file: " << filename << endl;
                        TFile probfile(DataFile(filename).path().c_str());
                        TDirectory *td = (TDirectory*) probfile.Get(systematics::shift_to_name.at(systematics::SysShifts::NOSYS).c_str());
			solver_.Init(td, false, true, true, true, true); //probfile, btag, nusolv,massdis,anghad,anglep

//			cut_tight_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.tightb"));
//			cut_loose_b_ = IDJet::tag(URParser::instance().getCfgPar<string>("best_permutation.looseb"));
                };

        //This method is called once per job at the beginning of the analysis
        //book here your histograms/tree and run every initialization needed
	virtual void begin()
	{
		Logger::log().debug() << "Beginning of begin() " << evt_idx_ << endl;
		outFile_.cd();

		opts::variables_map &values = URParser::instance().values();
		string output_file = values["output"].as<std::string>();
                string sample = systematics::get_sample(output_file);
		Logger::log().debug() << "		" << sample << endl;


	//Histograms
		mass_bins_ = 500;
		pt_bins_ = 200;
		DR_bins_ = 200;
		eta_bins_ = 200;

		DR_min_ = 0.;
		DR_max_ = 8.;
		pt_min_ = 0.;
		pt_max_ = 2000.;
		eta_min_ = -8.;
		eta_max_ = 8.;

		if( sample == "ttJetsM0" ){
			nbins = 22;
//                	m_bins[nbins+1] = {250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
			mass_min_ = 250.;
			mass_max_ = 2000.;
		}
		if( sample == "ttJetsM700" ){
			nbins = 10;
//                	m_bins[nbins+1] = {700, 800, 900, 1000};
			mass_min_ = 700.;
			mass_max_ = 1000.;
		}
		if( sample == "ttJetsM1000" ){
			nbins = 10;
//                	m_bins[nbins+1] = {1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
			mass_min_ = 1000.;
			mass_max_ = 2000.;
		}

		//Gen Plots
			//DR
		DR_LepBHad_vs_Mtt_hist = new TH2D("DR_LepBHad_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_LepBLep_vs_Mtt_hist = new TH2D("DR_LepBLep_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_LepWJa_vs_Mtt_hist = new TH2D("DR_LepWJa_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_LepWJb_vs_Mtt_hist = new TH2D("DR_LepWJb_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_BHadBLep_vs_Mtt_hist = new TH2D("DR_BHadBLep_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_BHadWJa_vs_Mtt_hist = new TH2D("DR_BHadWJa_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_BHadWJb_vs_Mtt_hist = new TH2D("DR_BHadWJb_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_BLepWJa_vs_Mtt_hist = new TH2D("DR_BLepWJa_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_BLepWJb_vs_Mtt_hist = new TH2D("DR_BLepWJb_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
                DR_WJaWJb_vs_Mtt_hist = new TH2D("DR_WJaWJb_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);

                DR_LepBHad_hist = new TH1F("DR_LepBHad_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_LepBLep_hist = new TH1F("DR_LepBLep_hist", "", DR_bins_, DR_min_, DR_max_) ;
                DR_LepWJa_hist = new TH1F("DR_LepWJa_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_LepWJb_hist = new TH1F("DR_LepWJb_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_BHadBLep_hist = new TH1F("DR_BHadBLep_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_BHadWJa_hist = new TH1F("DR_BHadWJa_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_BHadWJb_hist = new TH1F("DR_BHadWJb_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_BLepWJa_hist = new TH1F("DR_BLepWJa_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_BLepWJb_hist = new TH1F("DR_BLepWJb_hist", "", DR_bins_, DR_min_, DR_max_);
                DR_WJaWJb_hist = new TH1F("DR_WJaWJb_hist", "", DR_bins_, DR_min_, DR_max_);

	        DRmin_thad_vs_mttbar_hist = new TH2D("DRmin_thad_vs_mttbar_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
	        DRmin_tlep_vs_mttbar_hist = new TH2D("DRmin_tlep_vs_mttbar_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
	        DRmin_thad_vs_ptthad_hist = new TH2D("DRmin_thad_vs_ptthad_hist", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);
	        DRmin_tlep_vs_pttlep_hist = new TH2D("DRmin_tlep_vs_pttlep_hist", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);

		DRmin_thad_hist = new TH1F("DRmin_thad_hist", "", DR_bins_, DR_min_, DR_max_);
		DRmin_tlep_hist = new TH1F("DRmin_tlep_hist", "", DR_bins_, DR_min_, DR_max_);

		DRmin_thad_lp4_vs_mttbar_hist = new TH1F("DRmin_thad_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DRmin_thad_lp4_vs_ptthad_hist = new TH1F("DRmin_thad_lp4_vs_ptthad_hist", "", nbins, pt_min_, pt_max_);
		DRmin_tlep_lp4_vs_mttbar_hist = new TH1F("DRmin_tlep_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DRmin_tlep_lp4_vs_pttlep_hist = new TH1F("DRmin_tlep_lp4_vs_pttlep_hist", "", nbins, pt_min_, pt_max_);
		DRmin_thad_lp4_vs_mttbar_hist->Sumw2(); 
		DRmin_thad_lp4_vs_ptthad_hist->Sumw2();
		DRmin_tlep_lp4_vs_mttbar_hist->Sumw2();
		DRmin_tlep_lp4_vs_pttlep_hist->Sumw2();

		DRmin_thad_gp4_vs_mttbar_hist = new TH1F("DRmin_thad_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DRmin_thad_gp4_vs_ptthad_hist = new TH1F("DRmin_thad_gp4_vs_ptthad_hist", "", nbins, pt_min_, pt_max_);
		DRmin_tlep_gp4_vs_mttbar_hist = new TH1F("DRmin_tlep_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DRmin_tlep_gp4_vs_pttlep_hist = new TH1F("DRmin_tlep_gp4_vs_pttlep_hist", "", nbins, pt_min_, pt_max_);
		DRmin_thad_gp4_vs_mttbar_hist->Sumw2();
		DRmin_thad_gp4_vs_ptthad_hist->Sumw2();
		DRmin_tlep_gp4_vs_mttbar_hist->Sumw2();
		DRmin_tlep_gp4_vs_pttlep_hist->Sumw2();

	        DRmax_thad_vs_mttbar_hist = new TH2D("DRmax_thad_vs_mttbar_hist", "", mass_bins_, mass_min_, mass_max_, DR_bins_, DR_min_, DR_max_);
	        DRmax_thad_vs_ptthad_hist = new TH2D("DRmax_thad_vs_ptthad_hist", "", pt_bins_, pt_min_, pt_max_, DR_bins_, DR_min_, DR_max_);

		DRmax_thad_hist = new TH1F("DRmax_thad_hist", "", DR_bins_, DR_min_, DR_max_);

		DRmax_thad_lp4_vs_mttbar_hist = new TH1F("DRmax_thad_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DRmax_thad_lp4_vs_ptthad_hist = new TH1F("DRmax_thad_lp4_vs_ptthad_hist", "", nbins, pt_min_, pt_max_);
		DRmax_thad_lp4_vs_mttbar_hist->Sumw2(); 
		DRmax_thad_lp4_vs_ptthad_hist->Sumw2();

		DRmax_thad_gp4_vs_mttbar_hist = new TH1F("DRmax_thad_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DRmax_thad_gp4_vs_ptthad_hist = new TH1F("DRmax_thad_gp4_vs_ptthad_hist", "", nbins, pt_min_, pt_max_);
		DRmax_thad_gp4_vs_mttbar_hist->Sumw2();
		DRmax_thad_gp4_vs_ptthad_hist->Sumw2();

		DR_LepBHad_lp4_vs_mttbar_hist = new TH1F("DR_LepBHad_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepBLep_lp4_vs_mttbar_hist = new TH1F("DR_LepBLep_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepWJa_lp4_vs_mttbar_hist = new TH1F("DR_LepWJa_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepWJb_lp4_vs_mttbar_hist = new TH1F("DR_LepWJb_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BHadBLep_lp4_vs_mttbar_hist = new TH1F("DR_BHadBLep_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BHadWJa_lp4_vs_mttbar_hist = new TH1F("DR_BHadWJa_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BHadWJb_lp4_vs_mttbar_hist = new TH1F("DR_BHadWJb_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BLepWJa_lp4_vs_mttbar_hist = new TH1F("DR_BLepWJa_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BLepWJb_lp4_vs_mttbar_hist = new TH1F("DR_BLepWJb_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_WJaWJb_lp4_vs_mttbar_hist = new TH1F("DR_WJaWJb_lp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepBHad_lp4_vs_mttbar_hist->Sumw2();
		DR_LepBLep_lp4_vs_mttbar_hist->Sumw2(); 
		DR_LepWJa_lp4_vs_mttbar_hist->Sumw2();  
		DR_LepWJb_lp4_vs_mttbar_hist->Sumw2();  
		DR_BHadBLep_lp4_vs_mttbar_hist->Sumw2();
		DR_BHadWJa_lp4_vs_mttbar_hist->Sumw2();
		DR_BHadWJb_lp4_vs_mttbar_hist->Sumw2();
		DR_BLepWJa_lp4_vs_mttbar_hist->Sumw2();
		DR_BLepWJb_lp4_vs_mttbar_hist->Sumw2();
		DR_WJaWJb_lp4_vs_mttbar_hist->Sumw2();  
		
		DR_LepBHad_gp4_vs_mttbar_hist = new TH1F("DR_LepBHad_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepBLep_gp4_vs_mttbar_hist = new TH1F("DR_LepBLep_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepWJa_gp4_vs_mttbar_hist = new TH1F("DR_LepWJa_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepWJb_gp4_vs_mttbar_hist = new TH1F("DR_LepWJb_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BHadBLep_gp4_vs_mttbar_hist = new TH1F("DR_BHadBLep_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BHadWJa_gp4_vs_mttbar_hist = new TH1F("DR_BHadWJa_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BHadWJb_gp4_vs_mttbar_hist = new TH1F("DR_BHadWJb_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BLepWJa_gp4_vs_mttbar_hist = new TH1F("DR_BLepWJa_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_BLepWJb_gp4_vs_mttbar_hist = new TH1F("DR_BLepWJb_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_WJaWJb_gp4_vs_mttbar_hist = new TH1F("DR_WJaWJb_gp4_vs_mttbar_hist", "", nbins, mass_min_, mass_max_);
		DR_LepBHad_gp4_vs_mttbar_hist->Sumw2();
		DR_LepBLep_gp4_vs_mttbar_hist->Sumw2(); 
		DR_LepWJa_gp4_vs_mttbar_hist->Sumw2();  
		DR_LepWJb_gp4_vs_mttbar_hist->Sumw2();  
		DR_BHadBLep_gp4_vs_mttbar_hist->Sumw2();
		DR_BHadWJa_gp4_vs_mttbar_hist->Sumw2();
		DR_BHadWJb_gp4_vs_mttbar_hist->Sumw2();
		DR_BLepWJa_gp4_vs_mttbar_hist->Sumw2();
		DR_BLepWJb_gp4_vs_mttbar_hist->Sumw2();
		DR_WJaWJb_gp4_vs_mttbar_hist->Sumw2();  
			//Pt
		Pt_Lep_vs_Mtt_hist =  new TH2D("Pt_Lep_vs_Mtt_hist", "m_{t#bar t} [GeV]; l p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                Pt_BLep_vs_Mtt_hist = new TH2D("Pt_BLep_vs_Mtt_hist", "m_{t#bar t} [GeV]; b_{l} p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                Pt_BHad_vs_Mtt_hist = new TH2D("Pt_BHad_vs_Mtt_hist", "m_{t#bar t} [GeV]; b_{h} p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                Pt_WJa_vs_Mtt_hist = new TH2D("Pt_WJa_vs_Mtt_hist", "m_{t#bar t} [GeV]; WJa p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);
                Pt_WJb_vs_Mtt_hist = new TH2D("Pt_WJb_vs_Mtt_hist", "m_{t#bar t} [GeV]; WJb p_{t} ", mass_bins_, mass_min_, mass_max_, pt_bins_, pt_min_, pt_max_);

                Pt_Lep_hist = new TH1F("Pt_Lep_hist", "", pt_bins_, pt_min_, pt_max_);
                Pt_BLep_hist = new TH1F("Pt_BLep_hist", "", pt_bins_, pt_min_, pt_max_);
                Pt_BHad_hist = new TH1F("Pt_BHad_hist", "", pt_bins_, pt_min_, pt_max_);
                Pt_WJa_hist = new TH1F("Pt_WJa_hist", "", pt_bins_, pt_min_, pt_max_);
                Pt_WJb_hist = new TH1F("Pt_WJb_hist", "", pt_bins_, pt_min_, pt_max_);

		Pt_ttbar_hist = new TH1F("Pt_ttbar_hist", "", pt_bins_, pt_min_, pt_max_);
		Pt_thad_hist = new TH1F("Pt_thad_hist", "", pt_bins_, pt_min_, pt_max_);
		Pt_tlep_hist = new TH1F("Pt_tlep_hist", "", pt_bins_, pt_min_, pt_max_);
			//Eta
		Eta_Lep_vs_Mtt_hist = new TH2D("Eta_Lep_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                Eta_BLep_vs_Mtt_hist = new TH2D("Eta_BLep_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                Eta_BHad_vs_Mtt_hist = new TH2D("Eta_BHad_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                Eta_WJa_vs_Mtt_hist = new TH2D("Eta_WJa_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);
                Eta_WJb_vs_Mtt_hist = new TH2D("Eta_WJb_vs_Mtt_hist", "", mass_bins_, mass_min_, mass_max_, eta_bins_, eta_min_, eta_max_);

                Eta_Lep_hist = new TH1F("Eta_Lep_hist", "", eta_bins_, eta_min_, eta_max_);
                Eta_BLep_hist = new TH1F("Eta_BLep_hist", "", eta_bins_, eta_min_, eta_max_);
                Eta_BHad_hist = new TH1F("Eta_BHad_hist", "", eta_bins_, eta_min_, eta_max_);
                Eta_WJa_hist = new TH1F("Eta_WJa_hist", "", eta_bins_, eta_min_, eta_max_);
                Eta_WJb_hist = new TH1F("Eta_WJb_hist", "", eta_bins_, eta_min_, eta_max_);

		//System Plots
		Mass_ttbar_hist = new TH1F("Mass_ttbar_hist", "", mass_bins_, mass_min_, mass_max_);
		Mass_thad_hist = new TH1F("Mass_thad_hist", "", mass_bins_, 125., 250.);
		Mass_tlep_hist = new TH1F("Mass_tlep_hist", "", mass_bins_, 125., 250.);
		nJets_hist = new TH1F("nJets_hist", "", 13, 2.5, 15.5);
	    	nMatched_objects_3J_hist = new TH1F("nMatched_objects_3J_hist", "", 4, 0.5, 4.5);
	    	nMatched_objects_4J_hist = new TH1F("nMatched_objects_4J_hist", "", 4, 0.5, 4.5);
	    	nMatched_objects_5PJ_hist = new TH1F("nMatched_objects_5PJ_hist", "", 4, 0.5, 4.5);

	//reco plots
		//DR=0.4
			//Pt
		Matched_perm_BHad_pt_DRP4_hist = new TH1F("Matched_perm_BHad_pt_DRP4", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_BLep_pt_DRP4_hist = new TH1F("Matched_perm_BLep_pt_DRP4", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJa_pt_DRP4_hist = new TH1F("Matched_perm_WJa_pt_DRP4", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJb_pt_DRP4_hist = new TH1F("Matched_perm_WJb_pt_DRP4", "", pt_bins_, pt_min_, pt_max_);

			//Eta
		Matched_perm_BHad_eta_DRP4_hist = new TH1F("Matched_perm_BHad_eta_DRP4", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_BLep_eta_DRP4_hist = new TH1F("Matched_perm_BLep_eta_DRP4", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJa_eta_DRP4_hist = new TH1F("Matched_perm_WJa_eta_DRP4", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJb_eta_DRP4_hist = new TH1F("Matched_perm_WJb_eta_DRP4", "", eta_bins_, eta_min_, eta_max_);

			//break ttJetsM0 into 700-1000 and > 1000
	    	Matched_perm_WJa_ttbarM700_pt_DRP4_hist = new TH1F("Matched_perm_WJa_ttbarM700_pt_DRP4", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM700_pt_DRP4_hist = new TH1F("Matched_perm_WJb_ttbarM700_pt_DRP4", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_perm_WJa_ttbarM1000_pt_DRP4_hist = new TH1F("Matched_perm_WJa_ttbarM1000_pt_DRP4", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM1000_pt_DRP4_hist = new TH1F("Matched_perm_WJb_ttbarM1000_pt_DRP4", "", pt_bins_, pt_min_, pt_max_);

	    	Matched_perm_WJa_ttbarM700_frac_p_DRP4_hist = new TH1F("Matched_perm_WJa_ttbarM700_frac_p_DRP4", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM700_frac_p_DRP4_hist = new TH1F("Matched_perm_WJb_ttbarM700_frac_p_DRP4", "", 10, 0., 5.);
	    	Matched_perm_WJa_ttbarM1000_frac_p_DRP4_hist = new TH1F("Matched_perm_WJa_ttbarM1000_frac_p_DRP4", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM1000_frac_p_DRP4_hist = new TH1F("Matched_perm_WJb_ttbarM1000_frac_p_DRP4", "", 10, 0., 5.);

			//Matched objects
	    	Matched_BHadWJa_ptthad_3J_DRP4_hist = new TH1F("Matched_BHadWJa_ptthad_3J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_3J_DRP4_hist = new TH1F("Matched_BHadWJb_ptthad_3J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_3J_DRP4_hist = new TH1F("Matched_WJaWJb_ptthad_3J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_3J_DRP4_hist = new TH1F("All_Matched_BHadWJa_ptthad_3J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_3J_DRP4_hist = new TH1F("All_Matched_BHadWJb_ptthad_3J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_3J_DRP4_hist = new TH1F("All_Matched_WJaWJb_ptthad_3J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_3J_DRP4_hist->Sumw2();	
		Matched_BHadWJb_ptthad_3J_DRP4_hist->Sumw2();	
		Matched_WJaWJb_ptthad_3J_DRP4_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_3J_DRP4_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_3J_DRP4_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_3J_DRP4_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_4J_DRP4_hist = new TH1F("Matched_BHadWJa_ptthad_4J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_4J_DRP4_hist = new TH1F("Matched_BHadWJb_ptthad_4J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_4J_DRP4_hist = new TH1F("Matched_WJaWJb_ptthad_4J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_4J_DRP4_hist = new TH1F("All_Matched_BHadWJa_ptthad_4J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_4J_DRP4_hist = new TH1F("All_Matched_BHadWJb_ptthad_4J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_4J_DRP4_hist = new TH1F("All_Matched_WJaWJb_ptthad_4J_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_4J_DRP4_hist->Sumw2();	
		Matched_BHadWJb_ptthad_4J_DRP4_hist->Sumw2();	
		Matched_WJaWJb_ptthad_4J_DRP4_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_4J_DRP4_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_4J_DRP4_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_4J_DRP4_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_5PJ_DRP4_hist = new TH1F("Matched_BHadWJa_ptthad_5PJ_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_5PJ_DRP4_hist = new TH1F("Matched_BHadWJb_ptthad_5PJ_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_5PJ_DRP4_hist = new TH1F("Matched_WJaWJb_ptthad_5PJ_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_5PJ_DRP4_hist = new TH1F("All_Matched_BHadWJa_ptthad_5PJ_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_5PJ_DRP4_hist = new TH1F("All_Matched_BHadWJb_ptthad_5PJ_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_5PJ_DRP4_hist = new TH1F("All_Matched_WJaWJb_ptthad_5PJ_DRP4_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_5PJ_DRP4_hist->Sumw2();	
		Matched_BHadWJb_ptthad_5PJ_DRP4_hist->Sumw2();	
		Matched_WJaWJb_ptthad_5PJ_DRP4_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_5PJ_DRP4_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_5PJ_DRP4_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_5PJ_DRP4_hist->Sumw2();	

		BTag_lp4_BHad_loose_pass_hist = new TH1F("BTag_lp4_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BHad_loose_fail_hist = new TH1F("BTag_lp4_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BHad_medium_pass_hist = new TH1F("BTag_lp4_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BHad_medium_fail_hist = new TH1F("BTag_lp4_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BHad_tight_pass_hist = new TH1F("BTag_lp4_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BHad_tight_fail_hist = new TH1F("BTag_lp4_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BHad_loose_pass_hist->Sumw2(); 
        	BTag_lp4_BHad_loose_fail_hist->Sumw2();
        	BTag_lp4_BHad_medium_pass_hist->Sumw2();
        	BTag_lp4_BHad_medium_fail_hist->Sumw2();
        	BTag_lp4_BHad_tight_pass_hist->Sumw2(); 
		BTag_lp4_BHad_tight_fail_hist->Sumw2(); 

		BTag_lp4_BLep_loose_pass_hist = new TH1F("BTag_lp4_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BLep_loose_fail_hist = new TH1F("BTag_lp4_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BLep_medium_pass_hist = new TH1F("BTag_lp4_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BLep_medium_fail_hist = new TH1F("BTag_lp4_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BLep_tight_pass_hist = new TH1F("BTag_lp4_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BLep_tight_fail_hist = new TH1F("BTag_lp4_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp4_BLep_loose_pass_hist->Sumw2(); 
        	BTag_lp4_BLep_loose_fail_hist->Sumw2();
        	BTag_lp4_BLep_medium_pass_hist->Sumw2();
        	BTag_lp4_BLep_medium_fail_hist->Sumw2();
        	BTag_lp4_BLep_tight_pass_hist->Sumw2(); 
		BTag_lp4_BLep_tight_fail_hist->Sumw2(); 

		BTag_gp4_BHad_loose_pass_hist = new TH1F("BTag_gp4_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BHad_loose_fail_hist = new TH1F("BTag_gp4_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BHad_medium_pass_hist = new TH1F("BTag_gp4_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BHad_medium_fail_hist = new TH1F("BTag_gp4_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BHad_tight_pass_hist = new TH1F("BTag_gp4_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BHad_tight_fail_hist = new TH1F("BTag_gp4_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BHad_loose_pass_hist->Sumw2(); 
        	BTag_gp4_BHad_loose_fail_hist->Sumw2();
        	BTag_gp4_BHad_medium_pass_hist->Sumw2();
        	BTag_gp4_BHad_medium_fail_hist->Sumw2();
        	BTag_gp4_BHad_tight_pass_hist->Sumw2(); 
		BTag_gp4_BHad_tight_fail_hist->Sumw2(); 

		BTag_gp4_BLep_loose_pass_hist = new TH1F("BTag_gp4_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BLep_loose_fail_hist = new TH1F("BTag_gp4_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BLep_medium_pass_hist = new TH1F("BTag_gp4_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BLep_medium_fail_hist = new TH1F("BTag_gp4_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BLep_tight_pass_hist = new TH1F("BTag_gp4_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BLep_tight_fail_hist = new TH1F("BTag_gp4_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp4_BLep_loose_pass_hist->Sumw2(); 
        	BTag_gp4_BLep_loose_fail_hist->Sumw2();
        	BTag_gp4_BLep_medium_pass_hist->Sumw2();
        	BTag_gp4_BLep_medium_fail_hist->Sumw2();
        	BTag_gp4_BLep_tight_pass_hist->Sumw2(); 
		BTag_gp4_BLep_tight_fail_hist->Sumw2(); 

			//Merged jets
		Merged_BHadWJa_perm_and_WJb_mass_DRP4_hist = new TH1F("Merged_BHadWJa_perm_and_WJb_mass_DRP4", "", 100, 0., 300.);
		Merged_BHadWJa_perm_mass_DRP4_hist = new TH1F("Merged_BHadWJa_perm_mass_DRP4", "", 100, 0., 300.);
		Merged_WJb_mass_DRP4_hist = new TH1F("Merged_WJb_mass_DRP4", "", 100, 0., 200.);
		Merged_BHadWJa_perm_DRP4_LepBLep_hist = new TH1F("Merged_BHadWJa_perm_DRP4_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJa_massDivpt_DRP4_hist = new TH1F("Merged_BHadWJa_massDivpt_DRP4", "", 50, 0., 0.5);

		Merged_BHadWJb_perm_and_WJa_mass_DRP4_hist = new TH1F("Merged_BHadWJb_perm_and_WJa_mass_DRP4", "", 100, 0., 300.);
		Merged_BHadWJb_perm_mass_DRP4_hist = new TH1F("Merged_BHadWJb_perm_mass_DRP4", "", 100, 0., 300.);
		Merged_WJa_mass_DRP4_hist = new TH1F("Merged_WJa_mass_DRP4", "", 100, 0., 200.);
		Merged_BHadWJb_perm_DRP4_LepBLep_hist = new TH1F("Merged_BHadWJb_perm_DRP4_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJb_massDivpt_DRP4_hist = new TH1F("Merged_BHadWJb_massDivpt_DRP4", "", 50, 0., 0.5);

	        Merged_BHadWJet_mass_vs_BHad_mass_DRP4_hist = new TH2D("Merged_BHadWJet_mass_vs_BHad_mass_DRP4","",100, 0.,300.,100,0., 300.);

		Merged_BLepWJa_perm_and_WJb_mass_DRP4_hist = new TH1F("Merged_BLepWJa_perm_and_WJb_mass_DRP4", "", 100, 0., 300.);
		Merged_BLepWJa_perm_mass_DRP4_hist = new TH1F("Merged_BLepWJa_perm_mass_DRP4", "", 100, 0., 300.);
		Merged_BLepWJa_perm_DRP4_LepBLep_hist = new TH1F("Merged_BLepWJa_perm_DRP4_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJa_massDivpt_DRP4_hist = new TH1F("Merged_BLepWJa_massDivpt_DRP4", "", 50, 0., 0.5);

		Merged_BLepWJb_perm_and_WJa_mass_DRP4_hist = new TH1F("Merged_BLepWJb_perm_and_WJa_mass_DRP4", "", 100, 0., 300.);
		Merged_BLepWJb_perm_mass_DRP4_hist = new TH1F("Merged_BLepWJb_perm_mass_DRP4", "", 100, 0., 300.);
		Merged_BLepWJb_perm_DRP4_LepBLep_hist = new TH1F("Merged_BLepWJb_perm_DRP4_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJb_massDivpt_DRP4_hist = new TH1F("Merged_BLepWJb_massDivpt_DRP4", "", 50, 0., 0.5);

	        Merged_BLepWJet_mass_vs_BLep_mass_DRP4_hist = new TH2D("Merged_BLepWJet_mass_vs_BLep_mass_DRP4","",100, 0.,300.,100,0., 300.);


			//Unmerged jets
		Unmerged_BHad_mass_DRP4_hist = new TH1F("Unmerged_BHad_mass_DRP4", "", 50, 0., 100.);
		Unmerged_WJa_mass_DRP4_hist = new TH1F("Unmerged_WJa_mass_DRP4", "", 50, 0., 100.);
		Unmerged_WJb_mass_DRP4_hist = new TH1F("Unmerged_WJb_mass_DRP4", "", 50, 0., 100.);
		Unmerged_BHad_LepBLep_DRP4_hist = new TH1F("Unmerged_BHad_LepBLep_DRP4", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJa_LepBLep_DRP4_hist = new TH1F("Unmerged_WJa_LepBLep_DRP4", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJb_LepBLep_DRP4_hist = new TH1F("Unmerged_WJb_LepBLep_DRP4", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BHad_massDivpt_DRP4_hist = new TH1F("Unmerged_BHad_massDivpt_DRP4", "", 50, 0., 0.5);
	        Unmerged_WJa_massDivpt_DRP4_hist = new TH1F("Unmerged_WJa_massDivpt_DRP4", "", 50, 0., 0.5);
	        Unmerged_WJb_massDivpt_DRP4_hist = new TH1F("Unmerged_WJb_massDivpt_DRP4", "", 50, 0., 0.5);

	        Unmerged_BHadWJet_mass_vs_BHad_mass_DRP4_hist = new TH2D("Unmerged_BHadWJet_mass_vs_BHad_mass_DRP4","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP4_hist = new TH2D("Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP4", "", 100, 0., 300., 100, 0., 300.);

		Unmerged_BLep_mass_DRP4_hist = new TH1F("Unmerged_BLep_mass_DRP4", "", 50, 0., 100.);
		Unmerged_BLep_LepBLep_DRP4_hist = new TH1F("Unmerged_BLep_LepBLep_DRP4", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BLep_massDivpt_DRP4_hist = new TH1F("Unmerged_BLep_massDivpt_DRP4", "", 50, 0., 0.5);

	        Unmerged_BLepWJet_mass_vs_BLep_mass_DRP4_hist = new TH2D("Unmerged_BLepWJet_mass_vs_BLep_mass_DRP4","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP4_hist = new TH2D("Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP4", "", 100, 0., 300., 100, 0., 300.);


		//DR=0.5
			//Pt
		Matched_perm_BHad_pt_DRP5_hist = new TH1F("Matched_perm_BHad_pt_DRP5", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_BLep_pt_DRP5_hist = new TH1F("Matched_perm_BLep_pt_DRP5", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJa_pt_DRP5_hist = new TH1F("Matched_perm_WJa_pt_DRP5", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJb_pt_DRP5_hist = new TH1F("Matched_perm_WJb_pt_DRP5", "", pt_bins_, pt_min_, pt_max_);

			//Eta
		Matched_perm_BHad_eta_DRP5_hist = new TH1F("Matched_perm_BHad_eta_DRP5", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_BLep_eta_DRP5_hist = new TH1F("Matched_perm_BLep_eta_DRP5", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJa_eta_DRP5_hist = new TH1F("Matched_perm_WJa_eta_DRP5", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJb_eta_DRP5_hist = new TH1F("Matched_perm_WJb_eta_DRP5", "", eta_bins_, eta_min_, eta_max_);

			//break ttJetsM0 into 700-1000 and > 1000
	    	Matched_perm_WJa_ttbarM700_pt_DRP5_hist = new TH1F("Matched_perm_WJa_ttbarM700_pt_DRP5", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM700_pt_DRP5_hist = new TH1F("Matched_perm_WJb_ttbarM700_pt_DRP5", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_perm_WJa_ttbarM1000_pt_DRP5_hist = new TH1F("Matched_perm_WJa_ttbarM1000_pt_DRP5", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM1000_pt_DRP5_hist = new TH1F("Matched_perm_WJb_ttbarM1000_pt_DRP5", "", pt_bins_, pt_min_, pt_max_);

	    	Matched_perm_WJa_ttbarM700_frac_p_DRP5_hist = new TH1F("Matched_perm_WJa_ttbarM700_frac_p_DRP5", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM700_frac_p_DRP5_hist = new TH1F("Matched_perm_WJb_ttbarM700_frac_p_DRP5", "", 10, 0., 5.);
	    	Matched_perm_WJa_ttbarM1000_frac_p_DRP5_hist = new TH1F("Matched_perm_WJa_ttbarM1000_frac_p_DRP5", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM1000_frac_p_DRP5_hist = new TH1F("Matched_perm_WJb_ttbarM1000_frac_p_DRP5", "", 10, 0., 5.);

			//Matched objects
	    	Matched_BHadWJa_ptthad_3J_DRP5_hist = new TH1F("Matched_BHadWJa_ptthad_3J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_3J_DRP5_hist = new TH1F("Matched_BHadWJb_ptthad_3J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_3J_DRP5_hist = new TH1F("Matched_WJaWJb_ptthad_3J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_3J_DRP5_hist = new TH1F("All_Matched_BHadWJa_ptthad_3J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_3J_DRP5_hist = new TH1F("All_Matched_BHadWJb_ptthad_3J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_3J_DRP5_hist = new TH1F("All_Matched_WJaWJb_ptthad_3J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_3J_DRP5_hist->Sumw2();	
		Matched_BHadWJb_ptthad_3J_DRP5_hist->Sumw2();	
		Matched_WJaWJb_ptthad_3J_DRP5_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_3J_DRP5_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_3J_DRP5_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_3J_DRP5_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_4J_DRP5_hist = new TH1F("Matched_BHadWJa_ptthad_4J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_4J_DRP5_hist = new TH1F("Matched_BHadWJb_ptthad_4J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_4J_DRP5_hist = new TH1F("Matched_WJaWJb_ptthad_4J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_4J_DRP5_hist = new TH1F("All_Matched_BHadWJa_ptthad_4J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_4J_DRP5_hist = new TH1F("All_Matched_BHadWJb_ptthad_4J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_4J_DRP5_hist = new TH1F("All_Matched_WJaWJb_ptthad_4J_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_4J_DRP5_hist->Sumw2();	
		Matched_BHadWJb_ptthad_4J_DRP5_hist->Sumw2();	
		Matched_WJaWJb_ptthad_4J_DRP5_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_4J_DRP5_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_4J_DRP5_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_4J_DRP5_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_5PJ_DRP5_hist = new TH1F("Matched_BHadWJa_ptthad_5PJ_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_5PJ_DRP5_hist = new TH1F("Matched_BHadWJb_ptthad_5PJ_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_5PJ_DRP5_hist = new TH1F("Matched_WJaWJb_ptthad_5PJ_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_5PJ_DRP5_hist = new TH1F("All_Matched_BHadWJa_ptthad_5PJ_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_5PJ_DRP5_hist = new TH1F("All_Matched_BHadWJb_ptthad_5PJ_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_5PJ_DRP5_hist = new TH1F("All_Matched_WJaWJb_ptthad_5PJ_DRP5_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_5PJ_DRP5_hist->Sumw2();	
		Matched_BHadWJb_ptthad_5PJ_DRP5_hist->Sumw2();	
		Matched_WJaWJb_ptthad_5PJ_DRP5_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_5PJ_DRP5_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_5PJ_DRP5_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_5PJ_DRP5_hist->Sumw2();	

		BTag_lp5_BHad_loose_pass_hist = new TH1F("BTag_lp5_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BHad_loose_fail_hist = new TH1F("BTag_lp5_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BHad_medium_pass_hist = new TH1F("BTag_lp5_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BHad_medium_fail_hist = new TH1F("BTag_lp5_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BHad_tight_pass_hist = new TH1F("BTag_lp5_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BHad_tight_fail_hist = new TH1F("BTag_lp5_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BHad_loose_pass_hist->Sumw2(); 
        	BTag_lp5_BHad_loose_fail_hist->Sumw2();
        	BTag_lp5_BHad_medium_pass_hist->Sumw2();
        	BTag_lp5_BHad_medium_fail_hist->Sumw2();
        	BTag_lp5_BHad_tight_pass_hist->Sumw2(); 
		BTag_lp5_BHad_tight_fail_hist->Sumw2(); 

		BTag_lp5_BLep_loose_pass_hist = new TH1F("BTag_lp5_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BLep_loose_fail_hist = new TH1F("BTag_lp5_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BLep_medium_pass_hist = new TH1F("BTag_lp5_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BLep_medium_fail_hist = new TH1F("BTag_lp5_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BLep_tight_pass_hist = new TH1F("BTag_lp5_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BLep_tight_fail_hist = new TH1F("BTag_lp5_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp5_BLep_loose_pass_hist->Sumw2(); 
        	BTag_lp5_BLep_loose_fail_hist->Sumw2();
        	BTag_lp5_BLep_medium_pass_hist->Sumw2();
        	BTag_lp5_BLep_medium_fail_hist->Sumw2();
        	BTag_lp5_BLep_tight_pass_hist->Sumw2(); 
		BTag_lp5_BLep_tight_fail_hist->Sumw2(); 

		BTag_gp5_BHad_loose_pass_hist = new TH1F("BTag_gp5_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BHad_loose_fail_hist = new TH1F("BTag_gp5_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BHad_medium_pass_hist = new TH1F("BTag_gp5_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BHad_medium_fail_hist = new TH1F("BTag_gp5_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BHad_tight_pass_hist = new TH1F("BTag_gp5_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BHad_tight_fail_hist = new TH1F("BTag_gp5_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BHad_loose_pass_hist->Sumw2(); 
        	BTag_gp5_BHad_loose_fail_hist->Sumw2();
        	BTag_gp5_BHad_medium_pass_hist->Sumw2();
        	BTag_gp5_BHad_medium_fail_hist->Sumw2();
        	BTag_gp5_BHad_tight_pass_hist->Sumw2(); 
		BTag_gp5_BHad_tight_fail_hist->Sumw2(); 

		BTag_gp5_BLep_loose_pass_hist = new TH1F("BTag_gp5_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BLep_loose_fail_hist = new TH1F("BTag_gp5_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BLep_medium_pass_hist = new TH1F("BTag_gp5_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BLep_medium_fail_hist = new TH1F("BTag_gp5_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BLep_tight_pass_hist = new TH1F("BTag_gp5_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BLep_tight_fail_hist = new TH1F("BTag_gp5_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp5_BLep_loose_pass_hist->Sumw2(); 
        	BTag_gp5_BLep_loose_fail_hist->Sumw2();
        	BTag_gp5_BLep_medium_pass_hist->Sumw2();
        	BTag_gp5_BLep_medium_fail_hist->Sumw2();
        	BTag_gp5_BLep_tight_pass_hist->Sumw2(); 
		BTag_gp5_BLep_tight_fail_hist->Sumw2(); 

			//Merged jets
		Merged_BHadWJa_perm_and_WJb_mass_DRP5_hist = new TH1F("Merged_BHadWJa_perm_and_WJb_mass_DRP5", "", 100, 0., 500.);
		Merged_BHadWJa_perm_mass_DRP5_hist = new TH1F("Merged_BHadWJa_perm_mass_DRP5", "", 100, 0., 500.);
		Merged_WJb_mass_DRP5_hist = new TH1F("Merged_WJb_mass_DRP5", "", 100, 0., 300.);
		Merged_BHadWJa_perm_DRP5_LepBLep_hist = new TH1F("Merged_BHadWJa_perm_DRP5_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJa_massDivpt_DRP5_hist = new TH1F("Merged_BHadWJa_massDivpt_DRP5", "", 50, 0., 0.5);

		Merged_BHadWJb_perm_and_WJa_mass_DRP5_hist = new TH1F("Merged_BHadWJb_perm_and_WJa_mass_DRP5", "", 100, 0., 500.);
		Merged_BHadWJb_perm_mass_DRP5_hist = new TH1F("Merged_BHadWJb_perm_mass_DRP5", "", 100, 0., 500.);
		Merged_WJa_mass_DRP5_hist = new TH1F("Merged_WJa_mass_DRP5", "", 100, 0., 300.);
		Merged_BHadWJb_perm_DRP5_LepBLep_hist = new TH1F("Merged_BHadWJb_perm_DRP5_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJb_massDivpt_DRP5_hist = new TH1F("Merged_BHadWJb_massDivpt_DRP5", "", 50, 0., 0.5);

	        Merged_BHadWJet_mass_vs_BHad_mass_DRP5_hist = new TH2D("Merged_BHadWJet_mass_vs_BHad_mass_DRP5","",100, 0.,300.,100,0., 300.);

		Merged_BLepWJa_perm_and_WJb_mass_DRP5_hist = new TH1F("Merged_BLepWJa_perm_and_WJb_mass_DRP5", "", 100, 0., 300.);
		Merged_BLepWJa_perm_mass_DRP5_hist = new TH1F("Merged_BLepWJa_perm_mass_DRP5", "", 100, 0., 300.);
		Merged_BLepWJa_perm_DRP5_LepBLep_hist = new TH1F("Merged_BLepWJa_perm_DRP5_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJa_massDivpt_DRP5_hist = new TH1F("Merged_BLepWJa_massDivpt_DRP5", "", 50, 0., 0.5);

		Merged_BLepWJb_perm_and_WJa_mass_DRP5_hist = new TH1F("Merged_BLepWJb_perm_and_WJa_mass_DRP5", "", 100, 0., 300.);
		Merged_BLepWJb_perm_mass_DRP5_hist = new TH1F("Merged_BLepWJb_perm_mass_DRP5", "", 100, 0., 300.);
		Merged_BLepWJb_perm_DRP5_LepBLep_hist = new TH1F("Merged_BLepWJb_perm_DRP5_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJb_massDivpt_DRP5_hist = new TH1F("Merged_BLepWJb_massDivpt_DRP5", "", 50, 0., 0.5);

	        Merged_BLepWJet_mass_vs_BLep_mass_DRP5_hist = new TH2D("Merged_BLepWJet_mass_vs_BLep_mass_DRP5","",100, 0.,300.,100,0., 300.);

			//Unmerged jets
		Unmerged_BHad_mass_DRP5_hist = new TH1F("Unmerged_BHad_mass_DRP5", "", 50, 0., 100.);
		Unmerged_WJa_mass_DRP5_hist = new TH1F("Unmerged_WJa_mass_DRP5", "", 50, 0., 100.);
		Unmerged_WJb_mass_DRP5_hist = new TH1F("Unmerged_WJb_mass_DRP5", "", 50, 0., 100.);
		Unmerged_BHad_LepBLep_DRP5_hist = new TH1F("Unmerged_BHad_LepBLep_DRP5", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJa_LepBLep_DRP5_hist = new TH1F("Unmerged_WJa_LepBLep_DRP5", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJb_LepBLep_DRP5_hist = new TH1F("Unmerged_WJb_LepBLep_DRP5", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BHad_massDivpt_DRP5_hist = new TH1F("Unmerged_BHad_massDivpt_DRP5", "", 50, 0., 0.5);
	        Unmerged_WJa_massDivpt_DRP5_hist = new TH1F("Unmerged_WJa_massDivpt_DRP5", "", 50, 0., 0.5);
	        Unmerged_WJb_massDivpt_DRP5_hist = new TH1F("Unmerged_WJb_massDivpt_DRP5", "", 50, 0., 0.5);

	        Unmerged_BHadWJet_mass_vs_BHad_mass_DRP5_hist = new TH2D("Unmerged_BHadWJet_mass_vs_BHad_mass_DRP5","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP5_hist = new TH2D("Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP5", "", 100, 0., 300., 100, 0., 300.);

		Unmerged_BLep_mass_DRP5_hist = new TH1F("Unmerged_BLep_mass_DRP5", "", 50, 0., 100.);
		Unmerged_BLep_LepBLep_DRP5_hist = new TH1F("Unmerged_BLep_LepBLep_DRP5", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BLep_massDivpt_DRP5_hist = new TH1F("Unmerged_BLep_massDivpt_DRP5", "", 50, 0., 0.5);

	        Unmerged_BLepWJet_mass_vs_BLep_mass_DRP5_hist = new TH2D("Unmerged_BLepWJet_mass_vs_BLep_mass_DRP5","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP5_hist = new TH2D("Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP5", "", 100, 0., 300., 100, 0., 300.);

		//DR=0.6
			//Pt
		Matched_perm_BHad_pt_DRP6_hist = new TH1F("Matched_perm_BHad_pt_DRP6", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_BLep_pt_DRP6_hist = new TH1F("Matched_perm_BLep_pt_DRP6", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJa_pt_DRP6_hist = new TH1F("Matched_perm_WJa_pt_DRP6", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJb_pt_DRP6_hist = new TH1F("Matched_perm_WJb_pt_DRP6", "", pt_bins_, pt_min_, pt_max_);

			//Eta
		Matched_perm_BHad_eta_DRP6_hist = new TH1F("Matched_perm_BHad_eta_DRP6", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_BLep_eta_DRP6_hist = new TH1F("Matched_perm_BLep_eta_DRP6", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJa_eta_DRP6_hist = new TH1F("Matched_perm_WJa_eta_DRP6", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJb_eta_DRP6_hist = new TH1F("Matched_perm_WJb_eta_DRP6", "", eta_bins_, eta_min_, eta_max_);

			//break ttJetsM0 into 700-1000 and > 1000
	    	Matched_perm_WJa_ttbarM700_pt_DRP6_hist = new TH1F("Matched_perm_WJa_ttbarM700_pt_DRP6", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM700_pt_DRP6_hist = new TH1F("Matched_perm_WJb_ttbarM700_pt_DRP6", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_perm_WJa_ttbarM1000_pt_DRP6_hist = new TH1F("Matched_perm_WJa_ttbarM1000_pt_DRP6", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM1000_pt_DRP6_hist = new TH1F("Matched_perm_WJb_ttbarM1000_pt_DRP6", "", pt_bins_, pt_min_, pt_max_);

	    	Matched_perm_WJa_ttbarM700_frac_p_DRP6_hist = new TH1F("Matched_perm_WJa_ttbarM700_frac_p_DRP6", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM700_frac_p_DRP6_hist = new TH1F("Matched_perm_WJb_ttbarM700_frac_p_DRP6", "", 10, 0., 5.);
	    	Matched_perm_WJa_ttbarM1000_frac_p_DRP6_hist = new TH1F("Matched_perm_WJa_ttbarM1000_frac_p_DRP6", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM1000_frac_p_DRP6_hist = new TH1F("Matched_perm_WJb_ttbarM1000_frac_p_DRP6", "", 10, 0., 5.);

			//Matched objects
	    	Matched_BHadWJa_ptthad_3J_DRP6_hist = new TH1F("Matched_BHadWJa_ptthad_3J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_3J_DRP6_hist = new TH1F("Matched_BHadWJb_ptthad_3J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_3J_DRP6_hist = new TH1F("Matched_WJaWJb_ptthad_3J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_3J_DRP6_hist = new TH1F("All_Matched_BHadWJa_ptthad_3J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_3J_DRP6_hist = new TH1F("All_Matched_BHadWJb_ptthad_3J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_3J_DRP6_hist = new TH1F("All_Matched_WJaWJb_ptthad_3J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_3J_DRP6_hist->Sumw2();	
		Matched_BHadWJb_ptthad_3J_DRP6_hist->Sumw2();	
		Matched_WJaWJb_ptthad_3J_DRP6_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_3J_DRP6_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_3J_DRP6_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_3J_DRP6_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_4J_DRP6_hist = new TH1F("Matched_BHadWJa_ptthad_4J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_4J_DRP6_hist = new TH1F("Matched_BHadWJb_ptthad_4J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_4J_DRP6_hist = new TH1F("Matched_WJaWJb_ptthad_4J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_4J_DRP6_hist = new TH1F("All_Matched_BHadWJa_ptthad_4J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_4J_DRP6_hist = new TH1F("All_Matched_BHadWJb_ptthad_4J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_4J_DRP6_hist = new TH1F("All_Matched_WJaWJb_ptthad_4J_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_4J_DRP6_hist->Sumw2();	
		Matched_BHadWJb_ptthad_4J_DRP6_hist->Sumw2();	
		Matched_WJaWJb_ptthad_4J_DRP6_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_4J_DRP6_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_4J_DRP6_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_4J_DRP6_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_5PJ_DRP6_hist = new TH1F("Matched_BHadWJa_ptthad_5PJ_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_5PJ_DRP6_hist = new TH1F("Matched_BHadWJb_ptthad_5PJ_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_5PJ_DRP6_hist = new TH1F("Matched_WJaWJb_ptthad_5PJ_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_5PJ_DRP6_hist = new TH1F("All_Matched_BHadWJa_ptthad_5PJ_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_5PJ_DRP6_hist = new TH1F("All_Matched_BHadWJb_ptthad_5PJ_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_5PJ_DRP6_hist = new TH1F("All_Matched_WJaWJb_ptthad_5PJ_DRP6_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_5PJ_DRP6_hist->Sumw2();	
		Matched_BHadWJb_ptthad_5PJ_DRP6_hist->Sumw2();	
		Matched_WJaWJb_ptthad_5PJ_DRP6_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_5PJ_DRP6_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_5PJ_DRP6_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_5PJ_DRP6_hist->Sumw2();	

		BTag_lp6_BHad_loose_pass_hist = new TH1F("BTag_lp6_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BHad_loose_fail_hist = new TH1F("BTag_lp6_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BHad_medium_pass_hist = new TH1F("BTag_lp6_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BHad_medium_fail_hist = new TH1F("BTag_lp6_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BHad_tight_pass_hist = new TH1F("BTag_lp6_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BHad_tight_fail_hist = new TH1F("BTag_lp6_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BHad_loose_pass_hist->Sumw2(); 
        	BTag_lp6_BHad_loose_fail_hist->Sumw2();
        	BTag_lp6_BHad_medium_pass_hist->Sumw2();
        	BTag_lp6_BHad_medium_fail_hist->Sumw2();
        	BTag_lp6_BHad_tight_pass_hist->Sumw2(); 
		BTag_lp6_BHad_tight_fail_hist->Sumw2(); 

		BTag_lp6_BLep_loose_pass_hist = new TH1F("BTag_lp6_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BLep_loose_fail_hist = new TH1F("BTag_lp6_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BLep_medium_pass_hist = new TH1F("BTag_lp6_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BLep_medium_fail_hist = new TH1F("BTag_lp6_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BLep_tight_pass_hist = new TH1F("BTag_lp6_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BLep_tight_fail_hist = new TH1F("BTag_lp6_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp6_BLep_loose_pass_hist->Sumw2(); 
        	BTag_lp6_BLep_loose_fail_hist->Sumw2();
        	BTag_lp6_BLep_medium_pass_hist->Sumw2();
        	BTag_lp6_BLep_medium_fail_hist->Sumw2();
        	BTag_lp6_BLep_tight_pass_hist->Sumw2(); 
		BTag_lp6_BLep_tight_fail_hist->Sumw2(); 

		BTag_gp6_BHad_loose_pass_hist = new TH1F("BTag_gp6_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BHad_loose_fail_hist = new TH1F("BTag_gp6_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BHad_medium_pass_hist = new TH1F("BTag_gp6_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BHad_medium_fail_hist = new TH1F("BTag_gp6_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BHad_tight_pass_hist = new TH1F("BTag_gp6_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BHad_tight_fail_hist = new TH1F("BTag_gp6_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BHad_loose_pass_hist->Sumw2(); 
        	BTag_gp6_BHad_loose_fail_hist->Sumw2();
        	BTag_gp6_BHad_medium_pass_hist->Sumw2();
        	BTag_gp6_BHad_medium_fail_hist->Sumw2();
        	BTag_gp6_BHad_tight_pass_hist->Sumw2(); 
		BTag_gp6_BHad_tight_fail_hist->Sumw2(); 

		BTag_gp6_BLep_loose_pass_hist = new TH1F("BTag_gp6_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BLep_loose_fail_hist = new TH1F("BTag_gp6_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BLep_medium_pass_hist = new TH1F("BTag_gp6_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BLep_medium_fail_hist = new TH1F("BTag_gp6_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BLep_tight_pass_hist = new TH1F("BTag_gp6_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BLep_tight_fail_hist = new TH1F("BTag_gp6_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp6_BLep_loose_pass_hist->Sumw2(); 
        	BTag_gp6_BLep_loose_fail_hist->Sumw2();
        	BTag_gp6_BLep_medium_pass_hist->Sumw2();
        	BTag_gp6_BLep_medium_fail_hist->Sumw2();
        	BTag_gp6_BLep_tight_pass_hist->Sumw2(); 
		BTag_gp6_BLep_tight_fail_hist->Sumw2(); 

			//Merged jets
		Merged_BHadWJa_perm_and_WJb_mass_DRP6_hist = new TH1F("Merged_BHadWJa_perm_and_WJb_mass_DRP6", "", 100, 0., 500.);
		Merged_BHadWJa_perm_mass_DRP6_hist = new TH1F("Merged_BHadWJa_perm_mass_DRP6", "", 100, 0., 500.);
		Merged_WJb_mass_DRP6_hist = new TH1F("Merged_WJb_mass_DRP6", "", 100, 0., 300.);
		Merged_BHadWJa_perm_DRP6_LepBLep_hist = new TH1F("Merged_BHadWJa_perm_DRP6_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJa_massDivpt_DRP6_hist = new TH1F("Merged_BHadWJa_massDivpt_DRP6", "", 50, 0., 0.5);

		Merged_BHadWJb_perm_and_WJa_mass_DRP6_hist = new TH1F("Merged_BHadWJb_perm_and_WJa_mass_DRP6", "", 100, 0., 500.);
		Merged_BHadWJb_perm_mass_DRP6_hist = new TH1F("Merged_BHadWJb_perm_mass_DRP6", "", 100, 0., 500.);
		Merged_WJa_mass_DRP6_hist = new TH1F("Merged_WJa_mass_DRP6", "", 100, 0., 300.);
		Merged_BHadWJb_perm_DRP6_LepBLep_hist = new TH1F("Merged_BHadWJb_perm_DRP6_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJb_massDivpt_DRP6_hist = new TH1F("Merged_BHadWJb_massDivpt_DRP6", "", 50, 0., 0.5);

	        Merged_BHadWJet_mass_vs_BHad_mass_DRP6_hist = new TH2D("Merged_BHadWJet_mass_vs_BHad_mass_DRP6","",100, 0.,300.,100,0., 300.);

		Merged_BLepWJa_perm_and_WJb_mass_DRP6_hist = new TH1F("Merged_BLepWJa_perm_and_WJb_mass_DRP6", "", 100, 0., 300.);
		Merged_BLepWJa_perm_mass_DRP6_hist = new TH1F("Merged_BLepWJa_perm_mass_DRP6", "", 100, 0., 300.);
		Merged_BLepWJa_perm_DRP6_LepBLep_hist = new TH1F("Merged_BLepWJa_perm_DRP6_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJa_massDivpt_DRP6_hist = new TH1F("Merged_BLepWJa_massDivpt_DRP6", "", 50, 0., 0.5);

		Merged_BLepWJb_perm_and_WJa_mass_DRP6_hist = new TH1F("Merged_BLepWJb_perm_and_WJa_mass_DRP6", "", 100, 0., 300.);
		Merged_BLepWJb_perm_mass_DRP6_hist = new TH1F("Merged_BLepWJb_perm_mass_DRP6", "", 100, 0., 300.);
		Merged_BLepWJb_perm_DRP6_LepBLep_hist = new TH1F("Merged_BLepWJb_perm_DRP6_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJb_massDivpt_DRP6_hist = new TH1F("Merged_BLepWJb_massDivpt_DRP6", "", 50, 0., 0.5);

	        Merged_BLepWJet_mass_vs_BLep_mass_DRP6_hist = new TH2D("Merged_BLepWJet_mass_vs_BLep_mass_DRP6","",100, 0.,300.,100,0., 300.);

			//Unmerged jets
		Unmerged_BHad_mass_DRP6_hist = new TH1F("Unmerged_BHad_mass_DRP6", "", 50, 0., 100.);
		Unmerged_WJa_mass_DRP6_hist = new TH1F("Unmerged_WJa_mass_DRP6", "", 50, 0., 100.);
		Unmerged_WJb_mass_DRP6_hist = new TH1F("Unmerged_WJb_mass_DRP6", "", 50, 0., 100.);
		Unmerged_BHad_LepBLep_DRP6_hist = new TH1F("Unmerged_BHad_LepBLep_DRP6", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJa_LepBLep_DRP6_hist = new TH1F("Unmerged_WJa_LepBLep_DRP6", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJb_LepBLep_DRP6_hist = new TH1F("Unmerged_WJb_LepBLep_DRP6", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BHad_massDivpt_DRP6_hist = new TH1F("Unmerged_BHad_massDivpt_DRP6", "", 50, 0., 0.5);
	        Unmerged_WJa_massDivpt_DRP6_hist = new TH1F("Unmerged_WJa_massDivpt_DRP6", "", 50, 0., 0.5);
	        Unmerged_WJb_massDivpt_DRP6_hist = new TH1F("Unmerged_WJb_massDivpt_DRP6", "", 50, 0., 0.5);

	        Unmerged_BHadWJet_mass_vs_BHad_mass_DRP6_hist = new TH2D("Unmerged_BHadWJet_mass_vs_BHad_mass_DRP6","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP6_hist = new TH2D("Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP6", "", 100, 0., 300., 100, 0., 300.);

		Unmerged_BLep_mass_DRP6_hist = new TH1F("Unmerged_BLep_mass_DRP6", "", 50, 0., 100.);
		Unmerged_BLep_LepBLep_DRP6_hist = new TH1F("Unmerged_BLep_LepBLep_DRP6", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BLep_massDivpt_DRP6_hist = new TH1F("Unmerged_BLep_massDivpt_DRP6", "", 50, 0., 0.5);

	        Unmerged_BLepWJet_mass_vs_BLep_mass_DRP6_hist = new TH2D("Unmerged_BLepWJet_mass_vs_BLep_mass_DRP6","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP6_hist = new TH2D("Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP6", "", 100, 0., 300., 100, 0., 300.);

		//DR=0.8
			//Pt
		Matched_perm_BHad_pt_DRP8_hist = new TH1F("Matched_perm_BHad_pt_DRP8", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_BLep_pt_DRP8_hist = new TH1F("Matched_perm_BLep_pt_DRP8", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJa_pt_DRP8_hist = new TH1F("Matched_perm_WJa_pt_DRP8", "", pt_bins_, pt_min_, pt_max_);
		Matched_perm_WJb_pt_DRP8_hist = new TH1F("Matched_perm_WJb_pt_DRP8", "", pt_bins_, pt_min_, pt_max_);

			//Eta
		Matched_perm_BHad_eta_DRP8_hist = new TH1F("Matched_perm_BHad_eta_DRP8", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_BLep_eta_DRP8_hist = new TH1F("Matched_perm_BLep_eta_DRP8", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJa_eta_DRP8_hist = new TH1F("Matched_perm_WJa_eta_DRP8", "", eta_bins_, eta_min_, eta_max_);
		Matched_perm_WJb_eta_DRP8_hist = new TH1F("Matched_perm_WJb_eta_DRP8", "", eta_bins_, eta_min_, eta_max_);

			//break ttJetsM0 into 700-1000 and > 1000
	    	Matched_perm_WJa_ttbarM700_pt_DRP8_hist = new TH1F("Matched_perm_WJa_ttbarM700_pt_DRP8", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM700_pt_DRP8_hist = new TH1F("Matched_perm_WJb_ttbarM700_pt_DRP8", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_perm_WJa_ttbarM1000_pt_DRP8_hist = new TH1F("Matched_perm_WJa_ttbarM1000_pt_DRP8", "", pt_bins_, pt_min_, pt_max_); 
	    	Matched_perm_WJb_ttbarM1000_pt_DRP8_hist = new TH1F("Matched_perm_WJb_ttbarM1000_pt_DRP8", "", pt_bins_, pt_min_, pt_max_);

	    	Matched_perm_WJa_ttbarM700_frac_p_DRP8_hist = new TH1F("Matched_perm_WJa_ttbarM700_frac_p_DRP8", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM700_frac_p_DRP8_hist = new TH1F("Matched_perm_WJb_ttbarM700_frac_p_DRP8", "", 10, 0., 5.);
	    	Matched_perm_WJa_ttbarM1000_frac_p_DRP8_hist = new TH1F("Matched_perm_WJa_ttbarM1000_frac_p_DRP8", "", 10, 0., 5.);
	    	Matched_perm_WJb_ttbarM1000_frac_p_DRP8_hist = new TH1F("Matched_perm_WJb_ttbarM1000_frac_p_DRP8", "", 10, 0., 5.);

			//Matched objects
	    	Matched_BHadWJa_ptthad_3J_DRP8_hist = new TH1F("Matched_BHadWJa_ptthad_3J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_3J_DRP8_hist = new TH1F("Matched_BHadWJb_ptthad_3J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_3J_DRP8_hist = new TH1F("Matched_WJaWJb_ptthad_3J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_3J_DRP8_hist = new TH1F("All_Matched_BHadWJa_ptthad_3J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_3J_DRP8_hist = new TH1F("All_Matched_BHadWJb_ptthad_3J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_3J_DRP8_hist = new TH1F("All_Matched_WJaWJb_ptthad_3J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_3J_DRP8_hist->Sumw2();	
		Matched_BHadWJb_ptthad_3J_DRP8_hist->Sumw2();	
		Matched_WJaWJb_ptthad_3J_DRP8_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_3J_DRP8_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_3J_DRP8_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_3J_DRP8_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_4J_DRP8_hist = new TH1F("Matched_BHadWJa_ptthad_4J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_4J_DRP8_hist = new TH1F("Matched_BHadWJb_ptthad_4J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_4J_DRP8_hist = new TH1F("Matched_WJaWJb_ptthad_4J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_4J_DRP8_hist = new TH1F("All_Matched_BHadWJa_ptthad_4J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_4J_DRP8_hist = new TH1F("All_Matched_BHadWJb_ptthad_4J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_4J_DRP8_hist = new TH1F("All_Matched_WJaWJb_ptthad_4J_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_4J_DRP8_hist->Sumw2();	
		Matched_BHadWJb_ptthad_4J_DRP8_hist->Sumw2();	
		Matched_WJaWJb_ptthad_4J_DRP8_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_4J_DRP8_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_4J_DRP8_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_4J_DRP8_hist->Sumw2();	

	    	Matched_BHadWJa_ptthad_5PJ_DRP8_hist = new TH1F("Matched_BHadWJa_ptthad_5PJ_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_BHadWJb_ptthad_5PJ_DRP8_hist = new TH1F("Matched_BHadWJb_ptthad_5PJ_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	Matched_WJaWJb_ptthad_5PJ_DRP8_hist = new TH1F("Matched_WJaWJb_ptthad_5PJ_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJa_ptthad_5PJ_DRP8_hist = new TH1F("All_Matched_BHadWJa_ptthad_5PJ_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_BHadWJb_ptthad_5PJ_DRP8_hist = new TH1F("All_Matched_BHadWJb_ptthad_5PJ_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
	    	All_Matched_WJaWJb_ptthad_5PJ_DRP8_hist = new TH1F("All_Matched_WJaWJb_ptthad_5PJ_DRP8_hist", "", pt_bins_, pt_min_, pt_max_);
		Matched_BHadWJa_ptthad_5PJ_DRP8_hist->Sumw2();	
		Matched_BHadWJb_ptthad_5PJ_DRP8_hist->Sumw2();	
		Matched_WJaWJb_ptthad_5PJ_DRP8_hist->Sumw2();	
		All_Matched_BHadWJa_ptthad_5PJ_DRP8_hist->Sumw2();	
		All_Matched_BHadWJb_ptthad_5PJ_DRP8_hist->Sumw2();	
		All_Matched_WJaWJb_ptthad_5PJ_DRP8_hist->Sumw2();	

		BTag_lp8_BHad_loose_pass_hist = new TH1F("BTag_lp8_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BHad_loose_fail_hist = new TH1F("BTag_lp8_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BHad_medium_pass_hist = new TH1F("BTag_lp8_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BHad_medium_fail_hist = new TH1F("BTag_lp8_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BHad_tight_pass_hist = new TH1F("BTag_lp8_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BHad_tight_fail_hist = new TH1F("BTag_lp8_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BHad_loose_pass_hist->Sumw2(); 
        	BTag_lp8_BHad_loose_fail_hist->Sumw2();
        	BTag_lp8_BHad_medium_pass_hist->Sumw2();
        	BTag_lp8_BHad_medium_fail_hist->Sumw2();
        	BTag_lp8_BHad_tight_pass_hist->Sumw2(); 
		BTag_lp8_BHad_tight_fail_hist->Sumw2(); 

		BTag_lp8_BLep_loose_pass_hist = new TH1F("BTag_lp8_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BLep_loose_fail_hist = new TH1F("BTag_lp8_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BLep_medium_pass_hist = new TH1F("BTag_lp8_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BLep_medium_fail_hist = new TH1F("BTag_lp8_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BLep_tight_pass_hist = new TH1F("BTag_lp8_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BLep_tight_fail_hist = new TH1F("BTag_lp8_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_lp8_BLep_loose_pass_hist->Sumw2(); 
        	BTag_lp8_BLep_loose_fail_hist->Sumw2();
        	BTag_lp8_BLep_medium_pass_hist->Sumw2();
        	BTag_lp8_BLep_medium_fail_hist->Sumw2();
        	BTag_lp8_BLep_tight_pass_hist->Sumw2(); 
		BTag_lp8_BLep_tight_fail_hist->Sumw2(); 

		BTag_gp8_BHad_loose_pass_hist = new TH1F("BTag_gp8_BHad_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BHad_loose_fail_hist = new TH1F("BTag_gp8_BHad_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BHad_medium_pass_hist = new TH1F("BTag_gp8_BHad_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BHad_medium_fail_hist = new TH1F("BTag_gp8_BHad_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BHad_tight_pass_hist = new TH1F("BTag_gp8_BHad_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BHad_tight_fail_hist = new TH1F("BTag_gp8_BHad_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BHad_loose_pass_hist->Sumw2(); 
        	BTag_gp8_BHad_loose_fail_hist->Sumw2();
        	BTag_gp8_BHad_medium_pass_hist->Sumw2();
        	BTag_gp8_BHad_medium_fail_hist->Sumw2();
        	BTag_gp8_BHad_tight_pass_hist->Sumw2(); 
		BTag_gp8_BHad_tight_fail_hist->Sumw2(); 

		BTag_gp8_BLep_loose_pass_hist = new TH1F("BTag_gp8_BLep_loose_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BLep_loose_fail_hist = new TH1F("BTag_gp8_BLep_loose_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BLep_medium_pass_hist = new TH1F("BTag_gp8_BLep_medium_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BLep_medium_fail_hist = new TH1F("BTag_gp8_BLep_medium_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BLep_tight_pass_hist = new TH1F("BTag_gp8_BLep_tight_pass", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BLep_tight_fail_hist = new TH1F("BTag_gp8_BLep_tight_fail", "", nbins, mass_min_, mass_max_);
		BTag_gp8_BLep_loose_pass_hist->Sumw2(); 
        	BTag_gp8_BLep_loose_fail_hist->Sumw2();
        	BTag_gp8_BLep_medium_pass_hist->Sumw2();
        	BTag_gp8_BLep_medium_fail_hist->Sumw2();
        	BTag_gp8_BLep_tight_pass_hist->Sumw2(); 
		BTag_gp8_BLep_tight_fail_hist->Sumw2(); 

			//Merged jets
		Merged_BHadWJa_perm_and_WJb_mass_DRP8_hist = new TH1F("Merged_BHadWJa_perm_and_WJb_mass_DRP8", "", 100, 0., 500.);
		Merged_BHadWJa_perm_mass_DRP8_hist = new TH1F("Merged_BHadWJa_perm_mass_DRP8", "", 100, 0., 500.);
		Merged_WJb_mass_DRP8_hist = new TH1F("Merged_WJb_mass_DRP8", "", 100, 0., 300.);
		Merged_BHadWJa_perm_DRP8_LepBLep_hist = new TH1F("Merged_BHadWJa_perm_DRP8_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJa_massDivpt_DRP8_hist = new TH1F("Merged_BHadWJa_massDivpt_DRP8", "", 50, 0., 0.5);

		Merged_BHadWJb_perm_and_WJa_mass_DRP8_hist = new TH1F("Merged_BHadWJb_perm_and_WJa_mass_DRP8", "", 100, 0., 500.);
		Merged_BHadWJb_perm_mass_DRP8_hist = new TH1F("Merged_BHadWJb_perm_mass_DRP8", "", 100, 0., 500.);
		Merged_WJa_mass_DRP8_hist = new TH1F("Merged_WJa_mass_DRP8", "", 100, 0., 300.);
		Merged_BHadWJb_perm_DRP8_LepBLep_hist = new TH1F("Merged_BHadWJb_perm_DRP8_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BHadWJb_massDivpt_DRP8_hist = new TH1F("Merged_BHadWJb_massDivpt_DRP8", "", 50, 0., 0.5);

	        Merged_BHadWJet_mass_vs_BHad_mass_DRP8_hist = new TH2D("Merged_BHadWJet_mass_vs_BHad_mass_DRP8","",100, 0.,300.,100,0., 300.);

		Merged_BLepWJa_perm_and_WJb_mass_DRP8_hist = new TH1F("Merged_BLepWJa_perm_and_WJb_mass_DRP8", "", 100, 0., 300.);
		Merged_BLepWJa_perm_mass_DRP8_hist = new TH1F("Merged_BLepWJa_perm_mass_DRP8", "", 100, 0., 300.);
		Merged_BLepWJa_perm_DRP8_LepBLep_hist = new TH1F("Merged_BLepWJa_perm_DRP8_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJa_massDivpt_DRP8_hist = new TH1F("Merged_BLepWJa_massDivpt_DRP8", "", 50, 0., 0.5);

		Merged_BLepWJb_perm_and_WJa_mass_DRP8_hist = new TH1F("Merged_BLepWJb_perm_and_WJa_mass_DRP8", "", 100, 0., 300.);
		Merged_BLepWJb_perm_mass_DRP8_hist = new TH1F("Merged_BLepWJb_perm_mass_DRP8", "", 100, 0., 300.);
		Merged_BLepWJb_perm_DRP8_LepBLep_hist = new TH1F("Merged_BLepWJb_perm_DRP8_LepBLep", "", DR_bins_, DR_min_, DR_max_);
	        Merged_BLepWJb_massDivpt_DRP8_hist = new TH1F("Merged_BLepWJb_massDivpt_DRP8", "", 50, 0., 0.5);

	        Merged_BLepWJet_mass_vs_BLep_mass_DRP8_hist = new TH2D("Merged_BLepWJet_mass_vs_BLep_mass_DRP8","",100, 0.,300.,100,0., 300.);

			//Unmerged jets
		Unmerged_BHad_mass_DRP8_hist = new TH1F("Unmerged_BHad_mass_DRP8", "", 50, 0., 100.);
		Unmerged_WJa_mass_DRP8_hist = new TH1F("Unmerged_WJa_mass_DRP8", "", 50, 0., 100.);
		Unmerged_WJb_mass_DRP8_hist = new TH1F("Unmerged_WJb_mass_DRP8", "", 50, 0., 100.);
		Unmerged_BHad_LepBLep_DRP8_hist = new TH1F("Unmerged_BHad_LepBLep_DRP8", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJa_LepBLep_DRP8_hist = new TH1F("Unmerged_WJa_LepBLep_DRP8", "", DR_bins_, DR_min_, DR_max_);
		Unmerged_WJb_LepBLep_DRP8_hist = new TH1F("Unmerged_WJb_LepBLep_DRP8", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BHad_massDivpt_DRP8_hist = new TH1F("Unmerged_BHad_massDivpt_DRP8", "", 50, 0., 0.5);
	        Unmerged_WJa_massDivpt_DRP8_hist = new TH1F("Unmerged_WJa_massDivpt_DRP8", "", 50, 0., 0.5);
	        Unmerged_WJb_massDivpt_DRP8_hist = new TH1F("Unmerged_WJb_massDivpt_DRP8", "", 50, 0., 0.5);

	        Unmerged_BHadWJet_mass_vs_BHad_mass_DRP8_hist = new TH2D("Unmerged_BHadWJet_mass_vs_BHad_mass_DRP8","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP8_hist = new TH2D("Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP8", "", 100, 0., 300., 100, 0., 300.);

		Unmerged_BLep_mass_DRP8_hist = new TH1F("Unmerged_BLep_mass_DRP8", "", 50, 0., 100.);
		Unmerged_BLep_LepBLep_DRP8_hist = new TH1F("Unmerged_BLep_LepBLep_DRP8", "", DR_bins_, DR_min_, DR_max_);
	        Unmerged_BLep_massDivpt_DRP8_hist = new TH1F("Unmerged_BLep_massDivpt_DRP8", "", 50, 0., 0.5);

	        Unmerged_BLepWJet_mass_vs_BLep_mass_DRP8_hist = new TH2D("Unmerged_BLepWJet_mass_vs_BLep_mass_DRP8","",100, 0., 300., 100, 0., 300.);
	        Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP8_hist = new TH2D("Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP8", "", 100, 0., 300., 100, 0., 300.);


		Logger::log().debug() << "End of begin() " << evt_idx_ << endl;
        }

        //This method is called once every file, contains the event loop
        ///run your proper analysis here
	virtual void analyze()
	{
		Logger::log().debug() << "Beginning of analyze() " << evt_idx_ << endl;

                URStreamer event(tree_);

		while(event.next() /*&& evt_idx_ < 30*/)
                {
                        evt_idx_++;
                        if(evt_idx_ % 10000 == 0) Logger::log().debug() << "Beginning event " << evt_idx_ << endl;

			//generator selection
			bool selection = genp_selector_.select(event);
                        tracker_.track("gen selection");
                        if( !selection ){
                                Logger::log().debug() << "event has no selection " << endl;
                                continue;
                        }
			GenTTBar &ttbar = genp_selector_.ttbar_system();
			Pt_ttbar_hist->Fill(ttbar.Pt());
                        Mass_ttbar_hist->Fill(ttbar.M());

			//Gen level matching
			Permutation matched_perm;
	    		double nmatched_objects = 0;

			if( object_selector_.select(event) ){
				nJets_hist->Fill(object_selector_.clean_jets().size());
				if( ttbar.type == GenTTBar::DecayType::SEMILEP ){
					matched_perm = matcher_.match(
						genp_selector_.ttbar_final_system(),
						object_selector_.clean_jets(),
						object_selector_.loose_electrons(),
						object_selector_.loose_muons());
				}
				matched_perm.SetMET(object_selector_.met());
//				if( object_selector_.clean_jets().size() == 3 ){
//					for(vector<IDJet*>::iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){ 
////						cout << "jet E: " << (*jets)->E() << endl;
////						cout << "jet |p|: " << (*jets)->P() << endl;
////						cout << "jet M: " << (*jets)->M() << endl;
//					}
//				}
			}
//			if( matched_perm.L() || matched_perm.BHad() || matched_perm.BLep() || matched_perm.WJa() || matched_perm.WJb() ) cout << matched_perm << endl;
//			cout << matched_perm << endl;

			//initialize kinematic vars
			double DR_LepBHad = -1;
			double DR_LepBLep = -1;
			double DR_LepWJa = -1;
			double DR_LepWJb = -1;
			double DR_BHadBLep = -1;
			double DR_BHadWJa = -1;
			double DR_BHadWJb = -1;
			double DR_BLepWJa = -1;
			double DR_BLepWJb = -1;
			double DR_WJaWJb = -1;
			
			double Pt_Lep = -1;
			double Pt_BHad = -1;
			double Pt_BLep = -1;
			double Pt_WJa = -1;
			double Pt_WJb = -1;
			
			double Eta_Lep = 1e3;
			double Eta_BHad = 1e3;
			double Eta_BLep = 1e3;
			double Eta_WJa = 1e3;
			double Eta_WJb = 1e3;

			//initialize gen jets
			const GenObject* lepton = 0;
			const GenObject* BLep = 0;
			const GenObject* BHad = 0;
			const GenObject* WJa = 0;
			const GenObject* WJb = 0;
			const GenTop* tlep = 0;
			const GenTop* thad = 0;

			//initialize permutation objects
				//DR=0.4
			const IDJet* best_BHad_DRP4 = 0;
			const IDJet* best_BLep_DRP4 = 0;
			const IDJet* best_WJa_DRP4 = 0;
			const IDJet* best_WJb_DRP4 = 0;

				//DR=0.5
			const IDJet* best_BHad_DRP5 = 0;
			const IDJet* best_BLep_DRP5 = 0;
			const IDJet* best_WJa_DRP5 = 0;
			const IDJet* best_WJb_DRP5 = 0;

				//DR=0.6
			const IDJet* best_BHad_DRP6 = 0;
			const IDJet* best_BLep_DRP6 = 0;
			const IDJet* best_WJa_DRP6 = 0;
			const IDJet* best_WJb_DRP6 = 0;

				//DR=0.8
			const IDJet* best_BHad_DRP8 = 0;
			const IDJet* best_BLep_DRP8 = 0;
			const IDJet* best_WJa_DRP8 = 0;
			const IDJet* best_WJb_DRP8 = 0;

			//Generator level definitions
			if( ttbar.lep_top() ){ //check if lep top exists
				tlep = ttbar.lep_top();
				Pt_tlep_hist->Fill(tlep->Pt());
				Mass_tlep_hist->Fill(tlep->M());
			}
			if( ttbar.had_top() ){ //check if had top exists
				thad = ttbar.had_top();
				Pt_thad_hist->Fill(thad->Pt());
				Mass_thad_hist->Fill(thad->M());
			}
			if( ttbar.lepton() ){
				lepton = ttbar.lepton();
				Pt_Lep = lepton->Pt();
                        	Pt_Lep_vs_Mtt_hist->Fill(ttbar.M(), Pt_Lep);
                        	Pt_Lep_hist->Fill(Pt_Lep);

                        	Eta_Lep = lepton->Eta();
				Eta_Lep_vs_Mtt_hist->Fill(ttbar.M(), Eta_Lep);
                        	Eta_Lep_hist->Fill(Eta_Lep);
			}
			if( ttbar.lep_b() ){
				BLep = ttbar.lep_b();
				Pt_BLep = BLep->Pt();
                       		Pt_BLep_vs_Mtt_hist->Fill(ttbar.M(), Pt_BLep);
                        	Pt_BLep_hist->Fill(Pt_BLep);

                        	Eta_BLep = BLep->Eta();
                        	Eta_BLep_vs_Mtt_hist->Fill(ttbar.M(), Eta_BLep);
                        	Eta_BLep_hist->Fill(Eta_BLep);
			}
			if( ttbar.had_b() ){
				BHad = ttbar.had_b();
				Pt_BHad = BHad->Pt();
                        	Pt_BHad_vs_Mtt_hist->Fill(ttbar.M(), Pt_BHad);
                        	Pt_BHad_hist->Fill(Pt_BHad);

                        	Eta_BHad = BHad->Eta();
                        	Eta_BHad_vs_Mtt_hist->Fill(ttbar.M(), Eta_BHad);
                        	Eta_BHad_hist->Fill(Eta_BHad);
//				cout << "BHad Energy, |p|: " << BHad->E() << ", " << sqrt(pow(BHad->Px(), 2.0)+pow(BHad->Py(),2.0)+pow(BHad->Pz(),2.0)) << endl;
			}
			if( ttbar.had_W() ){
				if( ttbar.had_W()->first && !ttbar.had_W()->second ){
					WJa = ttbar.had_W()->first;
					Pt_WJa = WJa->Pt();
                        		Pt_WJa_vs_Mtt_hist->Fill(ttbar.M(), Pt_WJa);
                        		Pt_WJa_hist->Fill(Pt_WJa);

                        		Eta_WJa = WJa->Eta();
                        		Eta_WJa_vs_Mtt_hist->Fill(ttbar.M(), Eta_WJa);
                        		Eta_WJa_hist->Fill(Eta_WJa);
				}
				if( ttbar.had_W()->second && !ttbar.had_W()->first ){
					WJb = ttbar.had_W()->second;
					Pt_WJb = WJb->Pt();
                        		Pt_WJb_vs_Mtt_hist->Fill(ttbar.M(), Pt_WJb);
                        		Pt_WJb_hist->Fill(Pt_WJb);

                        		Eta_WJb = WJb->Eta();
                        		Eta_WJb_vs_Mtt_hist->Fill(ttbar.M(), Eta_WJb);
                        		Eta_WJb_hist->Fill(Eta_WJb);
				}
				if( ttbar.had_W()->first && ttbar.had_W()->second ){
					WJa = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->first : ttbar.had_W()->second;
					WJb = (ttbar.had_W()->first->E() > ttbar.had_W()->second->E() ) ? ttbar.had_W()->second : ttbar.had_W()->first;

					Pt_WJa = WJa->Pt();
                        		Pt_WJa_vs_Mtt_hist->Fill(ttbar.M(), Pt_WJa);
                        		Pt_WJa_hist->Fill(Pt_WJa);

                        		Eta_WJa = WJa->Eta();
                        		Eta_WJa_vs_Mtt_hist->Fill(ttbar.M(), Eta_WJa);
                        		Eta_WJa_hist->Fill(Eta_WJa);

					Pt_WJb = WJb->Pt();
                        		Pt_WJb_vs_Mtt_hist->Fill(ttbar.M(), Pt_WJb);
                        		Pt_WJb_hist->Fill(Pt_WJb);

                        		Eta_WJb = WJb->Eta();
                        		Eta_WJb_vs_Mtt_hist->Fill(ttbar.M(), Eta_WJb);
                        		Eta_WJb_hist->Fill(Eta_WJb);
				}
			}

				//DR
			if( !(lepton == 0) && !(BHad == 0) ){
				DR_LepBHad = lepton->DeltaR(*BHad);
				DR_LepBHad_vs_Mtt_hist->Fill(ttbar.M(), DR_LepBHad);
				DR_LepBHad_hist->Fill(DR_LepBHad);
				if( DR_LepBHad < 0.4 ){
				        DR_LepBHad_lp4_vs_mttbar_hist->Fill(ttbar.M());
				}
				if( DR_LepBHad >= 0.4 ){
				        DR_LepBHad_gp4_vs_mttbar_hist->Fill(ttbar.M());
				}
			}
                        if( !(lepton == 0) && !(BLep == 0) ){ //check if lepton and lep b exist
			        DR_LepBLep = lepton->DeltaR(*BLep);
			        DR_LepBLep_vs_Mtt_hist->Fill(ttbar.M(), DR_LepBLep);
			        DR_LepBLep_hist->Fill(DR_LepBLep);
			        DRmin_tlep_vs_mttbar_hist->Fill(ttbar.M(), DR_LepBLep);
			        DRmin_tlep_vs_pttlep_hist->Fill(tlep->Pt(), DR_LepBLep);
			        DRmin_tlep_hist->Fill(DR_LepBLep);
			        if( DR_LepBLep < 0.4 ){
			                DRmin_tlep_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
			                DRmin_tlep_lp4_vs_pttlep_hist->Fill(tlep->Pt());
			                DR_LepBLep_lp4_vs_mttbar_hist->Fill(ttbar.M());
			        }
			        if( DR_LepBLep >= 0.4 ){
			                DRmin_tlep_gp4_vs_mttbar_hist->Fill(ttbar.M());
			                DRmin_tlep_gp4_vs_pttlep_hist->Fill(tlep->Pt());
			                DR_LepBLep_gp4_vs_mttbar_hist->Fill(ttbar.M());
			        }
			}
			if( !(lepton == 0) && !(WJa == 0) ){ //check if lepton and WJa exist
			        DR_LepWJa = lepton->DeltaR(*WJa);
			        DR_LepWJa_vs_Mtt_hist->Fill(ttbar.M(), DR_LepWJa);
			        DR_LepWJa_hist->Fill(DR_LepWJa);
			        if( DR_LepWJa < 0.4 ){
			                DR_LepWJa_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			        if( DR_LepWJa >= 0.4 ){
			                DR_LepWJa_gp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			}
			  if( !(lepton == 0) && !(WJb == 0) ){
			        DR_LepWJb = lepton->DeltaR(*WJb);
			        DR_LepWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_LepWJb);
			        DR_LepWJb_hist->Fill(DR_LepWJb);
			        if( DR_LepWJb < 0.4 ){
			                DR_LepWJb_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			        if( DR_LepWJb >= 0.4 ){
			                DR_LepWJb_gp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			}
                        if( !(BHad == 0) && !(BLep == 0) ){
			        DR_BHadBLep = BHad->DeltaR(*BLep);
			        DR_BHadBLep_vs_Mtt_hist->Fill(ttbar.M(), DR_BHadBLep);
			        DR_BHadBLep_hist->Fill(DR_BHadBLep);
			        if( DR_BHadBLep < 0.4 ){
			                DR_BHadBLep_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			        if( DR_BHadBLep >= 0.4 ){
			                DR_BHadBLep_gp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			}
                        if( !(BHad == 0) && !(WJa == 0) ){
				DR_BHadWJa = BHad->DeltaR(*WJa);
				DR_BHadWJa_vs_Mtt_hist->Fill(ttbar.M(), DR_BHadWJa);
				DR_BHadWJa_hist->Fill(DR_BHadWJa);
				if( DR_BHadWJa < 0.4 ){
				        DR_BHadWJa_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
                        	}
				if( DR_BHadWJa >= 0.4 ){
				        DR_BHadWJa_gp4_vs_mttbar_hist->Fill(ttbar.M());
				}
			}
                        if( !(BHad == 0) && !(WJb == 0) ){
				DR_BHadWJb = BHad->DeltaR(*WJb);
				DR_BHadWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_BHadWJb);
				DR_BHadWJb_hist->Fill(DR_BHadWJb);
				if( DR_BHadWJb < 0.4 ){
				        DR_BHadWJb_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
                        	}
				if( DR_BHadWJb >= 0.4 ){
				        DR_BHadWJb_gp4_vs_mttbar_hist->Fill(ttbar.M());
				}
			}
                        if( !(BLep == 0) && !(WJa == 0) ){
			        DR_BLepWJa = BLep->DeltaR(*WJa);
			        DR_BLepWJa_vs_Mtt_hist->Fill(ttbar.M(), DR_BLepWJa);
			        DR_BLepWJa_hist->Fill(DR_BLepWJa);
			        if( DR_BLepWJa < 0.4 ){
			                DR_BLepWJa_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			        if( DR_BLepWJa >= 0.4 ){
			                DR_BLepWJa_gp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			}
			if( !(BLep == 0) && !(WJb == 0) ){
			        DR_BLepWJb = BLep->DeltaR(*WJb);
			        DR_BLepWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_BLepWJb);
			        DR_BLepWJb_hist->Fill(DR_BLepWJb);
			        if( DR_BLepWJb < 0.4 ){
			                DR_BLepWJb_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			        if( DR_BLepWJb >= 0.4 ){
			                DR_BLepWJb_gp4_vs_mttbar_hist->Fill(ttbar.M()); 
			        }
			}
                        if( !(WJa == 0) && !(WJb == 0) ){
			        DR_WJaWJb = WJa->DeltaR(*WJb);
			        DR_WJaWJb_vs_Mtt_hist->Fill(ttbar.M(), DR_WJaWJb);
			        DR_WJaWJb_hist->Fill(DR_WJaWJb);
			        if( DR_WJaWJb < 0.4 ){
			                DR_WJaWJb_lp4_vs_mttbar_hist->Fill(ttbar.M());   
			        }
			        if( DR_WJaWJb >= 0.4 ){
			                DR_WJaWJb_gp4_vs_mttbar_hist->Fill(ttbar.M());
			        }
			}

			//find min value of DR for thad objects
			list<double> thad_DR_list;
			list<double>::iterator thad_DR_it;
			if( !(BHad == 0) && !(WJa == 0) ) thad_DR_list.push_back(DR_BHadWJa);
			if( !(BHad == 0) && !(WJb == 0) ) thad_DR_list.push_back(DR_BHadWJb);
			if( !(WJa == 0) && !(WJb == 0) ) thad_DR_list.push_back(DR_WJaWJb);
			double thad_DRmin = 1e2;
			double thad_DRmax = 0;
			for( thad_DR_it = thad_DR_list.begin(); thad_DR_it != thad_DR_list.end(); ++thad_DR_it ){
				if( *thad_DR_it < thad_DRmin ){
					thad_DRmin = *thad_DR_it;
				}
				if( *thad_DR_it > thad_DRmax ){
					thad_DRmax = *thad_DR_it;
				}
			}
			if( thad_DRmin < 0.4 ){
				for( thad_DR_it = thad_DR_list.begin(); thad_DR_it != thad_DR_list.end(); ++thad_DR_it ){
					if( *thad_DR_it > thad_DRmax ){
						thad_DRmax = *thad_DR_it;
					}
				}

			}
			if( !(thad_DRmax == 0) ){
				DRmax_thad_vs_mttbar_hist->Fill(ttbar.M(), thad_DRmax);
				DRmax_thad_vs_ptthad_hist->Fill(thad->Pt(), thad_DRmax);
				DRmax_thad_hist->Fill(thad_DRmax);

				if( thad_DRmax < 0.4 ){
					DRmax_thad_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
					DRmax_thad_lp4_vs_ptthad_hist->Fill(thad->Pt());
				}
				if( thad_DRmax >= 0.4 ){
					DRmax_thad_gp4_vs_mttbar_hist->Fill(ttbar.M());
					DRmax_thad_gp4_vs_ptthad_hist->Fill(thad->Pt());
				}
			}
			if( !(thad_DRmin == 1e2) ){
				DRmin_thad_vs_mttbar_hist->Fill(ttbar.M(), thad_DRmin);
				DRmin_thad_vs_ptthad_hist->Fill(thad->Pt(), thad_DRmin);
				DRmin_thad_hist->Fill(thad_DRmin);

				if( thad_DRmin < 0.4 ){
					DRmin_thad_lp4_vs_mttbar_hist->Fill(ttbar.M()); 
					DRmin_thad_lp4_vs_ptthad_hist->Fill(thad->Pt());
				}
				if( thad_DRmin >= 0.4 ){
					DRmin_thad_gp4_vs_mttbar_hist->Fill(ttbar.M());
					DRmin_thad_gp4_vs_ptthad_hist->Fill(thad->Pt());
				}
			}


			//Perm Matching
				//DR=0.4
			list<IDJet*> bhad_DRP4_list;
			float bhad_DRP4 = 1e10;// itialize dr to high number
			if( object_selector_.clean_jets().size() < 3 ) continue;
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BHad == 0) ){
					if( (*jets)->DeltaR(*BHad) > 0.4 ) continue;
					if( (*jets)->DeltaR(*BHad) < bhad_DRP4 ){
						bhad_DRP4 = (*jets)->DeltaR(*BHad);
						bhad_DRP4_list.push_back(*jets);
						best_BHad_DRP4 = (bhad_DRP4_list.back());
					}
				}
			}
			list<IDJet*> blep_DRP4_list;
			float blep_DRP4 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BLep == 0) ){
					if( (*jets)->DeltaR(*BLep) > 0.4 ) continue;
					if( (*jets)->DeltaR(*BLep) < blep_DRP4 ){
						blep_DRP4 = (*jets)->DeltaR(*BLep);
						blep_DRP4_list.push_back(*jets);
						best_BLep_DRP4 = (blep_DRP4_list.back());
					}
				}
			}
			list<IDJet*> wja_DRP4_list;
			float wja_DRP4 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJa == 0) ){
					if( (*jets)->DeltaR(*WJa) > 0.4 ) continue;
					if( (*jets)->DeltaR(*WJa) < wja_DRP4 ){
						wja_DRP4 = (*jets)->DeltaR(*WJa);
						wja_DRP4_list.push_back(*jets);
						best_WJa_DRP4 = (wja_DRP4_list.back());
					}
				}
			}
			list<IDJet*> wjb_DRP4_list;
			float wjb_DRP4 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJb == 0) ){
					if( (*jets)->DeltaR(*WJb) > 0.4 ) continue;
					if( (*jets)->DeltaR(*WJb) < wjb_DRP4 ){
						wjb_DRP4 = (*jets)->DeltaR(*WJb);
						wjb_DRP4_list.push_back(*jets);
						best_WJb_DRP4 = (wjb_DRP4_list.back());
					}
				}
			}

			//Reco Obj. Kinematic Vars
			if( !(best_BHad_DRP4 == 0) ){//reco BHad exists
				Matched_perm_BHad_pt_DRP4_hist->Fill(best_BHad_DRP4->Pt());
				Matched_perm_BHad_eta_DRP4_hist->Fill(best_BHad_DRP4->Eta());
//						++nmatched_objects;
				if( best_BHad_DRP4 == best_BLep_DRP4 || best_BHad_DRP4 == best_WJa_DRP4 || best_BHad_DRP4 == best_WJb_DRP4 ){//reco BHad merged with something
					if( best_BHad_DRP4->BTagId(cut_loose_b_) ) BTag_lp4_BHad_loose_pass_hist->Fill(ttbar.M());
					if( best_BHad_DRP4->BTagId(cut_medium_b_) ) BTag_lp4_BHad_medium_pass_hist->Fill(ttbar.M());	
					if( best_BHad_DRP4->BTagId(cut_tight_b_) ) BTag_lp4_BHad_tight_pass_hist->Fill(ttbar.M());
					if( !best_BHad_DRP4->BTagId(cut_loose_b_) ) BTag_lp4_BHad_loose_fail_hist->Fill(ttbar.M());
					if( !best_BHad_DRP4->BTagId(cut_medium_b_) ) BTag_lp4_BHad_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BHad_DRP4->BTagId(cut_tight_b_) ) BTag_lp4_BHad_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_BLep_DRP4 == 0) ){//reco BLep exists
				Matched_perm_BLep_pt_DRP4_hist->Fill(best_BLep_DRP4->Pt());
				Matched_perm_BLep_eta_DRP4_hist->Fill(best_BLep_DRP4->Eta());
//				++nmatched_objects; 
				if( best_BLep_DRP4 == best_BHad_DRP4 || best_BLep_DRP4 == best_WJa_DRP4 || best_BLep_DRP4 == best_WJb_DRP4 ){//reco BLep merged with something
					if( best_BLep_DRP4->BTagId(cut_loose_b_) ) BTag_lp4_BLep_loose_pass_hist->Fill(ttbar.M());
					if( best_BLep_DRP4->BTagId(cut_medium_b_) ) BTag_lp4_BLep_medium_pass_hist->Fill(ttbar.M());	
					if( best_BLep_DRP4->BTagId(cut_tight_b_) ) BTag_lp4_BLep_tight_pass_hist->Fill(ttbar.M());
					if( !best_BLep_DRP4->BTagId(cut_loose_b_) ) BTag_lp4_BLep_loose_fail_hist->Fill(ttbar.M());
					if( !best_BLep_DRP4->BTagId(cut_medium_b_) ) BTag_lp4_BLep_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BLep_DRP4->BTagId(cut_tight_b_) ) BTag_lp4_BLep_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_WJa_DRP4 == 0) ){//reco WJa exists
//						++nmatched_objects; 
				Matched_perm_WJa_pt_DRP4_hist->Fill(best_WJa_DRP4->Pt());
				Matched_perm_WJa_eta_DRP4_hist->Fill(best_WJa_DRP4->Eta());
                        	if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                        	        Matched_perm_WJa_ttbarM700_pt_DRP4_hist->Fill(best_WJa_DRP4->Pt());
                        	        Matched_perm_WJa_ttbarM700_frac_p_DRP4_hist->Fill(best_WJa_DRP4->P()/(ttbar.M()/2));
                        	}
                        	if( ttbar.M() >= 1000 ){
                        	        Matched_perm_WJa_ttbarM1000_pt_DRP4_hist->Fill(best_WJa_DRP4->Pt());
                        	        Matched_perm_WJa_ttbarM1000_frac_p_DRP4_hist->Fill(best_WJa_DRP4->P()/(ttbar.M()/2));
                        	}
			}
			if( !(best_WJb_DRP4 == 0) ){//reco WJb exists
//						++nmatched_objects; 
				Matched_perm_WJb_pt_DRP4_hist->Fill(best_WJb_DRP4->Pt());
				Matched_perm_WJb_eta_DRP4_hist->Fill(best_WJb_DRP4->Eta());
				if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
				        Matched_perm_WJb_ttbarM700_pt_DRP4_hist->Fill(best_WJb_DRP4->Pt());
				        Matched_perm_WJb_ttbarM700_frac_p_DRP4_hist->Fill(best_WJb_DRP4->P()/(ttbar.M()/2));
				}
				if( ttbar.M() >= 1000 ){
				        Matched_perm_WJb_ttbarM1000_pt_DRP4_hist->Fill(best_WJb_DRP4->Pt());
				        Matched_perm_WJb_ttbarM1000_frac_p_DRP4_hist->Fill(best_WJb_DRP4->P()/(ttbar.M()/2));
				}
			}
			if( !(best_BLep_DRP4 == 0) && !(best_BLep_DRP4 == best_BHad_DRP4) && !(best_BLep_DRP4 == best_WJa_DRP4) && !(best_BLep_DRP4 == best_WJb_DRP4) ){//reco BLep not merged with anything
				if( best_BLep_DRP4->BTagId(cut_loose_b_) ) BTag_gp4_BLep_loose_pass_hist->Fill(ttbar.M());
				if( best_BLep_DRP4->BTagId(cut_medium_b_) ) BTag_gp4_BLep_medium_pass_hist->Fill(ttbar.M());	
				if( best_BLep_DRP4->BTagId(cut_tight_b_) ) BTag_gp4_BLep_tight_pass_hist->Fill(ttbar.M());
				if( !best_BLep_DRP4->BTagId(cut_loose_b_) ) BTag_gp4_BLep_loose_fail_hist->Fill(ttbar.M());
				if( !best_BLep_DRP4->BTagId(cut_medium_b_) ) BTag_gp4_BLep_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BLep_DRP4->BTagId(cut_tight_b_) ) BTag_gp4_BLep_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BLep_LepBLep_DRP4_hist->Fill(DR_LepBLep);
	
				Unmerged_BLep_massDivpt_DRP4_hist->Fill(best_BLep_DRP4->M()/best_BLep_DRP4->Pt());
				if( !(best_WJa_DRP4 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M(),(*best_BLep_DRP4+*best_WJa_DRP4).M());
				if( !(best_WJb_DRP4 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M(),(*best_BLep_DRP4+*best_WJb_DRP4).M());
				if( !(best_WJa_DRP4 == 0) && !(best_WJb_DRP4 == 0) ){
					if( best_WJa_DRP4->M() > best_WJb_DRP4->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M(),(*best_BLep_DRP4+*best_WJa_DRP4).M());
					if( best_WJa_DRP4->M() < best_WJb_DRP4->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M(),(*best_BLep_DRP4+*best_WJb_DRP4).M());
				}


			}
			if( !(best_BHad_DRP4 == 0) && !(best_BHad_DRP4 == best_WJa_DRP4) && !(best_BHad_DRP4 == best_WJb_DRP4) && !(best_BHad_DRP4 == best_BLep_DRP4) ){//reco BHad not merged with anything
				if( best_BHad_DRP4->BTagId(cut_loose_b_) ) BTag_gp4_BHad_loose_pass_hist->Fill(ttbar.M());
				if( best_BHad_DRP4->BTagId(cut_medium_b_) ) BTag_gp4_BHad_medium_pass_hist->Fill(ttbar.M());	
				if( best_BHad_DRP4->BTagId(cut_tight_b_) ) BTag_gp4_BHad_tight_pass_hist->Fill(ttbar.M());
				if( !best_BHad_DRP4->BTagId(cut_loose_b_) ) BTag_gp4_BHad_loose_fail_hist->Fill(ttbar.M());
				if( !best_BHad_DRP4->BTagId(cut_medium_b_) ) BTag_gp4_BHad_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BHad_DRP4->BTagId(cut_tight_b_) ) BTag_gp4_BHad_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BHad_LepBLep_DRP4_hist->Fill(DR_LepBLep);
	
				Unmerged_BHad_massDivpt_DRP4_hist->Fill(best_BHad_DRP4->M()/best_BHad_DRP4->Pt());
				if( !(best_WJa_DRP4 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M(),(*best_BHad_DRP4+*best_WJa_DRP4).M());
				if( !(best_WJb_DRP4 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M(),(*best_BHad_DRP4+*best_WJb_DRP4).M());
				if( !(best_WJa_DRP4 == 0) && !(best_WJb_DRP4 == 0) ){
					if( best_WJa_DRP4->M() > best_WJb_DRP4->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M(),(*best_BHad_DRP4+*best_WJa_DRP4).M());
					if( best_WJa_DRP4->M() < best_WJb_DRP4->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M(),(*best_BHad_DRP4+*best_WJb_DRP4).M());
				}
			}
			if( !(best_WJa_DRP4 == 0) && !(best_WJa_DRP4 == best_BHad_DRP4) && !(best_WJa_DRP4 == best_WJb_DRP4) && !(best_WJa_DRP4 == best_BLep_DRP4) ){//reco WJa not merged with anything
				Unmerged_WJa_massDivpt_DRP4_hist->Fill(best_WJa_DRP4->M()/best_WJa_DRP4->Pt());
				Unmerged_WJa_mass_DRP4_hist->Fill(best_WJa_DRP4->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJa_LepBLep_DRP4_hist->Fill(DR_LepBLep);
			}
			if( !(best_WJb_DRP4 == 0) && !(best_WJb_DRP4 == best_BHad_DRP4) && !(best_WJb_DRP4 == best_WJa_DRP4) && !(best_WJb_DRP4 == best_BLep_DRP4) ){//reco WJb not merged with anything
				Unmerged_WJb_massDivpt_DRP4_hist->Fill(best_WJb_DRP4->M()/best_WJb_DRP4->Pt());
				Unmerged_WJb_mass_DRP4_hist->Fill(best_WJb_DRP4->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJb_LepBLep_DRP4_hist->Fill(DR_LepBLep);
			}
			if( !(best_BHad_DRP4 == 0) && !(best_WJa_DRP4 == 0) ){ //reco BHad and WJa exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJa_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJa_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJa_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());

			}
			if( !(best_BHad_DRP4 == 0) && !(best_WJb_DRP4 == 0) ){ //reco BHad and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJb_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJb_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJb_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}
			if( !(best_WJa_DRP4 == 0) && !(best_WJb_DRP4 == 0) ){ //reco WJa and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_WJaWJb_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_WJaWJb_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_WJaWJb_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP4 == best_WJa_DRP4 && !(best_BHad_DRP4 == 0) && !(best_BHad_DRP4 == best_BLep_DRP4) && !(best_BHad_DRP4 == best_WJb_DRP4) ){//only reco BHad and WJa merged
				Merged_BHadWJa_massDivpt_DRP4_hist->Fill(best_BHad_DRP4->M()/best_BHad_DRP4->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJa_perm_DRP4_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJa_perm_mass_DRP4_hist->Fill(best_BHad_DRP4->M());//merged b and wjet mass
				if( !(best_WJb_DRP4 == 0) && !(best_WJb_DRP4 == best_BHad_DRP4) && !(best_WJb_DRP4 == best_BLep_DRP4) ){//reco WJb exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M(), (*best_BHad_DRP4+*best_WJb_DRP4).M());
					Merged_BHadWJa_perm_and_WJb_mass_DRP4_hist->Fill((*best_BHad_DRP4+*best_WJb_DRP4).M());//comb invariant mass of merged jet and other wjet
					Merged_WJb_mass_DRP4_hist->Fill(best_WJb_DRP4->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP4 == best_WJb_DRP4 && !(best_BHad_DRP4 == 0) && !(best_BHad_DRP4 == best_BLep_DRP4) && !(best_BHad_DRP4 == best_WJa_DRP4) ){//only reco BHad and WJb merged
				Merged_BHadWJb_massDivpt_DRP4_hist->Fill(best_BHad_DRP4->M()/best_BHad_DRP4->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJb_perm_DRP4_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJb_perm_mass_DRP4_hist->Fill(best_BHad_DRP4->M());//merged b and wjet mass
				if( !(best_WJa_DRP4 == 0) && !(best_WJa_DRP4 == best_BHad_DRP4) && !(best_WJa_DRP4 == best_BLep_DRP4) ){//reco WJa exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP4_hist->Fill(best_BHad_DRP4->M(), (*best_BHad_DRP4+*best_WJa_DRP4).M());
					Merged_BHadWJb_perm_and_WJa_mass_DRP4_hist->Fill((*best_BHad_DRP4+*best_WJa_DRP4).M());//comb mass of merged jet and other wjet
					Merged_WJa_mass_DRP4_hist->Fill(best_WJa_DRP4->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}
			if( best_WJa_DRP4 == best_WJb_DRP4 && !(best_WJa_DRP4 == 0) && !(best_WJa_DRP4 == best_BHad_DRP4) && !(best_WJa_DRP4 == best_BLep_DRP4) ){//only reco WJa and WJb merged
						if( object_selector_.clean_jets().size() == 3 ) Matched_WJaWJb_ptthad_3J_DRP4_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() == 4 ) Matched_WJaWJb_ptthad_4J_DRP4_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() > 4 ) Matched_WJaWJb_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}


			if( best_BLep_DRP4 == best_WJa_DRP4 && !(best_BLep_DRP4 == 0) && !(best_BHad_DRP4 == best_BLep_DRP4) && !(best_BLep_DRP4 == best_WJb_DRP4) ){//only reco BLep and WJa merged
				Merged_BLepWJa_massDivpt_DRP4_hist->Fill(best_BLep_DRP4->M()/best_BLep_DRP4->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJa_perm_DRP4_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJa_perm_mass_DRP4_hist->Fill(best_BLep_DRP4->M());//merged b and wjet mass
				if( !(best_WJb_DRP4 == 0) && !(best_WJb_DRP4 == best_BLep_DRP4) && !(best_WJb_DRP4 == best_BHad_DRP4) ){//reco WJb exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M(), (*best_BLep_DRP4+*best_WJb_DRP4).M());
					Merged_BLepWJa_perm_and_WJb_mass_DRP4_hist->Fill((*best_BLep_DRP4+*best_WJb_DRP4).M());//comb invariant mass of merged jet and other wjet
					//Merged_WJb_mass_DRP4_hist->Fill(best_WJb_DRP4->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}
			if( best_BLep_DRP4 == best_WJb_DRP4 && !(best_BLep_DRP4 == 0) && !(best_BHad_DRP4 == best_BLep_DRP4) && !(best_BLep_DRP4 == best_WJa_DRP4) ){//only reco BLep and WJb merged
				Merged_BLepWJb_massDivpt_DRP4_hist->Fill(best_BLep_DRP4->M()/best_BLep_DRP4->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJb_perm_DRP4_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJb_perm_mass_DRP4_hist->Fill(best_BLep_DRP4->M());//merged b and wjet mass
				if( !(best_WJa_DRP4 == 0) && !(best_WJa_DRP4 == best_BLep_DRP4) && !(best_WJa_DRP4 == best_BHad_DRP4) ){//reco WJa exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP4_hist->Fill(best_BLep_DRP4->M(), (*best_BLep_DRP4+*best_WJa_DRP4).M());
					Merged_BLepWJb_perm_and_WJa_mass_DRP4_hist->Fill((*best_BLep_DRP4+*best_WJa_DRP4).M());//comb mass of merged jet and other wjet
				//	Merged_WJa_mass_DRP4_hist->Fill(best_WJa_DRP4->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP4_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP4_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP4_hist->Fill(thad->Pt());
			}


				//DR=0.5
			list<IDJet*> bhad_DRP5_list;
			float bhad_DRP5 = 1e10;// itialize dr to high number
			if( object_selector_.clean_jets().size() < 3 ) continue;
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BHad == 0) ){
					if( (*jets)->DeltaR(*BHad) > 0.5 ) continue;
					if( (*jets)->DeltaR(*BHad) < bhad_DRP5 ){
						bhad_DRP5 = (*jets)->DeltaR(*BHad);
						bhad_DRP5_list.push_back(*jets);
						best_BHad_DRP5 = (bhad_DRP5_list.back());
					}
				}
			}
			list<IDJet*> blep_DRP5_list;
			float blep_DRP5 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BLep == 0) ){
					if( (*jets)->DeltaR(*BLep) > 0.5 ) continue;
					if( (*jets)->DeltaR(*BLep) < blep_DRP5 ){
						blep_DRP5 = (*jets)->DeltaR(*BLep);
						blep_DRP5_list.push_back(*jets);
						best_BLep_DRP5 = (blep_DRP5_list.back());
					}
				}
			}
			list<IDJet*> wja_DRP5_list;
			float wja_DRP5 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJa == 0) ){
					if( (*jets)->DeltaR(*WJa) > 0.5 ) continue;
					if( (*jets)->DeltaR(*WJa) < wja_DRP5 ){
						wja_DRP5 = (*jets)->DeltaR(*WJa);
						wja_DRP5_list.push_back(*jets);
						best_WJa_DRP5 = (wja_DRP5_list.back());
					}
				}
			}
			list<IDJet*> wjb_DRP5_list;
			float wjb_DRP5 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJb == 0) ){
					if( (*jets)->DeltaR(*WJb) > 0.5 ) continue;
					if( (*jets)->DeltaR(*WJb) < wjb_DRP5 ){
						wjb_DRP5 = (*jets)->DeltaR(*WJb);
						wjb_DRP5_list.push_back(*jets);
						best_WJb_DRP5 = (wjb_DRP5_list.back());
					}
				}
			}

			//Reco Obj. Kinematic Vars
			if( !(best_BHad_DRP5 == 0) ){//reco BHad exists
				Matched_perm_BHad_pt_DRP5_hist->Fill(best_BHad_DRP5->Pt());
				Matched_perm_BHad_eta_DRP5_hist->Fill(best_BHad_DRP5->Eta());
//						++nmatched_objects;
				if( best_BHad_DRP5 == best_BLep_DRP5 || best_BHad_DRP5 == best_WJa_DRP5 || best_BHad_DRP5 == best_WJb_DRP5 ){//reco BHad merged with something
					if( best_BHad_DRP5->BTagId(cut_loose_b_) ) BTag_lp5_BHad_loose_pass_hist->Fill(ttbar.M());
					if( best_BHad_DRP5->BTagId(cut_medium_b_) ) BTag_lp5_BHad_medium_pass_hist->Fill(ttbar.M());	
					if( best_BHad_DRP5->BTagId(cut_tight_b_) ) BTag_lp5_BHad_tight_pass_hist->Fill(ttbar.M());
					if( !best_BHad_DRP5->BTagId(cut_loose_b_) ) BTag_lp5_BHad_loose_fail_hist->Fill(ttbar.M());
					if( !best_BHad_DRP5->BTagId(cut_medium_b_) ) BTag_lp5_BHad_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BHad_DRP5->BTagId(cut_tight_b_) ) BTag_lp5_BHad_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_BLep_DRP5 == 0) ){//reco BLep exists
				Matched_perm_BLep_pt_DRP5_hist->Fill(best_BLep_DRP5->Pt());
				Matched_perm_BLep_eta_DRP5_hist->Fill(best_BLep_DRP5->Eta());
//				++nmatched_objects; 
				if( best_BLep_DRP5 == best_BHad_DRP5 || best_BLep_DRP5 == best_WJa_DRP5 || best_BLep_DRP5 == best_WJb_DRP5 ){//reco BLep merged with something
					if( best_BLep_DRP5->BTagId(cut_loose_b_) ) BTag_lp5_BLep_loose_pass_hist->Fill(ttbar.M());
					if( best_BLep_DRP5->BTagId(cut_medium_b_) ) BTag_lp5_BLep_medium_pass_hist->Fill(ttbar.M());	
					if( best_BLep_DRP5->BTagId(cut_tight_b_) ) BTag_lp5_BLep_tight_pass_hist->Fill(ttbar.M());
					if( !best_BLep_DRP5->BTagId(cut_loose_b_) ) BTag_lp5_BLep_loose_fail_hist->Fill(ttbar.M());
					if( !best_BLep_DRP5->BTagId(cut_medium_b_) ) BTag_lp5_BLep_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BLep_DRP5->BTagId(cut_tight_b_) ) BTag_lp5_BLep_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_WJa_DRP5 == 0) ){//reco WJa exists
//						++nmatched_objects; 
				Matched_perm_WJa_pt_DRP5_hist->Fill(best_WJa_DRP5->Pt());
				Matched_perm_WJa_eta_DRP5_hist->Fill(best_WJa_DRP5->Eta());
                        	if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                        	        Matched_perm_WJa_ttbarM700_pt_DRP5_hist->Fill(best_WJa_DRP5->Pt());
                        	        Matched_perm_WJa_ttbarM700_frac_p_DRP5_hist->Fill(best_WJa_DRP5->P()/(ttbar.M()/2));
                        	}
                        	if( ttbar.M() >= 1000 ){
                        	        Matched_perm_WJa_ttbarM1000_pt_DRP5_hist->Fill(best_WJa_DRP5->Pt());
                        	        Matched_perm_WJa_ttbarM1000_frac_p_DRP5_hist->Fill(best_WJa_DRP5->P()/(ttbar.M()/2));
                        	}
			}
			if( !(best_WJb_DRP5 == 0) ){//reco WJb exists
//						++nmatched_objects; 
				Matched_perm_WJb_pt_DRP5_hist->Fill(best_WJb_DRP5->Pt());
				Matched_perm_WJb_eta_DRP5_hist->Fill(best_WJb_DRP5->Eta());
				if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
				        Matched_perm_WJb_ttbarM700_pt_DRP5_hist->Fill(best_WJb_DRP5->Pt());
				        Matched_perm_WJb_ttbarM700_frac_p_DRP5_hist->Fill(best_WJb_DRP5->P()/(ttbar.M()/2));
				}
				if( ttbar.M() >= 1000 ){
				        Matched_perm_WJb_ttbarM1000_pt_DRP5_hist->Fill(best_WJb_DRP5->Pt());
				        Matched_perm_WJb_ttbarM1000_frac_p_DRP5_hist->Fill(best_WJb_DRP5->P()/(ttbar.M()/2));
				}
			}
			if( !(best_BLep_DRP5 == 0) && !(best_BLep_DRP5 == best_BHad_DRP5) && !(best_BLep_DRP5 == best_WJa_DRP5) && !(best_BLep_DRP5 == best_WJb_DRP5) ){//reco BLep not merged with anything
				if( best_BLep_DRP5->BTagId(cut_loose_b_) ) BTag_gp5_BLep_loose_pass_hist->Fill(ttbar.M());
				if( best_BLep_DRP5->BTagId(cut_medium_b_) ) BTag_gp5_BLep_medium_pass_hist->Fill(ttbar.M());	
				if( best_BLep_DRP5->BTagId(cut_tight_b_) ) BTag_gp5_BLep_tight_pass_hist->Fill(ttbar.M());
				if( !best_BLep_DRP5->BTagId(cut_loose_b_) ) BTag_gp5_BLep_loose_fail_hist->Fill(ttbar.M());
				if( !best_BLep_DRP5->BTagId(cut_medium_b_) ) BTag_gp5_BLep_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BLep_DRP5->BTagId(cut_tight_b_) ) BTag_gp5_BLep_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BLep_LepBLep_DRP5_hist->Fill(DR_LepBLep);
	
				Unmerged_BLep_massDivpt_DRP5_hist->Fill(best_BLep_DRP5->M()/best_BLep_DRP5->Pt());
				if( !(best_WJa_DRP5 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M(),(*best_BLep_DRP5+*best_WJa_DRP5).M());
				if( !(best_WJb_DRP5 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M(),(*best_BLep_DRP5+*best_WJb_DRP5).M());
				if( !(best_WJa_DRP5 == 0) && !(best_WJb_DRP5 == 0) ){
					if( best_WJa_DRP5->M() > best_WJb_DRP5->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M(),(*best_BLep_DRP5+*best_WJa_DRP5).M());
					if( best_WJa_DRP5->M() < best_WJb_DRP5->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M(),(*best_BLep_DRP5+*best_WJb_DRP5).M());
				}
			}
			if( !(best_BHad_DRP5 == 0) && !(best_BHad_DRP5 == best_WJa_DRP5) && !(best_BHad_DRP5 == best_WJb_DRP5) && !(best_BHad_DRP5 == best_BLep_DRP5) ){//reco BHad not merged with anything
				if( best_BHad_DRP5->BTagId(cut_loose_b_) ) BTag_gp5_BHad_loose_pass_hist->Fill(ttbar.M());
				if( best_BHad_DRP5->BTagId(cut_medium_b_) ) BTag_gp5_BHad_medium_pass_hist->Fill(ttbar.M());	
				if( best_BHad_DRP5->BTagId(cut_tight_b_) ) BTag_gp5_BHad_tight_pass_hist->Fill(ttbar.M());
				if( !best_BHad_DRP5->BTagId(cut_loose_b_) ) BTag_gp5_BHad_loose_fail_hist->Fill(ttbar.M());
				if( !best_BHad_DRP5->BTagId(cut_medium_b_) ) BTag_gp5_BHad_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BHad_DRP5->BTagId(cut_tight_b_) ) BTag_gp5_BHad_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BHad_LepBLep_DRP5_hist->Fill(DR_LepBLep);

				Unmerged_BHad_massDivpt_DRP5_hist->Fill(best_BHad_DRP5->M()/best_BHad_DRP5->Pt());
				if( !(best_WJa_DRP5 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M(),(*best_BHad_DRP5+*best_WJa_DRP5).M());
				if( !(best_WJb_DRP5 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M(),(*best_BHad_DRP5+*best_WJb_DRP5).M());
				if( !(best_WJa_DRP5 == 0) && !(best_WJb_DRP5 == 0) ){
					if( best_WJa_DRP5->M() > best_WJb_DRP5->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M(),(*best_BHad_DRP5+*best_WJa_DRP5).M());
					if( best_WJa_DRP5->M() < best_WJb_DRP5->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M(),(*best_BHad_DRP5+*best_WJb_DRP5).M());
				}
			}
			if( !(best_WJa_DRP5 == 0) && !(best_WJa_DRP5 == best_BHad_DRP5) && !(best_WJa_DRP5 == best_WJb_DRP5) && !(best_WJa_DRP5 == best_BLep_DRP5) ){//reco WJa not merged with anything
				Unmerged_WJa_massDivpt_DRP5_hist->Fill(best_WJa_DRP5->M()/best_WJa_DRP5->Pt());
				Unmerged_WJa_mass_DRP5_hist->Fill(best_WJa_DRP5->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJa_LepBLep_DRP5_hist->Fill(DR_LepBLep);
			}
			if( !(best_WJb_DRP5 == 0) && !(best_WJb_DRP5 == best_BHad_DRP5) && !(best_WJb_DRP5 == best_WJa_DRP5) && !(best_WJb_DRP5 == best_BLep_DRP5) ){//reco WJb not merged with anything
				Unmerged_WJb_massDivpt_DRP5_hist->Fill(best_WJb_DRP5->M()/best_WJb_DRP5->Pt());
				Unmerged_WJb_mass_DRP5_hist->Fill(best_WJb_DRP5->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJb_LepBLep_DRP5_hist->Fill(DR_LepBLep);
			}
			if( !(best_BHad_DRP5 == 0) && !(best_WJa_DRP5 == 0) ){ //reco BHad and WJa exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJa_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJa_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJa_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}
			if( !(best_BHad_DRP5 == 0) && !(best_WJb_DRP5 == 0) ){ //reco BHad and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJb_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJb_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJb_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}
			if( !(best_WJa_DRP5 == 0) && !(best_WJb_DRP5 == 0) ){ //reco WJa and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_WJaWJb_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_WJaWJb_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_WJaWJb_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP5 == best_WJa_DRP5 && !(best_BHad_DRP5 == 0) && !(best_BHad_DRP5 == best_BLep_DRP5) && !(best_BHad_DRP5 == best_WJb_DRP5) ){//only reco BHad and WJa merged
				Merged_BHadWJa_massDivpt_DRP5_hist->Fill(best_BHad_DRP5->M()/best_BHad_DRP5->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJa_perm_DRP5_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJa_perm_mass_DRP5_hist->Fill(best_BHad_DRP5->M());//merged b and wjet mass
				if( !(best_WJb_DRP5 == 0) && !(best_WJb_DRP5 == best_BHad_DRP5) && !(best_WJb_DRP5 == best_BLep_DRP5) ){//reco WJb exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M(), (*best_BHad_DRP5+*best_WJb_DRP5).M());
					Merged_BHadWJa_perm_and_WJb_mass_DRP5_hist->Fill((*best_BHad_DRP5+*best_WJb_DRP5).M());//comb invariant mass of merged jet and other wjet
					Merged_WJb_mass_DRP5_hist->Fill(best_WJb_DRP5->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP5 == best_WJb_DRP5 && !(best_BHad_DRP5 == 0) && !(best_BHad_DRP5 == best_BLep_DRP5) && !(best_BHad_DRP5 == best_WJa_DRP5) ){//only reco BHad and WJb merged
				Merged_BHadWJb_massDivpt_DRP5_hist->Fill(best_BHad_DRP5->M()/best_BHad_DRP5->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJb_perm_DRP5_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJb_perm_mass_DRP5_hist->Fill(best_BHad_DRP5->M());//merged b and wjet mass
				if( !(best_WJa_DRP5 == 0) && !(best_WJa_DRP5 == best_BHad_DRP5) && !(best_WJa_DRP5 == best_BLep_DRP5) ){//reco WJa exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP5_hist->Fill(best_BHad_DRP5->M(), (*best_BHad_DRP5+*best_WJa_DRP5).M());
					Merged_BHadWJb_perm_and_WJa_mass_DRP5_hist->Fill((*best_BHad_DRP5+*best_WJa_DRP5).M());//comb mass of merged jet and other wjet
					Merged_WJa_mass_DRP5_hist->Fill(best_WJa_DRP5->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}
			if( best_WJa_DRP5 == best_WJb_DRP5 && !(best_WJa_DRP5 == 0) && !(best_WJa_DRP5 == best_BHad_DRP5) && !(best_WJa_DRP5 == best_BLep_DRP5) ){//only reco WJa and WJb merged
						if( object_selector_.clean_jets().size() == 3 ) Matched_WJaWJb_ptthad_3J_DRP5_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() == 4 ) Matched_WJaWJb_ptthad_4J_DRP5_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() > 4 ) Matched_WJaWJb_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}

			if( best_BLep_DRP5 == best_WJa_DRP5 && !(best_BLep_DRP5 == 0) && !(best_BHad_DRP5 == best_BLep_DRP5) && !(best_BLep_DRP5 == best_WJb_DRP5) ){//only reco BLep and WJa merged
				Merged_BLepWJa_massDivpt_DRP5_hist->Fill(best_BLep_DRP5->M()/best_BLep_DRP5->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJa_perm_DRP5_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJa_perm_mass_DRP5_hist->Fill(best_BLep_DRP5->M());//merged b and wjet mass
				if( !(best_WJb_DRP5 == 0) && !(best_WJb_DRP5 == best_BLep_DRP5) && !(best_WJb_DRP5 == best_BHad_DRP5) ){//reco WJb exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M(), (*best_BLep_DRP5+*best_WJb_DRP5).M());
					Merged_BLepWJa_perm_and_WJb_mass_DRP5_hist->Fill((*best_BLep_DRP5+*best_WJb_DRP5).M());//comb invariant mass of merged jet and other wjet
					//Merged_WJb_mass_DRP5_hist->Fill(best_WJb_DRP5->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}
			if( best_BLep_DRP5 == best_WJb_DRP5 && !(best_BLep_DRP5 == 0) && !(best_BHad_DRP5 == best_BLep_DRP5) && !(best_BLep_DRP5 == best_WJa_DRP5) ){//only reco BLep and WJb merged
				Merged_BLepWJb_massDivpt_DRP5_hist->Fill(best_BLep_DRP5->M()/best_BLep_DRP5->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJb_perm_DRP5_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJb_perm_mass_DRP5_hist->Fill(best_BLep_DRP5->M());//merged b and wjet mass
				if( !(best_WJa_DRP5 == 0) && !(best_WJa_DRP5 == best_BLep_DRP5) && !(best_WJa_DRP5 == best_BHad_DRP5) ){//reco WJa exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP5_hist->Fill(best_BLep_DRP5->M(), (*best_BLep_DRP5+*best_WJa_DRP5).M());
					Merged_BLepWJb_perm_and_WJa_mass_DRP5_hist->Fill((*best_BLep_DRP5+*best_WJa_DRP5).M());//comb mass of merged jet and other wjet
				//	Merged_WJa_mass_DRP4_hist->Fill(best_WJa_DRP5->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP5_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP5_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP5_hist->Fill(thad->Pt());
			}


				//DR=0.6
			list<IDJet*> bhad_DRP6_list;
			float bhad_DRP6 = 1e10;// itialize dr to high number
			if( object_selector_.clean_jets().size() < 3 ) continue;
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BHad == 0) ){
					if( (*jets)->DeltaR(*BHad) > 0.6 ) continue;
					if( (*jets)->DeltaR(*BHad) < bhad_DRP6 ){
						bhad_DRP6 = (*jets)->DeltaR(*BHad);
						bhad_DRP6_list.push_back(*jets);
						best_BHad_DRP6 = (bhad_DRP6_list.back());
					}
				}
			}
			list<IDJet*> blep_DRP6_list;
			float blep_DRP6 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BLep == 0) ){
					if( (*jets)->DeltaR(*BLep) > 0.6 ) continue;
					if( (*jets)->DeltaR(*BLep) < blep_DRP6 ){
						blep_DRP6 = (*jets)->DeltaR(*BLep);
						blep_DRP6_list.push_back(*jets);
						best_BLep_DRP6 = (blep_DRP6_list.back());
					}
				}
			}
			list<IDJet*> wja_DRP6_list;
			float wja_DRP6 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJa == 0) ){
					if( (*jets)->DeltaR(*WJa) > 0.6 ) continue;
					if( (*jets)->DeltaR(*WJa) < wja_DRP6 ){
						wja_DRP6 = (*jets)->DeltaR(*WJa);
						wja_DRP6_list.push_back(*jets);
						best_WJa_DRP6 = (wja_DRP6_list.back());
					}
				}
			}
			list<IDJet*> wjb_DRP6_list;
			float wjb_DRP6 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJb == 0) ){
					if( (*jets)->DeltaR(*WJb) > 0.6 ) continue;
					if( (*jets)->DeltaR(*WJb) < wjb_DRP6 ){
						wjb_DRP6 = (*jets)->DeltaR(*WJb);
						wjb_DRP6_list.push_back(*jets);
						best_WJb_DRP6 = (wjb_DRP6_list.back());
					}
				}
			}

			//Reco Obj. Kinematic Vars
			if( !(best_BHad_DRP6 == 0) ){//reco BHad exists
				Matched_perm_BHad_pt_DRP6_hist->Fill(best_BHad_DRP6->Pt());
				Matched_perm_BHad_eta_DRP6_hist->Fill(best_BHad_DRP6->Eta());
//						++nmatched_objects;
				if( best_BHad_DRP6 == best_BLep_DRP6 || best_BHad_DRP6 == best_WJa_DRP6 || best_BHad_DRP6 == best_WJb_DRP6 ){//reco BHad merged with something
					if( best_BHad_DRP6->BTagId(cut_loose_b_) ) BTag_lp6_BHad_loose_pass_hist->Fill(ttbar.M());
					if( best_BHad_DRP6->BTagId(cut_medium_b_) ) BTag_lp6_BHad_medium_pass_hist->Fill(ttbar.M());	
					if( best_BHad_DRP6->BTagId(cut_tight_b_) ) BTag_lp6_BHad_tight_pass_hist->Fill(ttbar.M());
					if( !best_BHad_DRP6->BTagId(cut_loose_b_) ) BTag_lp6_BHad_loose_fail_hist->Fill(ttbar.M());
					if( !best_BHad_DRP6->BTagId(cut_medium_b_) ) BTag_lp6_BHad_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BHad_DRP6->BTagId(cut_tight_b_) ) BTag_lp6_BHad_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_BLep_DRP6 == 0) ){//reco BLep exists
				Matched_perm_BLep_pt_DRP6_hist->Fill(best_BLep_DRP6->Pt());
				Matched_perm_BLep_eta_DRP6_hist->Fill(best_BLep_DRP6->Eta());
//				++nmatched_objects; 
				if( best_BLep_DRP6 == best_BHad_DRP6 || best_BLep_DRP6 == best_WJa_DRP6 || best_BLep_DRP6 == best_WJb_DRP6 ){//reco BLep merged with something
					if( best_BLep_DRP6->BTagId(cut_loose_b_) ) BTag_lp6_BLep_loose_pass_hist->Fill(ttbar.M());
					if( best_BLep_DRP6->BTagId(cut_medium_b_) ) BTag_lp6_BLep_medium_pass_hist->Fill(ttbar.M());	
					if( best_BLep_DRP6->BTagId(cut_tight_b_) ) BTag_lp6_BLep_tight_pass_hist->Fill(ttbar.M());
					if( !best_BLep_DRP6->BTagId(cut_loose_b_) ) BTag_lp6_BLep_loose_fail_hist->Fill(ttbar.M());
					if( !best_BLep_DRP6->BTagId(cut_medium_b_) ) BTag_lp6_BLep_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BLep_DRP6->BTagId(cut_tight_b_) ) BTag_lp6_BLep_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_WJa_DRP6 == 0) ){//reco WJa exists
//						++nmatched_objects; 
				Matched_perm_WJa_pt_DRP6_hist->Fill(best_WJa_DRP6->Pt());
				Matched_perm_WJa_eta_DRP6_hist->Fill(best_WJa_DRP6->Eta());
                        	if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                        	        Matched_perm_WJa_ttbarM700_pt_DRP6_hist->Fill(best_WJa_DRP6->Pt());
                        	        Matched_perm_WJa_ttbarM700_frac_p_DRP6_hist->Fill(best_WJa_DRP6->P()/(ttbar.M()/2));
                        	}
                        	if( ttbar.M() >= 1000 ){
                        	        Matched_perm_WJa_ttbarM1000_pt_DRP6_hist->Fill(best_WJa_DRP6->Pt());
                        	        Matched_perm_WJa_ttbarM1000_frac_p_DRP6_hist->Fill(best_WJa_DRP6->P()/(ttbar.M()/2));
                        	}
			}
			if( !(best_WJb_DRP6 == 0) ){//reco WJb exists
//						++nmatched_objects; 
				Matched_perm_WJb_pt_DRP6_hist->Fill(best_WJb_DRP6->Pt());
				Matched_perm_WJb_eta_DRP6_hist->Fill(best_WJb_DRP6->Eta());
				if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
				        Matched_perm_WJb_ttbarM700_pt_DRP6_hist->Fill(best_WJb_DRP6->Pt());
				        Matched_perm_WJb_ttbarM700_frac_p_DRP6_hist->Fill(best_WJb_DRP6->P()/(ttbar.M()/2));
				}
				if( ttbar.M() >= 1000 ){
				        Matched_perm_WJb_ttbarM1000_pt_DRP6_hist->Fill(best_WJb_DRP6->Pt());
				        Matched_perm_WJb_ttbarM1000_frac_p_DRP6_hist->Fill(best_WJb_DRP6->P()/(ttbar.M()/2));
				}
			}
			if( !(best_BLep_DRP6 == 0) && !(best_BLep_DRP6 == best_BHad_DRP6) && !(best_BLep_DRP6 == best_WJa_DRP6) && !(best_BLep_DRP6 == best_WJb_DRP6) ){//reco BLep not merged with anything
				if( best_BLep_DRP6->BTagId(cut_loose_b_) ) BTag_gp6_BLep_loose_pass_hist->Fill(ttbar.M());
				if( best_BLep_DRP6->BTagId(cut_medium_b_) ) BTag_gp6_BLep_medium_pass_hist->Fill(ttbar.M());	
				if( best_BLep_DRP6->BTagId(cut_tight_b_) ) BTag_gp6_BLep_tight_pass_hist->Fill(ttbar.M());
				if( !best_BLep_DRP6->BTagId(cut_loose_b_) ) BTag_gp6_BLep_loose_fail_hist->Fill(ttbar.M());
				if( !best_BLep_DRP6->BTagId(cut_medium_b_) ) BTag_gp6_BLep_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BLep_DRP6->BTagId(cut_tight_b_) ) BTag_gp6_BLep_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BLep_LepBLep_DRP6_hist->Fill(DR_LepBLep);
	
				Unmerged_BLep_massDivpt_DRP6_hist->Fill(best_BLep_DRP6->M()/best_BLep_DRP6->Pt());
				if( !(best_WJa_DRP6 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M(),(*best_BLep_DRP6+*best_WJa_DRP6).M());
				if( !(best_WJb_DRP6 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M(),(*best_BLep_DRP6+*best_WJb_DRP6).M());
				if( !(best_WJa_DRP6 == 0) && !(best_WJb_DRP6 == 0) ){
					if( best_WJa_DRP6->M() > best_WJb_DRP6->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M(),(*best_BLep_DRP6+*best_WJa_DRP6).M());
					if( best_WJa_DRP6->M() < best_WJb_DRP6->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M(),(*best_BLep_DRP6+*best_WJb_DRP6).M());
				}

			}
			if( !(best_BHad_DRP6 == 0) && !(best_BHad_DRP6 == best_WJa_DRP6) && !(best_BHad_DRP6 == best_WJb_DRP6) && !(best_BHad_DRP6 == best_BLep_DRP6) ){//reco BHad not merged with anything
				if( best_BHad_DRP6->BTagId(cut_loose_b_) ) BTag_gp6_BHad_loose_pass_hist->Fill(ttbar.M());
				if( best_BHad_DRP6->BTagId(cut_medium_b_) ) BTag_gp6_BHad_medium_pass_hist->Fill(ttbar.M());	
				if( best_BHad_DRP6->BTagId(cut_tight_b_) ) BTag_gp6_BHad_tight_pass_hist->Fill(ttbar.M());
				if( !best_BHad_DRP6->BTagId(cut_loose_b_) ) BTag_gp6_BHad_loose_fail_hist->Fill(ttbar.M());
				if( !best_BHad_DRP6->BTagId(cut_medium_b_) ) BTag_gp6_BHad_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BHad_DRP6->BTagId(cut_tight_b_) ) BTag_gp6_BHad_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BHad_LepBLep_DRP6_hist->Fill(DR_LepBLep);

				Unmerged_BHad_massDivpt_DRP6_hist->Fill(best_BHad_DRP6->M()/best_BHad_DRP6->Pt());
				if( !(best_WJa_DRP6 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M(),(*best_BHad_DRP6+*best_WJa_DRP6).M());
				if( !(best_WJb_DRP6 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M(),(*best_BHad_DRP6+*best_WJb_DRP6).M());
				if( !(best_WJa_DRP6 == 0) && !(best_WJb_DRP6 == 0) ){
					if( best_WJa_DRP6->M() > best_WJb_DRP6->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M(),(*best_BHad_DRP6+*best_WJa_DRP6).M());
					if( best_WJa_DRP6->M() < best_WJb_DRP6->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M(),(*best_BHad_DRP6+*best_WJb_DRP6).M());
				}
			}
			if( !(best_WJa_DRP6 == 0) && !(best_WJa_DRP6 == best_BHad_DRP6) && !(best_WJa_DRP6 == best_WJb_DRP6) && !(best_WJa_DRP6 == best_BLep_DRP6) ){//reco WJa not merged with anything
				Unmerged_WJa_massDivpt_DRP6_hist->Fill(best_WJa_DRP6->M()/best_WJa_DRP6->Pt());
				Unmerged_WJa_mass_DRP6_hist->Fill(best_WJa_DRP6->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJa_LepBLep_DRP6_hist->Fill(DR_LepBLep);
			}
			if( !(best_WJb_DRP6 == 0) && !(best_WJb_DRP6 == best_BHad_DRP6) && !(best_WJb_DRP6 == best_WJa_DRP6) && !(best_WJb_DRP6 == best_BLep_DRP6) ){//reco WJb not merged with anything
				Unmerged_WJb_massDivpt_DRP6_hist->Fill(best_WJb_DRP6->M()/best_WJb_DRP6->Pt());
				Unmerged_WJb_mass_DRP6_hist->Fill(best_WJb_DRP6->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJb_LepBLep_DRP6_hist->Fill(DR_LepBLep);
			}
			if( !(best_BHad_DRP6 == 0) && !(best_WJa_DRP6 == 0) ){ //reco BHad and WJa exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJa_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJa_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJa_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());

			}
			if( !(best_BHad_DRP6 == 0) && !(best_WJb_DRP6 == 0) ){ //reco BHad and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJb_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJb_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJb_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}
			if( !(best_WJa_DRP6 == 0) && !(best_WJb_DRP6 == 0) ){ //reco WJa and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_WJaWJb_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_WJaWJb_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_WJaWJb_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP6 == best_WJa_DRP6 && !(best_BHad_DRP6 == 0) && !(best_BHad_DRP6 == best_BLep_DRP6) && !(best_BHad_DRP6 == best_WJb_DRP6) ){//only reco BHad and WJa merged
				Merged_BHadWJa_massDivpt_DRP6_hist->Fill(best_BHad_DRP6->M()/best_BHad_DRP6->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJa_perm_DRP6_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJa_perm_mass_DRP6_hist->Fill(best_BHad_DRP6->M());//merged b and wjet mass
				if( !(best_WJb_DRP6 == 0) && !(best_WJb_DRP6 == best_BHad_DRP6) && !(best_WJb_DRP6 == best_BLep_DRP6) ){//reco WJb exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M(), (*best_BHad_DRP6+*best_WJb_DRP6).M());
					Merged_BHadWJa_perm_and_WJb_mass_DRP6_hist->Fill((*best_BHad_DRP6+*best_WJb_DRP6).M());//comb invariant mass of merged jet and other wjet
					Merged_WJb_mass_DRP6_hist->Fill(best_WJb_DRP6->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP6 == best_WJb_DRP6 && !(best_BHad_DRP6 == 0) && !(best_BHad_DRP6 == best_BLep_DRP6) && !(best_BHad_DRP6 == best_WJa_DRP6) ){//only reco BHad and WJb merged
				Merged_BHadWJb_massDivpt_DRP6_hist->Fill(best_BHad_DRP6->M()/best_BHad_DRP6->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJb_perm_DRP6_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJb_perm_mass_DRP6_hist->Fill(best_BHad_DRP6->M());//merged b and wjet mass
				if( !(best_WJa_DRP6 == 0) && !(best_WJa_DRP6 == best_BHad_DRP6) && !(best_WJa_DRP6 == best_BLep_DRP6) ){//reco WJa exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP6_hist->Fill(best_BHad_DRP6->M(), (*best_BHad_DRP6+*best_WJa_DRP6).M());
					Merged_BHadWJb_perm_and_WJa_mass_DRP6_hist->Fill((*best_BHad_DRP6+*best_WJa_DRP6).M());//comb mass of merged jet and other wjet
					Merged_WJa_mass_DRP6_hist->Fill(best_WJa_DRP6->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}
			if( best_WJa_DRP6 == best_WJb_DRP6 && !(best_WJa_DRP6 == 0) && !(best_WJa_DRP6 == best_BHad_DRP6) && !(best_WJa_DRP6 == best_BLep_DRP6) ){//only reco WJa and WJb merged
						if( object_selector_.clean_jets().size() == 3 ) Matched_WJaWJb_ptthad_3J_DRP6_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() == 4 ) Matched_WJaWJb_ptthad_4J_DRP6_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() > 4 ) Matched_WJaWJb_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}

			if( best_BLep_DRP6 == best_WJa_DRP6 && !(best_BLep_DRP6 == 0) && !(best_BHad_DRP6 == best_BLep_DRP6) && !(best_BLep_DRP6 == best_WJb_DRP6) ){//only reco BLep and WJa merged
				Merged_BLepWJa_massDivpt_DRP6_hist->Fill(best_BLep_DRP6->M()/best_BLep_DRP6->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJa_perm_DRP6_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJa_perm_mass_DRP6_hist->Fill(best_BLep_DRP6->M());//merged b and wjet mass
				if( !(best_WJb_DRP6 == 0) && !(best_WJb_DRP6 == best_BLep_DRP6) && !(best_WJb_DRP6 == best_BHad_DRP6) ){//reco WJb exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M(), (*best_BLep_DRP6+*best_WJb_DRP6).M());
					Merged_BLepWJa_perm_and_WJb_mass_DRP6_hist->Fill((*best_BLep_DRP6+*best_WJb_DRP6).M());//comb invariant mass of merged jet and other wjet
					//Merged_WJb_mass_DRP6_hist->Fill(best_WJb_DRP6->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}
			if( best_BLep_DRP6 == best_WJb_DRP6 && !(best_BLep_DRP6 == 0) && !(best_BHad_DRP6 == best_BLep_DRP6) && !(best_BLep_DRP6 == best_WJa_DRP6) ){//only reco BLep and WJb merged
				Merged_BLepWJb_massDivpt_DRP6_hist->Fill(best_BLep_DRP6->M()/best_BLep_DRP6->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJb_perm_DRP6_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJb_perm_mass_DRP6_hist->Fill(best_BLep_DRP6->M());//merged b and wjet mass
				if( !(best_WJa_DRP6 == 0) && !(best_WJa_DRP6 == best_BLep_DRP6) && !(best_WJa_DRP6 == best_BHad_DRP6) ){//reco WJa exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP6_hist->Fill(best_BLep_DRP6->M(), (*best_BLep_DRP6+*best_WJa_DRP6).M());
					Merged_BLepWJb_perm_and_WJa_mass_DRP6_hist->Fill((*best_BLep_DRP6+*best_WJa_DRP6).M());//comb mass of merged jet and other wjet
				//	Merged_WJa_mass_DRP6_hist->Fill(best_WJa_DRP6->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP6_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP6_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP6_hist->Fill(thad->Pt());
			}


				//DR=0.8
			list<IDJet*> bhad_DRP8_list;
			float bhad_DRP8 = 1e10;// itialize dr to high number
			if( object_selector_.clean_jets().size() < 3 ) continue;
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BHad == 0) ){
					if( (*jets)->DeltaR(*BHad) > 0.8 ) continue;
					if( (*jets)->DeltaR(*BHad) < bhad_DRP8 ){
						bhad_DRP8 = (*jets)->DeltaR(*BHad);
						bhad_DRP8_list.push_back(*jets);
						best_BHad_DRP8 = (bhad_DRP8_list.back());
					}
				}
			}
			list<IDJet*> blep_DRP8_list;
			float blep_DRP8 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(BLep == 0) ){
					if( (*jets)->DeltaR(*BLep) > 0.8 ) continue;
					if( (*jets)->DeltaR(*BLep) < blep_DRP8 ){
						blep_DRP8 = (*jets)->DeltaR(*BLep);
						blep_DRP8_list.push_back(*jets);
						best_BLep_DRP8 = (blep_DRP8_list.back());
					}
				}
			}
			list<IDJet*> wja_DRP8_list;
			float wja_DRP8 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJa == 0) ){
					if( (*jets)->DeltaR(*WJa) > 0.8 ) continue;
					if( (*jets)->DeltaR(*WJa) < wja_DRP8 ){
						wja_DRP8 = (*jets)->DeltaR(*WJa);
						wja_DRP8_list.push_back(*jets);
						best_WJa_DRP8 = (wja_DRP8_list.back());
					}
				}
			}
			list<IDJet*> wjb_DRP8_list;
			float wjb_DRP8 = 1e10;// itialize dr to high number
			for(vector<IDJet*>::const_iterator jets = object_selector_.clean_jets().begin(); jets != object_selector_.clean_jets().end(); ++jets){
				if( !(WJb == 0) ){
					if( (*jets)->DeltaR(*WJb) > 0.8 ) continue;
					if( (*jets)->DeltaR(*WJb) < wjb_DRP8 ){
						wjb_DRP8 = (*jets)->DeltaR(*WJb);
						wjb_DRP8_list.push_back(*jets);
						best_WJb_DRP8 = (wjb_DRP8_list.back());
					}
				}
			}

			//Reco Obj. Kinematic Vars
			if( !(best_BHad_DRP8 == 0) ){//reco BHad exists
				Matched_perm_BHad_pt_DRP8_hist->Fill(best_BHad_DRP8->Pt());
				Matched_perm_BHad_eta_DRP8_hist->Fill(best_BHad_DRP8->Eta());
//						++nmatched_objects;
				if( best_BHad_DRP8 == best_BLep_DRP8 || best_BHad_DRP8 == best_WJa_DRP8 || best_BHad_DRP8 == best_WJb_DRP8 ){//reco BHad merged with something
					if( best_BHad_DRP8->BTagId(cut_loose_b_) ) BTag_lp8_BHad_loose_pass_hist->Fill(ttbar.M());
					if( best_BHad_DRP8->BTagId(cut_medium_b_) ) BTag_lp8_BHad_medium_pass_hist->Fill(ttbar.M());	
					if( best_BHad_DRP8->BTagId(cut_tight_b_) ) BTag_lp8_BHad_tight_pass_hist->Fill(ttbar.M());
					if( !best_BHad_DRP8->BTagId(cut_loose_b_) ) BTag_lp8_BHad_loose_fail_hist->Fill(ttbar.M());
					if( !best_BHad_DRP8->BTagId(cut_medium_b_) ) BTag_lp8_BHad_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BHad_DRP8->BTagId(cut_tight_b_) ) BTag_lp8_BHad_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_BLep_DRP8 == 0) ){//reco BLep exists
				Matched_perm_BLep_pt_DRP8_hist->Fill(best_BLep_DRP8->Pt());
				Matched_perm_BLep_eta_DRP8_hist->Fill(best_BLep_DRP8->Eta());
//				++nmatched_objects; 
				if( best_BLep_DRP8 == best_BHad_DRP8 || best_BLep_DRP8 == best_WJa_DRP8 || best_BLep_DRP8 == best_WJb_DRP8 ){//reco BLep merged with something
					if( best_BLep_DRP8->BTagId(cut_loose_b_) ) BTag_lp8_BLep_loose_pass_hist->Fill(ttbar.M());
					if( best_BLep_DRP8->BTagId(cut_medium_b_) ) BTag_lp8_BLep_medium_pass_hist->Fill(ttbar.M());	
					if( best_BLep_DRP8->BTagId(cut_tight_b_) ) BTag_lp8_BLep_tight_pass_hist->Fill(ttbar.M());
					if( !best_BLep_DRP8->BTagId(cut_loose_b_) ) BTag_lp8_BLep_loose_fail_hist->Fill(ttbar.M());
					if( !best_BLep_DRP8->BTagId(cut_medium_b_) ) BTag_lp8_BLep_medium_fail_hist->Fill(ttbar.M());	
					if( !best_BLep_DRP8->BTagId(cut_tight_b_) ) BTag_lp8_BLep_tight_fail_hist->Fill(ttbar.M());
				}
			}
			if( !(best_WJa_DRP8 == 0) ){//reco WJa exists
//						++nmatched_objects; 
				Matched_perm_WJa_pt_DRP8_hist->Fill(best_WJa_DRP8->Pt());
				Matched_perm_WJa_eta_DRP8_hist->Fill(best_WJa_DRP8->Eta());
                        	if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
                        	        Matched_perm_WJa_ttbarM700_pt_DRP8_hist->Fill(best_WJa_DRP8->Pt());
                        	        Matched_perm_WJa_ttbarM700_frac_p_DRP8_hist->Fill(best_WJa_DRP8->P()/(ttbar.M()/2));
                        	}
                        	if( ttbar.M() >= 1000 ){
                        	        Matched_perm_WJa_ttbarM1000_pt_DRP8_hist->Fill(best_WJa_DRP8->Pt());
                        	        Matched_perm_WJa_ttbarM1000_frac_p_DRP8_hist->Fill(best_WJa_DRP8->P()/(ttbar.M()/2));
                        	}
			}
			if( !(best_WJb_DRP8 == 0) ){//reco WJb exists
//						++nmatched_objects; 
				Matched_perm_WJb_pt_DRP8_hist->Fill(best_WJb_DRP8->Pt());
				Matched_perm_WJb_eta_DRP8_hist->Fill(best_WJb_DRP8->Eta());
				if( ttbar.M() >= 700 && ttbar.M() < 1000 ){
				        Matched_perm_WJb_ttbarM700_pt_DRP8_hist->Fill(best_WJb_DRP8->Pt());
				        Matched_perm_WJb_ttbarM700_frac_p_DRP8_hist->Fill(best_WJb_DRP8->P()/(ttbar.M()/2));
				}
				if( ttbar.M() >= 1000 ){
				        Matched_perm_WJb_ttbarM1000_pt_DRP8_hist->Fill(best_WJb_DRP8->Pt());
				        Matched_perm_WJb_ttbarM1000_frac_p_DRP8_hist->Fill(best_WJb_DRP8->P()/(ttbar.M()/2));
				}
			}
			if( !(best_BLep_DRP8 == 0) && !(best_BLep_DRP8 == best_BHad_DRP8) && !(best_BLep_DRP8 == best_WJa_DRP8) && !(best_BLep_DRP8 == best_WJb_DRP8) ){//reco BLep not merged with anything
				if( best_BLep_DRP8->BTagId(cut_loose_b_) ) BTag_gp8_BLep_loose_pass_hist->Fill(ttbar.M());
				if( best_BLep_DRP8->BTagId(cut_medium_b_) ) BTag_gp8_BLep_medium_pass_hist->Fill(ttbar.M());	
				if( best_BLep_DRP8->BTagId(cut_tight_b_) ) BTag_gp8_BLep_tight_pass_hist->Fill(ttbar.M());
				if( !best_BLep_DRP8->BTagId(cut_loose_b_) ) BTag_gp8_BLep_loose_fail_hist->Fill(ttbar.M());
				if( !best_BLep_DRP8->BTagId(cut_medium_b_) ) BTag_gp8_BLep_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BLep_DRP8->BTagId(cut_tight_b_) ) BTag_gp8_BLep_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BLep_LepBLep_DRP8_hist->Fill(DR_LepBLep);
	
				Unmerged_BLep_massDivpt_DRP8_hist->Fill(best_BLep_DRP8->M()/best_BLep_DRP8->Pt());
				if( !(best_WJa_DRP8 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M(),(*best_BLep_DRP8+*best_WJa_DRP8).M());
				if( !(best_WJb_DRP8 == 0) ) Unmerged_BLepWJet_mass_vs_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M(),(*best_BLep_DRP8+*best_WJb_DRP8).M());
				if( !(best_WJa_DRP8 == 0) && !(best_WJb_DRP8 == 0) ){
					if( best_WJa_DRP8->M() > best_WJb_DRP8->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M(),(*best_BLep_DRP8+*best_WJa_DRP8).M());
					if( best_WJa_DRP8->M() < best_WJb_DRP8->M() ) Unmerged_BLepWJet_highest_mass_vs_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M(),(*best_BLep_DRP8+*best_WJb_DRP8).M());
				}

			}
			if( !(best_BHad_DRP8 == 0) && !(best_BHad_DRP8 == best_WJa_DRP8) && !(best_BHad_DRP8 == best_WJb_DRP8) && !(best_BHad_DRP8 == best_BLep_DRP8) ){//reco BHad not merged with anything
				if( best_BHad_DRP8->BTagId(cut_loose_b_) ) BTag_gp8_BHad_loose_pass_hist->Fill(ttbar.M());
				if( best_BHad_DRP8->BTagId(cut_medium_b_) ) BTag_gp8_BHad_medium_pass_hist->Fill(ttbar.M());	
				if( best_BHad_DRP8->BTagId(cut_tight_b_) ) BTag_gp8_BHad_tight_pass_hist->Fill(ttbar.M());
				if( !best_BHad_DRP8->BTagId(cut_loose_b_) ) BTag_gp8_BHad_loose_fail_hist->Fill(ttbar.M());
				if( !best_BHad_DRP8->BTagId(cut_medium_b_) ) BTag_gp8_BHad_medium_fail_hist->Fill(ttbar.M());	
				if( !best_BHad_DRP8->BTagId(cut_tight_b_) ) BTag_gp8_BHad_tight_fail_hist->Fill(ttbar.M());

				Unmerged_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M());
				if( !(DR_LepBLep == -1) ) Unmerged_BHad_LepBLep_DRP8_hist->Fill(DR_LepBLep);

				Unmerged_BHad_massDivpt_DRP8_hist->Fill(best_BHad_DRP8->M()/best_BHad_DRP8->Pt());
				if( !(best_WJa_DRP8 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M(),(*best_BHad_DRP8+*best_WJa_DRP8).M());
				if( !(best_WJb_DRP8 == 0) ) Unmerged_BHadWJet_mass_vs_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M(),(*best_BHad_DRP8+*best_WJb_DRP8).M());
				if( !(best_WJa_DRP8 == 0) && !(best_WJb_DRP8 == 0) ){
					if( best_WJa_DRP8->M() > best_WJb_DRP8->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M(),(*best_BHad_DRP8+*best_WJa_DRP8).M());
					if( best_WJa_DRP8->M() < best_WJb_DRP8->M() ) Unmerged_BHadWJet_highest_mass_vs_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M(),(*best_BHad_DRP8+*best_WJb_DRP8).M());
				}
			}
			if( !(best_WJa_DRP8 == 0) && !(best_WJa_DRP8 == best_BHad_DRP8) && !(best_WJa_DRP8 == best_WJb_DRP8) && !(best_WJa_DRP8 == best_BLep_DRP8) ){//reco WJa not merged with anything
				Unmerged_WJa_massDivpt_DRP8_hist->Fill(best_WJa_DRP8->M()/best_WJa_DRP8->Pt());
				Unmerged_WJa_mass_DRP8_hist->Fill(best_WJa_DRP8->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJa_LepBLep_DRP8_hist->Fill(DR_LepBLep);
			}
			if( !(best_WJb_DRP8 == 0) && !(best_WJb_DRP8 == best_BHad_DRP8) && !(best_WJb_DRP8 == best_WJa_DRP8) && !(best_WJb_DRP8 == best_BLep_DRP8) ){//reco WJb not merged with anything
				Unmerged_WJb_massDivpt_DRP8_hist->Fill(best_WJb_DRP8->M()/best_WJb_DRP8->Pt());
				Unmerged_WJb_mass_DRP8_hist->Fill(best_WJb_DRP8->M());
				if( !(DR_LepBLep == -1) ) Unmerged_WJb_LepBLep_DRP8_hist->Fill(DR_LepBLep);
			}
			if( !(best_BHad_DRP8 == 0) && !(best_WJa_DRP8 == 0) ){ //reco BHad and WJa exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJa_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJa_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJa_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());

			}
			if( !(best_BHad_DRP8 == 0) && !(best_WJb_DRP8 == 0) ){ //reco BHad and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_BHadWJb_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_BHadWJb_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_BHadWJb_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}
			if( !(best_WJa_DRP8 == 0) && !(best_WJb_DRP8 == 0) ){ //reco WJa and WJb exist
				if( object_selector_.clean_jets().size() == 3 ) All_Matched_WJaWJb_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) All_Matched_WJaWJb_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) All_Matched_WJaWJb_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP8 == best_WJa_DRP8 && !(best_BHad_DRP8 == 0) && !(best_BHad_DRP8 == best_BLep_DRP8) && !(best_BHad_DRP8 == best_WJb_DRP8) ){//only reco BHad and WJa merged
				Merged_BHadWJa_massDivpt_DRP8_hist->Fill(best_BHad_DRP8->M()/best_BHad_DRP8->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJa_perm_DRP8_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJa_perm_mass_DRP8_hist->Fill(best_BHad_DRP8->M());//merged b and wjet mass
				if( !(best_WJb_DRP8 == 0) && !(best_WJb_DRP8 == best_BHad_DRP8) && !(best_WJb_DRP8 == best_BLep_DRP8) ){//reco WJb exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M(), (*best_BHad_DRP8+*best_WJb_DRP8).M());
					Merged_BHadWJa_perm_and_WJb_mass_DRP8_hist->Fill((*best_BHad_DRP8+*best_WJb_DRP8).M());//comb invariant mass of merged jet and other wjet
					Merged_WJb_mass_DRP8_hist->Fill(best_WJb_DRP8->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}
			if( best_BHad_DRP8 == best_WJb_DRP8 && !(best_BHad_DRP8 == 0) && !(best_BHad_DRP8 == best_BLep_DRP8) && !(best_BHad_DRP8 == best_WJa_DRP8) ){//only reco BHad and WJb merged
				Merged_BHadWJb_massDivpt_DRP8_hist->Fill(best_BHad_DRP8->M()/best_BHad_DRP8->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BHadWJb_perm_DRP8_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BHadWJb_perm_mass_DRP8_hist->Fill(best_BHad_DRP8->M());//merged b and wjet mass
				if( !(best_WJa_DRP8 == 0) && !(best_WJa_DRP8 == best_BHad_DRP8) && !(best_WJa_DRP8 == best_BLep_DRP8) ){//reco WJa exists and not merged
					Merged_BHadWJet_mass_vs_BHad_mass_DRP8_hist->Fill(best_BHad_DRP8->M(), (*best_BHad_DRP8+*best_WJa_DRP8).M());
					Merged_BHadWJb_perm_and_WJa_mass_DRP8_hist->Fill((*best_BHad_DRP8+*best_WJa_DRP8).M());//comb mass of merged jet and other wjet
					Merged_WJa_mass_DRP8_hist->Fill(best_WJa_DRP8->M());
				}
				if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}
			if( best_WJa_DRP8 == best_WJb_DRP8 && !(best_WJa_DRP8 == 0) && !(best_WJa_DRP8 == best_BHad_DRP8) && !(best_WJa_DRP8 == best_BLep_DRP8) ){//only reco WJa and WJb merged
						if( object_selector_.clean_jets().size() == 3 ) Matched_WJaWJb_ptthad_3J_DRP8_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() == 4 ) Matched_WJaWJb_ptthad_4J_DRP8_hist->Fill(thad->Pt());
						if( object_selector_.clean_jets().size() > 4 ) Matched_WJaWJb_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}

			if( best_BLep_DRP8 == best_WJa_DRP8 && !(best_BLep_DRP8 == 0) && !(best_BHad_DRP8 == best_BLep_DRP8) && !(best_BLep_DRP8 == best_WJb_DRP8) ){//only reco BLep and WJa merged
				Merged_BLepWJa_massDivpt_DRP8_hist->Fill(best_BLep_DRP8->M()/best_BLep_DRP8->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJa_perm_DRP8_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJa_perm_mass_DRP8_hist->Fill(best_BLep_DRP8->M());//merged b and wjet mass
				if( !(best_WJb_DRP8 == 0) && !(best_WJb_DRP8 == best_BLep_DRP8) && !(best_WJb_DRP8 == best_BHad_DRP8) ){//reco WJb exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M(), (*best_BLep_DRP8+*best_WJb_DRP8).M());
					Merged_BLepWJa_perm_and_WJb_mass_DRP8_hist->Fill((*best_BLep_DRP8+*best_WJb_DRP8).M());//comb invariant mass of merged jet and other wjet
					//Merged_WJb_mass_DRP8_hist->Fill(best_WJb_DRP8->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJa_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJa_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJa_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}
			if( best_BLep_DRP8 == best_WJb_DRP8 && !(best_BLep_DRP8 == 0) && !(best_BHad_DRP8 == best_BLep_DRP8) && !(best_BLep_DRP8 == best_WJa_DRP8) ){//only reco BLep and WJb merged
				Merged_BLepWJb_massDivpt_DRP8_hist->Fill(best_BLep_DRP8->M()/best_BLep_DRP8->Pt());
				if( !(DR_LepBLep == -1) ) Merged_BLepWJb_perm_DRP8_LepBLep_hist->Fill(DR_LepBLep);
				Merged_BLepWJb_perm_mass_DRP8_hist->Fill(best_BLep_DRP8->M());//merged b and wjet mass
				if( !(best_WJa_DRP8 == 0) && !(best_WJa_DRP8 == best_BLep_DRP8) && !(best_WJa_DRP8 == best_BHad_DRP8) ){//reco WJa exists and not merged
					Merged_BLepWJet_mass_vs_BLep_mass_DRP8_hist->Fill(best_BLep_DRP8->M(), (*best_BLep_DRP8+*best_WJa_DRP8).M());
					Merged_BLepWJb_perm_and_WJa_mass_DRP8_hist->Fill((*best_BLep_DRP8+*best_WJa_DRP8).M());//comb mass of merged jet and other wjet
				//	Merged_WJa_mass_DRP8_hist->Fill(best_WJa_DRP8->M());
				}
				//if( object_selector_.clean_jets().size() == 3 ) Matched_BHadWJb_ptthad_3J_DRP8_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() == 4 ) Matched_BHadWJb_ptthad_4J_DRP8_hist->Fill(thad->Pt());
				//if( object_selector_.clean_jets().size() > 4 ) Matched_BHadWJb_ptthad_5PJ_DRP8_hist->Fill(thad->Pt());
			}



//			if( !(nmatched_objects == 0) ){
//				if( object_selector_.clean_jets().size() == 3 ) nMatched_objects_3J_hist->Fill(nmatched_objects);
//				if( object_selector_.clean_jets().size() == 4 ) nMatched_objects_4J_hist->Fill(nmatched_objects);
//				if( object_selector_.clean_jets().size() > 4 ) nMatched_objects_5PJ_hist->Fill(nmatched_objects);
////				cout << "Number of matched jets: " << nmatched_objects << endl;
////				cout << "" << endl;
//			}

		}
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

	static void setOptions()
	{
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
	URDriver<jet_effs> test;
	int thing = test.run();
	return thing;
}
