#ifndef TTBARAN_H
#define TTBARAN_H
#include <iostream>
#include <list>
#include <set>
#include "AnalyzerBase.h"
#include "URStreamer.h"
#include "URDriver.h"
#include "IDMuon.h"
#include "IDElectron.h"
#include "IDJet.h"
#include "IDMet.h"
#include "GenObject.h"
#include "TTBarGenPlots.h"
#include "TTBarPlots.h"
#include "TTBarSolver.h"
#include "TTBarResponse.h"
#include "Permutation.h"
#include "BtagEff.h"
#include "JetScale.h"
#include "JetScaler.h"

using namespace std;
class PDFuncertainty;

class ttbar : public AnalyzerBase
{
	friend class TTBarGenPlots;
    friend class TTBarPlots;
    friend class TTBarResponse;

	private:
		bool isMC = false;
		map<int, set<int> >  runinfo;
		double selectionprob;
		PDFuncertainty* pdfunc;
		//Collections
		//Gen:
		bool FULLHAD;
		bool SEMILEP;
		bool FULLLEP;
		bool SEMILEPACC;
		bool SEMILEPACCLOOSE;
		bool skew_pt_distro;
		bool is_ttbar;
		bool is_data;
		list<GenObject> sgenparticles;
		vector<GenObject*> genwpartons;
		vector<GenObject*> gencls;
		vector<GenObject*> gennls;
		vector<GenObject*> genfincls;
		GenObject* gent;
		GenObject* gentbar;
		GenObject* genb;
		GenObject* genbbar;
		GenObject* genbl;
		GenObject* genbh;
		TLorentzVector gentoplep;
		TLorentzVector gentophad;

		list<Genjet> sgenjets;
		vector<Genjet*> genaddjets;

		//matched
		//vector<Jet*> recbjets;
		//vector<Jet*> recwjets;
		//Jet* recbjet;
		//Jet* recbbarjet;
		//Jet* recbhjet;
		//Jet* recbljet;
		//int nttjets;
		Permutation rightper;

		//reco
		list<IDJet> sjets;
		vector<IDJet*> cleanedjets;
		vector<IDJet*> reducedjets;
		list<IDMuon> smuons;
		vector<IDMuon*> loosemuons;
		vector<IDMuon*> tightmuons;
		list<IDElectron> selectrons;
		vector<IDElectron*> looseelectrons;
		vector<IDElectron*> mediumelectrons;
		IDMet met;

		//hists
		TH1DCollection gen1d;
		TH2DCollection gen2d;
		TH1DCollection reco1d;
		TH2DCollection reco2d;
		TH1DCollection truth1d;
		TH2DCollection truth2d;

		TTBarGenPlots ttp_genall;
		TTBarGenPlots ttp_genacc;

		TTBarPlots ttp_truth;
		TTBarPlots ttp_right;
		TTBarPlots ttp_wrong;
		TTBarPlots ttp_semi;
		TTBarPlots ttp_other;
		TTBarPlots ttp_all;

		TTBarPlots ttp_jetspos_right;
		TTBarPlots ttp_jetspos_wrong;
	//	TTBarPlots ttp_hadjets_right;
	//	TTBarPlots ttp_hadjets_wrong;
	//	TTBarPlots ttp_jets_right;
	//	TTBarPlots ttp_jets_wrong;
	//	TTBarPlots ttp_blep_right;
	//	TTBarPlots ttp_blep_wrong;
		TTBarPlots ttp_whad_right;
		TTBarPlots ttp_whad_wrong;

		TTBarPlots ttp_tlepthad_right;
		TTBarPlots ttp_tlep_right;
		TTBarPlots ttp_thad_right;
		TTBarPlots ttp_nn_right;
		TTBarPlots ttp_nsemi_right;

		TTBarResponse response;

		TTBarPlots semilep_visible_right;
		TTBarPlots semilep_right_thad;
		TTBarPlots semilep_right_tlep;
		TTBarPlots semilep_wrong;
		TTBarPlots other_tt_decay;


		BtagEff btageff;

		//ttbar solver
		TTBarSolver ttsolver;

		JetScaler jetscaler;

		//configuration
		bool DATASIM;
		bool PSEUDOTOP;
		bool BTAGMODE;
		bool JETSCALEMODE;
		bool MUONS;
		bool ELECTRONS;
		double B_TIGHT;
		double B_MEDIUM;
		int cnbtag;
		size_t cnusedjets;
		double cwjetptsoft;
		double cwjetpthard;
		double cbjetptsoft;
		double cbjetpthard;
		double cjetetamax;
		double clptmin;
		double cletamax;
		// For the fiducial loose phase space
		double cpwjetptsoft;
		double cpwjetpthard;
		double cpbjetptsoft;
		double cpbjetpthard;
		double cpjetetamax;
		double cplptmin;
		double cpletamax;
		//uncertainties
		double csigmajet;
		double csigmamet;
		double ctopptweight;
		int cfacscale;
		int crenscale;
		int crandomseed;
		double jetptmin;
		double cpjetsep;
		double cjetres;
		double cttptweight;
		bool HERWIGPP;
	
		double weight;

		//binning vectors
		vector<double> topptbins;
		vector<double> topybins;
		vector<double> ttmbins;
		vector<double> ttybins;
		vector<double> ttptbins;
		vector<double> metbins;
		vector<double> jetbins;
		vector<double> nobins;

		JetScale jetscale;

		TH1D* puhist;
		TH1D* musfhist;
		TH1D* elsfhist;
		
	public:

		ttbar(const std::string output_filename);
		~ttbar();

		//This method is called once per job at the beginning of the analysis
		//book here your histograms/tree and run every initialization needed
		virtual void begin();
		//virtual void end();
		virtual void analyze();

		void SelectGenParticles(URStreamer& event);
		void SelectRecoParticles(URStreamer& event);
		void ttanalysis(URStreamer& event);

		static void setOptions() {}
};

#endif
