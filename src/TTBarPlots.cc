#include <TTBarPlots.h>
#include <URStreamer.h>
#include <Permutation.h>
#include "Logger.h"
#include "TTObjectSelector.h"

#include <sstream>
#include <iostream>

using namespace std;
using namespace TMath;

TTBarPlots::TTBarPlots(string prefix) : TTBarPlotsBase(prefix)
{

}

TTBarPlots::~TTBarPlots()
{

}

void TTBarPlots::Init(const vector<double>& topptbins, const vector<double>& topybins, const vector<double>& ttmbins,
                      const vector<double>& ttybins  , const vector<double>& ttptbins, const vector<double>& metbins,
                      const vector<double>& jetbins  , const vector<double>& nobins)
{
	TTBarPlotsBase::Init();
	int ta = 40.;
	double tamin = 6.;
	double tamax = 18.;
	int tb = 40.;
	//double tbmin = -20.;
	//double tbmax = 10.;
	double tbmin = 0.;
	double tbmax = 80.;
	
	plot1d.AddHist("MET", 500, 0, 2000, "MET", "Events");
	plot1d.AddHist("weight", 200, -5, 5, "weight", "Events");
	plot1d.AddHist("njets", 15, 0, 15, "n-jets", "Events");
	plot1d.AddHist("nvtx", 50, 0, 50, "n-vertices", "Events");
	plot1d.AddHist("rho", 30, 0, 30, "rho", "Events");
	plot1d.AddHist("trigger", 9, 0, 9, "trigger", "Events");
  TAxis *xax = plot1d["trigger"]->GetXaxis();
  xax->SetBinLabel(1, "HLT_IsoMu24_eta2p1");
  xax->SetBinLabel(2, "HLT_IsoMu20_eta2p1");
  xax->SetBinLabel(3, "HLT_IsoMu22");
  xax->SetBinLabel(4, "HLT_IsoMu20");
  xax->SetBinLabel(5, "HLT_IsoMu27");
	plot1d.AddHist("ptaddjets", 200, 0, 400, "p_{T}(add. jets) [GeV]", "Events");
	plot1d.AddHist("DRminWjets", 200, 0, 10, "#DeltaR_{min W-jet}", "Events");
	plot1d.AddHist("DRminbjets", 200, 0, 10, "#DeltaR_{min b-jet}", "Events");
	plot1d.AddHist("DPhiMET_Nu", 100, 0, 3, "#Delta#Phi(#nu, MET)", "Events");
	plot2d.AddHist("METvsDPhiMET_Nu", 120, 0, 1200, 100, 0, 3, "MET [GeV]", "#Delta#Phi(#nu, MET)");
	plot2d.AddHist("METvsChi", 120, 0, 1200, 25, 0., 100., "MET [GeV]", "#chi");
	plot1d.AddHist("Mt_W", 500, 0, 500, "M_{t}(W) [GeV]", "Events");
	plot2d.AddHist("METunc", 100, 0, 0.5, 100, 0., .5, "#sigma(MET_{x})/MET_{x}", "#sigma(MET_{y})/MET_{y}");
	plot1d.AddHist("massDiscr", ta, tamin, tamax, "#lambda_{m}", "p_{T}(t_{had}) [GeV]");
	for(int jn : jetbins_)
	{
		stringstream jb;
		if(jn != -1) jb << jn << "_";
		plot2d.AddHist("massDiscr_"+jb.str()+"thadpt", ta, tamin, tamax, topptbins, "#lambda_{m}", "p_{T}(t_{had}) [GeV]");
		plot2d.AddHist("nsDiscr_"+jb.str()+"thadpt", tb, tbmin, tbmax, topptbins, "D_{min}", "p_{T}(t_{had}) [GeV]");
		plot2d.AddHist("massDiscr_"+jb.str()+"nobin", ta, tamin, tamax, nobins, "#lambda_{m}", "all");
		plot2d.AddHist("nsDiscr_"+jb.str()+"nobin", tb, tbmin, tbmax, nobins, "D_{min}", "all");
		plot2d.AddHist("massDiscr_"+jb.str()+"tleppt", ta, tamin, tamax, topptbins, "#lambda_{m}", "p_{T}(t_{lep}) [GeV]");
		plot2d.AddHist("nsDiscr_"+jb.str()+"tleppt", tb, tbmin, tbmax, topptbins, "D_{min}", "p_{T}(t_{lep}) [GeV]");
		plot2d.AddHist("massDiscr_"+jb.str()+"thady", ta, tamin, tamax, topybins, "#lambda_{m}", "|y(t_{had})|");
		plot2d.AddHist("nsDiscr_"+jb.str()+"thady", tb, tbmin, tbmax, topybins, "D_{min}", "|y(t_{had})|");
		plot2d.AddHist("massDiscr_"+jb.str()+"tlepy", ta, tamin, tamax, topybins, "#lambda_{m}", "|y(t_{lep})|");
		plot2d.AddHist("nsDiscr_"+jb.str()+"tlepy", tb, tbmin, tbmax, topybins, "D_{min}", "|y(t_{lep})|");
		plot2d.AddHist("massDiscr_"+jb.str()+"ttm", ta, tamin, tamax, ttmbins, "#lambda_{m}", "M(tt) [GeV]");
		plot2d.AddHist("nsDiscr_"+jb.str()+"ttm", tb, tbmin, tbmax, ttmbins, "D_{min}", "M(tt) [GeV]");
		plot2d.AddHist("massDiscr_"+jb.str()+"tty", ta, tamin, tamax, ttybins, "#lambda_{m}", "y(tt)");
		plot2d.AddHist("nsDiscr_"+jb.str()+"tty", tb, tbmin, tbmax, ttybins, "D_{min}", "y(tt)");
		plot2d.AddHist("massDiscr_"+jb.str()+"ttpt", ta, tamin, tamax, ttptbins, "#lambda_{m}", "p_{T}(tt) [GeV]");
		plot2d.AddHist("nsDiscr_"+jb.str()+"ttpt", tb, tbmin, tbmax, ttptbins, "D_{min}", "p_{T}(tt) [GeV]");
		plot2d.AddHist("massDiscr_"+jb.str()+"costhetastar", ta, tamin, tamax, 10, -1., 1., "#lambda_{m}", "cos(#Theta*)");
		plot2d.AddHist("nsDiscr_"+jb.str()+"costhetastar", tb, tbmin, tbmax, 10, -1., 1., "D_{min}", "cos(#Theta*)");
		plot2d.AddHist("massDiscr_"+jb.str()+"njet", ta, tamin, tamax, jetbins, "#lambda_{m}", "n-jets");
		plot2d.AddHist("nsDiscr_"+jb.str()+"njet", tb, tbmin, tbmax, jetbins, "D_{min}", "n-jets");
		plot2d.AddHist("massDiscr_"+jb.str()+"met", ta, tamin, tamax, metbins, "#lambda_{m}", "MET [GeV]");
		plot2d.AddHist("nsDiscr_"+jb.str()+"met", tb, tbmin, tbmax, metbins, "D_{min}", "MET [GeV]");
		plot2d.AddHist("massDiscr_"+jb.str()+"nsDiscr", ta, tamin, tamax, tb, tbmin, tbmax, "#lambda_{m}", "nsDiscr");
	}
}

void TTBarPlots::Fill(Permutation& per, TTObjectSelector& objects, URStreamer &evt, double weight)
{
	TLorentzVector nu(per.Nu());
	TTBarPlotsBase::Fill(per.BHad(), per.WJa(), per.WJb(), per.BLep(), per.L(), &nu, objects.lepton_charge(), weight);
	double massDiscr = per.MassDiscr();
	//double testb = per.BDiscr();
	double nsDiscr = Sqrt(per.NuChisq());
	//double testb = (*per.BLep() + *per.L()).Mt(); 
	//double full_disc = per.Prob();
	
	if(massDiscr == numeric_limits<double>::max()) {massDiscr = 0; nsDiscr = 0;}
	plot1d["MET"]->Fill(objects.met()->Pt(), weight);
	plot1d["njets"]->Fill(objects.clean_jets().size()-4, weight);
	plot1d["weight"]->Fill(weight);
	double drminw = 100.;
	double drminb = 100.;
	for(size_t j = 0 ; j < objects.clean_jets().size() ; ++j)
	{
		if(per.IsJetIn(objects.clean_jets()[j])) continue;
		if(drminw > per.WJa()->DeltaR(*objects.clean_jets()[j])) {drminw = per.WJa()->DeltaR(*objects.clean_jets()[j]);}
		if(drminw > per.WJb()->DeltaR(*objects.clean_jets()[j])) {drminw = per.WJb()->DeltaR(*objects.clean_jets()[j]);}
		if(drminb > per.BHad()->DeltaR(*objects.clean_jets()[j])) {drminb = per.BHad()->DeltaR(*objects.clean_jets()[j]);}
		if(drminb > per.BLep()->DeltaR(*objects.clean_jets()[j])) {drminb = per.BLep()->DeltaR(*objects.clean_jets()[j]);}
		plot1d["ptaddjets"]->Fill(objects.clean_jets()[j]->Pt(), weight);
	}
	plot1d["nvtx"]->Fill(evt.rho().value(), weight);
	plot1d["rho" ]->Fill(evt.vertexs().size(), weight);

  if(evt.trigger().HLT_IsoMu24_eta2p1() == 1)  plot1d["trigger"]->Fill(0.5, weight);
  if(evt.trigger().HLT_IsoMu20_eta2p1() == 1)  plot1d["trigger"]->Fill(1.5, weight);
  if(evt.trigger().HLT_IsoMu22() == 1)  plot1d["trigger"]->Fill(2.5, weight);
  if(evt.trigger().HLT_IsoMu20() == 1)  plot1d["trigger"]->Fill(3.5, weight);
  if(evt.trigger().HLT_IsoMu27() == 1)  plot1d["trigger"]->Fill(4.5, weight);


	plot1d["massDiscr"]->Fill(massDiscr, weight);
	plot1d["DRminWjets"]->Fill(drminw, weight);
	plot1d["DRminbjets"]->Fill(drminb, weight);
	plot2d["METvsChi"]->Fill(objects.met()->Pt(), nsDiscr, weight);
	plot1d["DPhiMET_Nu"]->Fill(Abs(nu.DeltaPhi(*objects.met())), weight);
	plot2d["METvsDPhiMET_Nu"]->Fill(objects.met()->Pt(), Abs(nu.DeltaPhi(*objects.met())), weight);
	double Mt_W = Sqrt(2.*objects.met()->Pt()*per.L()->Pt()-2.*(objects.met()->Px()*per.L()->Px() + objects.met()->Py()*per.L()->Py()));
	plot1d["Mt_W"]->Fill(Mt_W, weight);
	plot2d["METunc"]->Fill(objects.met()->pxunctot()/objects.met()->Px(), objects.met()->pyunctot()/objects.met()->Py(), weight);
	for(int jn : jetbins_)
	{
		stringstream jb;
		if(jn != -1) jb << jn << "_";
		if((jn == -1) || (per.permutating_jets() - 4 == size_t(jn)) || (jn == jetbins_.back() && per.permutating_jets() - 4 > size_t(jn)))
		{
			plot2d["massDiscr_"+jb.str()+"thadpt"]->Fill(massDiscr, thad.Pt(), weight);
			plot2d["nsDiscr_"+jb.str()+"thadpt"]->Fill(nsDiscr, thad.Pt(), weight);
			plot2d["massDiscr_"+jb.str()+"nobin"]->Fill(massDiscr, thad.Pt(), weight);
			plot2d["nsDiscr_"+jb.str()+"nobin"]->Fill(nsDiscr, thad.Pt(), weight);
			plot2d["massDiscr_"+jb.str()+"tleppt"]->Fill(massDiscr, tlep.Pt(), weight);
			plot2d["nsDiscr_"+jb.str()+"tleppt"]->Fill(nsDiscr, tlep.Pt(), weight);
			plot2d["massDiscr_"+jb.str()+"thady"]->Fill(massDiscr, Abs(thad.Rapidity()), weight);
			plot2d["nsDiscr_"+jb.str()+"thady"]->Fill(nsDiscr, Abs(thad.Rapidity()), weight);
			plot2d["massDiscr_"+jb.str()+"tlepy"]->Fill(massDiscr, Abs(tlep.Rapidity()), weight);
			plot2d["nsDiscr_"+jb.str()+"tlepy"]->Fill(nsDiscr, Abs(tlep.Rapidity()), weight);
			plot2d["massDiscr_"+jb.str()+"ttm"]->Fill(massDiscr, tt.M(), weight);
			plot2d["nsDiscr_"+jb.str()+"ttm"]->Fill(nsDiscr, tt.M(), weight);
			plot2d["massDiscr_"+jb.str()+"tty"]->Fill(massDiscr, Abs(tt.Rapidity()), weight);
			plot2d["nsDiscr_"+jb.str()+"tty"]->Fill(nsDiscr, Abs(tt.Rapidity()), weight);
			plot2d["massDiscr_"+jb.str()+"ttpt"]->Fill(massDiscr, tt.Pt(), weight);
			plot2d["nsDiscr_"+jb.str()+"ttpt"]->Fill(nsDiscr, tt.Pt(), weight);
			plot2d["massDiscr_"+jb.str()+"costhetastar"]->Fill(massDiscr, tCMS.CosTheta(), weight);
			plot2d["nsDiscr_"+jb.str()+"costhetastar"]->Fill(nsDiscr, tCMS.CosTheta(), weight);
			plot2d["massDiscr_"+jb.str()+"njet"]->Fill(massDiscr, objects.clean_jets().size()-4, weight);
			plot2d["nsDiscr_"+jb.str()+"njet"]->Fill(nsDiscr, objects.clean_jets().size()-4, weight);
			plot2d["massDiscr_"+jb.str()+"met"]->Fill(massDiscr, objects.met()->Pt(), weight);
			plot2d["nsDiscr_"+jb.str()+"met"]->Fill(nsDiscr, objects.met()->Pt(), weight);
			plot2d["massDiscr_"+jb.str()+"nsDiscr"]->Fill(massDiscr, nsDiscr, weight);
		}
	}
}


