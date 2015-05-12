#include <TTBarPlots.h>
#include <URStreamer.h>
#include <Permutation.h>
#include <ttbarxsec.h>
#include "Logger.h"

using namespace std;
using namespace TMath;

TTBarPlots::TTBarPlots(string prefix) : TTBarPlotsBase(prefix)
{

}

TTBarPlots::~TTBarPlots()
{

}

void TTBarPlots::Init(ttbar* analysis)
{
	TTBarPlotsBase::Init(analysis);
	int ta = 60.;
	double tamin = 0.;
	double tamax = 12.;
	int tb = 60.;
	double tbmin = -20.;
	double tbmax = 10.;

	plot2d.AddHist("test_ptthad", ta, tamin, tamax, an->topptbins, "test", "p_{T}(t_{had}) (GeV)");
	plot2d.AddHist("testb_ptthad", tb, tbmin, tbmax, an->topptbins, "testb", "p_{T}(t_{had}) (GeV)");
	plot2d.AddHist("test_pttlep", ta, tamin, tamax, an->topptbins, "test", "p_{T}(t_{lep}) (GeV)");
	plot2d.AddHist("testb_pttlep", tb, tbmin, tbmax, an->topptbins, "testb", "p_{T}(t_{lep}) (GeV)");
	plot2d.AddHist("test_etathad", ta, tamin, tamax, an->topetabins, "test", "#eta(t_{had})");
	plot2d.AddHist("testb_etathad", tb, tbmin, tbmax, an->topetabins, "testb", "#eta(t_{had})");
	plot2d.AddHist("test_etatlep", ta, tamin, tamax, an->topetabins, "test", "#eta(t_{lep})");
	plot2d.AddHist("testb_etatlep", tb, tbmin, tbmax, an->topetabins, "testb", "#eta(t_{lep})");
	plot2d.AddHist("test_ttm", ta, tamin, tamax, an->ttmbins, "test", "M(tt) (GeV)");
	plot2d.AddHist("testb_ttm", tb, tbmin, tbmax, an->ttmbins, "testb", "M(tt) (GeV)");
	plot2d.AddHist("test_tty", ta, tamin, tamax, an->ttybins, "test", "y(tt)");
	plot2d.AddHist("testb_tty", tb, tbmin, tbmax, an->ttybins, "testb", "y(tt)");
	plot2d.AddHist("test_ttpt", ta, tamin, tamax, an->ttptbins, "test", "p_{T}(tt) (GeV)");
	plot2d.AddHist("testb_ttpt", tb, tbmin, tbmax, an->ttptbins, "testb", "p_{T}(tt) (GeV)");
	plot2d.AddHist("test_costhetastar", ta, tamin, tamax, 10, -1., 1., "test", "cos(#Theta*)");
	plot2d.AddHist("testb_costhetastar", tb, tbmin, tbmax, 10, -1., 1., "testb", "cos(#Theta*)");
	plot2d.AddHist("test_njet", ta, tamin, tamax, 20, 0., 20., "test", "n-jets");
	plot2d.AddHist("testb_njet", tb, tbmin, tbmax, 20, 0., 20., "testb", "n-jets");
	plot2d.AddHist("test_testb", ta, tamin, tamax, tb, tbmin, tbmax, "test", "testb");

	int full_disc_nbins = 160;
	double full_disc_min=-30, full_disc_max=20;
	plot2d.AddHist("fullDiscr_ptthad" , full_disc_nbins, full_disc_min, full_disc_max, an->topptbins, "test", "p_{T}(t_{had}) (GeV)");
	plot2d.AddHist("fullDiscr_pttlep" , full_disc_nbins, full_disc_min, full_disc_max, an->topptbins, "test", "p_{T}(t_{lep}) (GeV)");
	plot2d.AddHist("fullDiscr_etathad", full_disc_nbins, full_disc_min, full_disc_max, an->topetabins, "test", "#eta(t_{had})");
	plot2d.AddHist("fullDiscr_etatlep", full_disc_nbins, full_disc_min, full_disc_max, an->topetabins, "test", "#eta(t_{lep})");
	plot2d.AddHist("fullDiscr_ttm" 		, full_disc_nbins, full_disc_min, full_disc_max, an->ttmbins, "test", "M(tt) (GeV)");
	plot2d.AddHist("fullDiscr_tty" 		, full_disc_nbins, full_disc_min, full_disc_max, an->ttybins, "test", "y(tt)");
	plot2d.AddHist("fullDiscr_ttpt"		, full_disc_nbins, full_disc_min, full_disc_max, an->ttptbins, "test", "p_{T}(tt) (GeV)");
	plot2d.AddHist("fullDiscr_njet"   , full_disc_nbins, full_disc_min, full_disc_max, 20, 0., 20., "test", "n-jets");
	plot2d.AddHist("fullDiscr_costhetastar", full_disc_nbins, full_disc_min, full_disc_max, 10, -1., 1., "test", "cos(#Theta*)");

	plot1d.AddHist("ptthad", an->topptbins, "test", "p_{T}(t_{had}) (GeV)");
	plot1d.AddHist("pttlep", an->topptbins, "test", "p_{T}(t_{lep}) (GeV)");
	plot1d.AddHist("etathad", an->topetabins, "test", "#eta(t_{had})");
	plot1d.AddHist("etatlep", an->topetabins, "test", "#eta(t_{lep})");
	plot1d.AddHist("ttm", an->ttmbins, "test", "M(tt) (GeV)");
	plot1d.AddHist("tty", an->ttybins, "test", "y(tt)");
	plot1d.AddHist("ttpt", an->ttptbins, "test", "p_{T}(tt) (GeV)");
	plot1d.AddHist("costhetastar", 10, -1., 1., "test", "cos(#Theta*)");
	plot1d.AddHist("njet", 20, 0., 20., "test", "n-jets");
	plot1d.AddHist("fullDiscr", full_disc_nbins, full_disc_min, full_disc_max, "test", "test");
}

void TTBarPlots::Fill(Permutation& per, int lepcharge, double weight)
{
	TLorentzVector nu(per.Nu());
	TTBarPlotsBase::Fill(per.BHad(), per.WJa(), per.WJb(), per.BLep(), per.L(), &nu, lepcharge, weight);
	double test = per.MassDiscr();
	double testb = per.BDiscr();
	double full_disc = per.Prob();
	plot2d["test_testb"]->Fill(test, testb, weight);
	plot2d["test_ptthad"]->Fill(test, thad.Pt(), weight);
	plot2d["testb_ptthad"]->Fill(testb, thad.Pt(), weight);
	plot2d["test_pttlep"]->Fill(test, tlep.Pt(), weight);
	plot2d["testb_pttlep"]->Fill(testb, tlep.Pt(), weight);
	plot2d["test_etathad"]->Fill(test, Abs(thad.Eta()), weight);
	plot2d["testb_etathad"]->Fill(testb, Abs(thad.Eta()), weight);
	plot2d["test_etatlep"]->Fill(test, Abs(tlep.Eta()), weight);
	plot2d["testb_etatlep"]->Fill(testb, Abs(tlep.Eta()), weight);
	plot2d["test_ttm"]->Fill(test, tt.M(), weight);
	plot2d["testb_ttm"]->Fill(testb, tt.M(), weight);
	plot2d["test_tty"]->Fill(test, Abs(tt.Rapidity()), weight);
	plot2d["testb_tty"]->Fill(testb, Abs(tt.Rapidity()), weight);
	plot2d["test_ttpt"]->Fill(test, tt.Pt(), weight);
	plot2d["testb_ttpt"]->Fill(testb, tt.Pt(), weight);
	plot2d["test_costhetastar"]->Fill(test, tCMS.CosTheta(), weight);
	plot2d["testb_costhetastar"]->Fill(testb, tCMS.CosTheta(), weight);
	plot2d["test_njet"]->Fill(test, an->reducedjets.size(), weight);
	plot2d["testb_njet"]->Fill(test, an->reducedjets.size(), weight);

	plot2d["fullDiscr_ptthad"	]->Fill(full_disc, thad.Pt(), weight);
	plot2d["fullDiscr_pttlep"	]->Fill(full_disc, tlep.Pt(), weight);
	plot2d["fullDiscr_etathad"]->Fill(full_disc, Abs(thad.Eta()), weight);
	plot2d["fullDiscr_etatlep"]->Fill(full_disc, Abs(tlep.Eta()), weight);
	plot2d["fullDiscr_ttm"		]->Fill(full_disc, tt.M(), weight);
	plot2d["fullDiscr_tty"		]->Fill(full_disc, Abs(tt.Rapidity()), weight);
	plot2d["fullDiscr_ttpt"   ]->Fill(full_disc, tt.Pt(), weight);
	plot2d["fullDiscr_njet"   ]->Fill(full_disc, an->reducedjets.size(), weight);
	plot2d["fullDiscr_costhetastar"]->Fill(full_disc, tCMS.CosTheta(), weight);

	plot1d["ptthad"	]->Fill(thad.Pt(), weight);
	plot1d["pttlep"	]->Fill(tlep.Pt(), weight);
	plot1d["etathad"]->Fill(Abs(thad.Eta()), weight);
	plot1d["etatlep"]->Fill(Abs(tlep.Eta()), weight);
	plot1d["ttm"		]->Fill(tt.M(), weight);
	plot1d["tty"		]->Fill(Abs(tt.Rapidity()), weight);
	plot1d["ttpt"		]->Fill(tt.Pt(), weight);
	plot1d["njet"		]->Fill(an->reducedjets.size(), weight);
	plot1d["costhetastar"]->Fill(tCMS.CosTheta(), weight);
	plot1d["fullDiscr"]->Fill(full_disc, weight);
}


