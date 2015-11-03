#include <TTBarResponse.h>
#include <sstream>
#include "TMath.h"
#include <TDirectory.h>

TTBarResponse::TTBarResponse(string prefix) : prefix_(prefix), dir(0), plot1d(""), plot2d(""), recojets(-1), genjets(-1)
{
}

TTBarResponse::~TTBarResponse() {}

void TTBarResponse::AddMatrix(string name, const vector<double>& Mbins, const vector<double>& Tbins, string label)
{
	TDirectory* olddir = gDirectory;
	if(dir == 0)
	{
		dir = olddir->mkdir(prefix_.c_str());
	}
	dir->cd();

	plot1d.AddHist(name + "_truth", Tbins, label, "Events");
	plot1d.AddHist(name + "_reco", Mbins, label, "Events");
	plot2d.AddHist(name + "_matrix", Tbins, Mbins, "gen " + label, "reco " + label);
	for(int njet = 0 ; njet < 4 ; ++njet)
	{
		stringstream hname;
		hname << name << "_truth_" << njet;
		plot1d.AddHist(hname.str(), Tbins, label, "Events");
		hname.str("");
		hname << name << "_reco_" << njet;
		plot1d.AddHist(hname.str(), Mbins, label, "Events");
		hname.str("");
		hname << name << "_matrix_" << njet;
		plot2d.AddHist(hname.str(), Tbins, Mbins, "gen " + label, "reco " + label);
	}

	olddir->cd();
	values_[name].first = Tbins.front() - 1.;
	values_[name].second = Mbins.front() - 1.;
}

void TTBarResponse::FillTruth(string name, double val, size_t njets, double weight)
{
	weight_ = weight;
	values_[name].first = val;
	genjets = TMath::Min(njets, size_t(3));
}

void TTBarResponse::FillReco(string name, double val, size_t njets, double weight)
{
	weight_ = weight;
	values_[name].second =val;
	recojets = TMath::Min(njets - 4, size_t(3));
}

void TTBarResponse::FillTruthReco(string name, double tval, double rval, double weight)
{
	weight_ = weight;
	values_[name].first = tval;
	values_[name].second = rval;
}

void TTBarResponse::Flush()
{
	stringstream hname;
	for(auto& entry : values_)
	{
		if(genjets != -1)
		{
			plot1d[entry.first + "_truth"]->Fill(entry.second.first, weight_);
			hname.str("");
			hname << entry.first << "_truth_" << genjets;
			plot1d[hname.str()]->Fill(entry.second.first, weight_);
		}

		if(recojets != -1)
		{
			plot1d[entry.first + "_reco" ]->Fill(entry.second.second, weight_);
			hname.str("");
			hname << entry.first << "_reco_" << recojets;
			plot1d[hname.str()]->Fill(entry.second.second, weight_);

			hname.str("");
			hname << entry.first << "_matrix_" << recojets;
			plot2d[hname.str()]->Fill(entry.second.first, entry.second.second, weight_);
		}

		plot2d[entry.first +"_matrix"]->Fill(entry.second.first, entry.second.second, weight_);
	}

	//reset defaults
	for(auto& entry : values_)
	{
		entry.second.first = plot2d[entry.first + "_matrix"]->GetXaxis()->GetXmin()-1.;
		entry.second.second = plot2d[entry.first + "_matrix"]->GetYaxis()->GetXmin()-1.;
	}
	genjets = -1;
	recojets = -1;
	weight_=1;
}
