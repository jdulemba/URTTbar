#include <TTBarResponse.h>
#include <sstream>

#include <ttbarxsec.h>
#include <TDirectory.h>

TTBarResponse::TTBarResponse(string prefix, ttbar* an) : prefix_(prefix), an_(an), dir(0), plot1d(""), plot2d("")
{
}

TTBarResponse::~TTBarResponse() {}

void TTBarResponse::AddMatrix(string name, double default_val, const vector<double>& Mbins, const vector<double>& Tbins, string label)
{
	TDirectory* olddir = gDirectory;
	if(dir == 0)
	{
		dir = olddir->mkdir(prefix_.c_str());
	}
	dir->cd();

	plot1d.AddHist(name + "_truth", Tbins, "Truth " + label, "Events");
	plot1d.AddHist(name + "_reco", Mbins, label, "Events");
	plot2d.AddHist(name + "_matrix", Tbins, Mbins, "gen " + label, "reco " + label);
	// for(int njet = 0 ; njet < 4 ; ++njet)
	// {
	// 	stringstream hname;
	// 	hname << name << "_reco_" << njet;
	// 	plot1d.AddHist(hname.str(), Mbins, label, "Events");
	// 	hname.str("");
	// 	hname << name << "_matrix_" << njet;
	// 	plot2d.AddHist(hname.str(), Tbins, Mbins, "gen " + label, "reco " + label);
	// }
	defaults_[name] = default_val;
	values_[name] = make_pair(default_val, default_val);

	olddir->cd();
}

void TTBarResponse::FillTruth(string name, double val, double weight)
{
	weight_ = weight;
	values_[name].first = val;
	// plot1d[name + "_truth"]->Fill(val, weight);
}

void TTBarResponse::FillReco(string name, double val, double weight)
{
	weight_ = weight;
	values_[name].second =val;
	// int njet = Min(an_->reducedjets.size() - 4, size_t(3));
	// stringstream hname;
	// hname << name << "_reco_" << njet;
	// plot1d[hname.str()]->Fill(val, weight);
	// plot1d[name+"_reco"]->Fill(val, weight);
}

void TTBarResponse::FillTruthReco(string name, double tval, double rval, double weight)
{
	weight_ = weight;
	values_[name].first = tval;
	values_[name].second = rval;
	// int njet = Min(an_->reducedjets.size() - 4, size_t(3));
	// stringstream hname;
	// hname << name << "_matrix_" << njet;
	// plot2d[hname.str()]->Fill(tval, rval, weight);
	// plot2d[name+"_matrix"]->Fill(tval, rval, weight);
}

void TTBarResponse::Flush()
{
	for(auto& entry : values_)
	{
		if(defaults_[entry.first] != entry.second.first)
			plot1d[entry.first + "_truth"]->Fill(entry.second.first, weight_);

		if(defaults_[entry.first] != entry.second.second)
			plot1d[entry.first + "_reco" ]->Fill(entry.second.second, weight_);

		plot2d[entry.first +"_matrix"]->Fill(entry.second.first, entry.second.second, weight_);
	}

	//reset defaults
	for(auto& entry : defaults_)
	{
		values_[entry.first] = make_pair(entry.second, entry.second);
	}
	weight_=1;
}
