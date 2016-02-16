#include "Permutation.h"
#include "TTBarSolver.h"
#include "IDMet.h"
#include "Logger.h"

using namespace std;

Permutation::Permutation(IDJet* wja, IDJet* wjb, IDJet* bjh, IDJet* bjl, TLorentzVector* lep, IDMet* met) :
	wja_(wja),
	wjb_(wjb),
	bjh_(bjh),
	bjl_(bjl),
	lep_(lep),
	met_(met)
{
}

void Permutation::Reset()
{
	prob_ = numeric_limits<double>::max();
	nu_chisq_          = numeric_limits<double>::max();
	nu_discriminant_   = numeric_limits<double>::max();
	btag_discriminant_ = numeric_limits<double>::max();
	mass_discriminant_ = numeric_limits<double>::max();
	wja_ = 0;
	wjb_ = 0;
	bjh_ = 0;
	bjl_ = 0;
	lep_ = 0;
	met_ = 0;
	kinfit_ = false;
	improvedobjects.clear();
}

double Permutation::Solve(TTBarSolver& ttsolver, bool kinfit)
{
  if(!IsComplete()) {                          
    Logger::log().fatal() << "The permutation you are trying to solve is not complete!" << std::endl;
    throw 42;
  }
	kinfit_ = kinfit;
	ttsolver.Solve(bjh_, wjb_, wja_, bjl_, lep_, met_);
	nu_ = ttsolver.Nu();
	prob_ = ttsolver.Res();
	nu_chisq_          = ttsolver.NSChi2();
	nu_discriminant_   = ttsolver.NSRes();
	btag_discriminant_ = ttsolver.BTagRes();
	mass_discriminant_ = ttsolver.MassRes();

	if(kinfit_)
	{
		improvedobjects.push_back(ttsolver.Wja());
		improvedobjects.push_back(ttsolver.Wjb());
		improvedobjects.push_back(ttsolver.BHad());
		improvedobjects.push_back(ttsolver.BLep());
		improvedobjects.push_back(ttsolver.L());
		improvedobjects.push_back(ttsolver.Nu());
	}
	return(prob_);
}

bool operator<(const Permutation& A, const Permutation& B)
{
	return(A.Prob() < B.Prob());
}

bool operator>(const Permutation& A, const Permutation& B)
{
	return(A.Prob() > B.Prob());
}

std::ostream & operator<<(std::ostream &os, const Permutation& p) {
  return os << "l: " << p.lep_ << ", b_l: " << p.bjl_ << ", b_h: " << p.bjh_ << ", j1: " << p.wja_ << ", j1: " << p.wjb_ << std::endl
            << ", full discr: " << p.prob_ << ", mass disc: " << p.mass_discriminant_;
}
