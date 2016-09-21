#include "Analyses/URTTbar/interface/Permutation.h"
#include "Analyses/URTTbar/interface/IDMet.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"

using namespace std;

Permutation::Permutation(IDJet* wja, IDJet* wjb, IDJet* bjh, IDJet* bjl, TLorentzVector* lep, IDMet* met, int lcharge) :
	wja_(wja),
	wjb_(wjb),
	bjh_(bjh),
	bjl_(bjl),
	lep_(lep),
	met_(met),
  lepcharge_(lcharge)
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
}

bool operator<(const Permutation& A, const Permutation& B)
{
	return(A.Prob() < B.Prob());
}

bool operator>(const Permutation& A, const Permutation& B)
{
	return(A.Prob() > B.Prob());
}

std::ostream & operator<<(std::ostream &os, const TLorentzVector& p) {
  return os << "LV(" << p.Px() << ", " << p.Py()<< ", "<< p.Pz()<< ", "<< p.E()<< ")";
}

std::ostream & operator<<(std::ostream &os, const Permutation& p) {
  return os << "l: " << p.lep_ << ", b_l: " << p.bjl_ << ", b_h: " << p.bjh_ << ", j1: " << p.wja_ << ", j1: " << p.wjb_ << std::endl
            << ", full discr: " << p.prob_ << ", mass disc: " << p.mass_discriminant_;
}
