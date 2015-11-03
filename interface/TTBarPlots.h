#ifndef RECOPLOTTTBAR
#define RECOPLOTTTBAR
#include <helper.h>
#include <string>
#include <iostream>
#include <TLorentzVector.h>
#include <TTBarPlotsBase.h>

using namespace std;

class Jet;
class Permutation;
class TTObjectSelector;

class TTBarPlots : public TTBarPlotsBase
{
	private:
		vector<int> jetbins_ = {-1, 0, 1, 2, 3};

	public:
		TTBarPlots(string prefix="");

		~TTBarPlots();
		void Init(const vector<double>& topptbins, const vector<double>& topybins, const vector<double>& ttmbins,
              const vector<double>& ttybins  , const vector<double>& ttptbins, const vector<double>& metbins,
              const vector<double>& jetbins  , const vector<double>& nobins);
		void Fill(Permutation& per, TTObjectSelector& objects, double weight);
};

#endif
