#ifndef IDELECTRON_H
#define IDELECTRON_H
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include <atomic>
#include <map>
#include <string>

class IDElectron : public Electron, public MCMatchable
{
private:
	double rho_;

public:
	IDElectron(const Electron el, double rho=-1.) : 
		Electron(el),
		MCMatchable(),
		rho_(rho)
	{
	}
	//FIXME
	//why this here? it is likely to break 
  //if we ever move to threaded running!
	static bool USEISO;
	enum IDS {FAIL, TIGHT_15, MEDIUM_15, LOOSE_15, VETO_15, TIGHT_15_NoECAL_Gap};
  static const std::map<std::string, IDElectron::IDS> id_names;

  static IDS id(std::string label);
  double rho() {return rho_;}
	double PFIsolationRho2015() const;
	double etaSC() const {return TVector3(x(), y(), z()).Eta();}

  bool LooseID25ns() const;
  bool MediumID25ns() const;
  bool TightID25ns() const;
	bool VetoID25ns() const {return (eidCutVeto() > 0.5);}

	bool ID(IDS idtyp);
};

#endif
