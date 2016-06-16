#ifndef IDMUON_H
#define IDMUON_H
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include <atomic>
#include <map>
#include <string>

class IDMuon : public Muon, public MCMatchable
{
private:
	double rho_;

public:
	static bool USEISO;
	enum IDS {TIGHT_12, LOOSE_12, TIGHT_12Db, LOOSE_12Db, TIGHT_15, LOOSE_15};
  static const std::map<std::string, IDS> id_names;

	IDMuon(const Muon mu, double rho=-1);
  static IDMuon::IDS id(const std::string label);
  double rho() {return rho_;}
	double PFIsoDb();
	double RelPFIsoDb() {return PFIsoDb()/Pt();}
	double CorPFIsolation2015();
	bool ID(IDS idtyp);
  bool isTight();
  bool isLoose();
};

#endif
