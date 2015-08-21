#ifndef IDMUON_H
#define IDMUON_H
#include "URStreamer.h"
#include "MCMatchable.h"
#include <atomic>

class IDMuon : public Muon, public MCMatchable
{
private:
	double rho_;

public:
	static bool USEISO;
	enum IDS {TIGHT_12, LOOSE_12, TIGHT_12Db, LOOSE_12Db};
	IDMuon(const Muon mu, double rho=-1);
	double PFIsoDb();
	double CorPFIsolation2015();
	bool ID(IDS idtyp);
};
#endif
