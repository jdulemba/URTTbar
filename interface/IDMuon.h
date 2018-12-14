#ifndef IDMUON_H
#define IDMUON_H
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include <atomic>
#include <map>
#include <string>

class IDMuon : public Muons, public MCMatchable
{
    private:
        double rho_;

    public:
        static bool USEISO;
        enum IDS {FAIL, TIGHT_12, LOOSE_12, TIGHT_12Db, LOOSE_12Db, TIGHT_15, LOOSE_15, TIGHT_15Db, LOOSE_15Db, TIGHT_NOISO, ANTILOOSE_15Db};
        static const std::map<std::string, IDS> id_names;

        IDMuon(const Muons mu, double rho=-1);
        static IDMuon::IDS id(const std::string label);
        double rho() {return rho_;}
        double PFIsoDb() const {return pfRelIso04_all()*Pt();}
        double RelPFIsoDb() const {return pfRelIso04_all();}
        double CorPFIsolation2015();
        bool ID(IDS idtyp);
        bool isTight() const {return tightId();} // CHECK
        bool isLoose() const {return softId();} // CHECK
        double etaSC() const {return 0.;} //dummy value just to woek in templates w/ electrons
};

#endif
