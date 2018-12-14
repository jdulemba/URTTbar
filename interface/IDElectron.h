#ifndef IDELECTRON_H
#define IDELECTRON_H
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include <atomic>
#include <map>
#include <string>

class IDElectron : public Electrons, public MCMatchable
{
    private:
        double rho_;

    public:
        IDElectron(const Electrons el, double rho=-1.) : 
            Electrons(el),
            MCMatchable(),
            rho_(rho)
    {
    }
        //FIXME
        //why this here? it is likely to break 
        //if we ever move to threaded running!
        static bool USEISO;
        enum IDS {FAIL, TIGHT_15, MEDIUM_15, LOOSE_15, VETO_15, TIGHT_15_NoECAL_Gap, NOTVETO_15, FAKES};
        static const std::map<std::string, IDElectron::IDS> id_names;

        static IDS id(std::string label);
        double rho() {return rho_;}
        double PFIsolationRho2015() const;
        double etaSC() const {return deltaEtaSC()+Eta();}

        bool LooseID25ns() const {return IPCuts() && (cutBased() == 2);}
        bool MediumID25ns() const {return IPCuts() && (cutBased() == 3);}
        bool TightID25ns() const {return IPCuts() && (cutBased() == 4);}
        bool VetoID25ns() const {return cutBased() == 1;}	
        bool FakeID() const;// {return IPCuts() && (eidCutNoIsoTight() > 0.5) ;TIGHT_15_NoECAL_Gap
        bool IPCuts() const {
            return (etaSC() < 1.479) ? (std::fabs(dxy()) < 0.05 && std::fabs(dz()) < 0.10) : (std::fabs(dxy()) < 0.10 && std::fabs(dz()) < 0.20);
        }

        bool ID(IDS idtyp);
};

#endif
