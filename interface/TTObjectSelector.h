#ifndef TTObjectSelector_h
#define TTObjectSelector_h

#include "Analyses/URTTbar/interface/URStreamer.h"
#include <list>
#include <vector>
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include "Analyses/URTTbar/interface/IDMet.h"
#include "Analyses/URTTbar/interface/JetScale.h"
#include "Analyses/URTTbar/interface/JetScaler.h"
#include <map>
#include <string>
#include "URAnalysis/AnalysisFW/interface/CutFlowTracker.h"
#include "Analyses/URTTbar/interface/systematics.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"

using namespace std;
class CutFlowTracker;

class TTObjectSelector {
    public:

        template<typename LEP>
            class LepSelector {
                public:
                    LepSelector(std::string name) {
                        URParser &parser = URParser::instance();
                        parser.addCfgParameter<std::string>(name, "id", "ID to be applied", "FAIL");
                        parser.addCfgParameter<float>      (name, "ptmin", "minimum pt", 9999.);
                        parser.addCfgParameter<float>      (name, "etamax", "maximum eta", 9999.);
                        parser.addCfgParameter<float>      (name, "etascmax", "maximum eta SC", 9999.);
                        parser.parseArguments();
                        id_ = LEP::id(
                                parser.getCfgPar<std::string>(name, "id"    ));
                        ptmin_  = parser.getCfgPar<float>(name, "ptmin" );
                        etamax_ = parser.getCfgPar<float>(name, "etamax");
                        etascmax_ = parser.getCfgPar<float>(name, "etascmax");
                    }
                    bool pass(LEP &el) {
                        return el.ID(id_) && el.Pt() > ptmin_ && fabs(el.Eta()) < etamax_ && fabs(el.etaSC()) < etascmax_;
                    }
                private:
                    typename LEP::IDS id_;
                    float ptmin_, etamax_, etascmax_;
            };

        enum EvtType {NOTSET, TIGHTMU, TIGHTEL, LOOSEMU, LOOSEEL};
        TTObjectSelector(int objsel=0); //-1 e only, 1 mu only, 0 any

        EvtType event_type() const {return evt_type_;}
        void lepton_type(int ltype) {objsel_=ltype;} //-1 e only, 1 mu only, 0 any
        int lepton_type() {return pass_lepton_;} //-1 passing e, 1 passing mu
        void reset();
        bool select(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS, bool sync=false);
        bool pass_through(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
        void allow_loose(bool val=true) {allow_loose_ = val;}

        //getters for selected objects
        vector<IDMuon*>& veto_muons() {return veto_muons_;}
        vector<IDMuon*>& loose_muons() {return loose_muons_;}
        vector<IDMuon*>& tight_muons() {return tight_muons_;}
        IDMuon* muon() {
            if(tight_muons_.size()) return tight_muons_[0];
            else if(loose_muons_.size()) return loose_muons_[0];
            return 0;
        }

        vector<IDElectron*>& veto_electrons() {return veto_electrons_;}
        vector<IDElectron*>& loose_electrons() {return loose_electrons_;}
        vector<IDElectron*>& tight_electrons() {return tight_electrons_;}
        IDElectron* electron() {
            if(tight_electrons_.size()) return tight_electrons_[0];
            else if(loose_electrons_.size()) return loose_electrons_[0];
            return 0;
        }

        vector<IDJet*>& clean_jets() {return clean_jets_;}
        IDMet * met() {return &met_;}

        TLorentzVector* lepton();
        int lepton_charge();

        void set_tracker(CutFlowTracker *t) {tracker_ = t;}
        bool pass_trig(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
        bool pass_filter(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
        bool pass_vertex(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);

    private:

        template<typename T>
            size_t nvetos(const vector<T*>& vetos, const T* val) {
                size_t ret = 0;
                for(T* i : vetos) {
                    if(i != val) ret++;
                }
                return ret;
            }

        bool select_muons(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
        bool select_electrons(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);
        bool select_jetmet(URStreamer &event, systematics::SysShifts shift=systematics::SysShifts::NOSYS);

        list<IDJet> sel_jets_;
        vector<IDJet*> clean_jets_;
        list<IDMuon> sel_muons_;
        vector<IDMuon*> veto_muons_;
        LepSelector<IDMuon> vmu_sel_; 
        vector<IDMuon*> loose_muons_;
        LepSelector<IDMuon> lmu_sel_; 
        vector<IDMuon*> tight_muons_;
        LepSelector<IDMuon> tmu_sel_; 
        list<IDElectron> sel_electrons_;
        vector<IDElectron*> veto_electrons_;
        LepSelector<IDElectron> vel_sel_; 
        vector<IDElectron*> loose_electrons_;
        LepSelector<IDElectron> lel_sel_; 
        vector<IDElectron*> tight_electrons_;
        LepSelector<IDElectron> tel_sel_; 
        IDMet met_;
        bool el_trg_, mu_trg_;
        EvtType evt_type_=NOTSET;
        //tracker
        CutFlowTracker *tracker_=0;

        //cuts
        bool is_configured_ = false;
        bool apply_jer_ = true;
        bool use_trg_ = true;
        bool use_filters_ = true;
        bool smear_met_ = true;
        bool allow_loose_ = true;
        int objsel_=0;
        int pass_lepton_ = 0;

        float cut_jet_ptmin_, cut_jet_etamax_;
        size_t cut_nminjets_;
        int cut_nmaxjets_;
        float cut_mt_;
};

#endif
