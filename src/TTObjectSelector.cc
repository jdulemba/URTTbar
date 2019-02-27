#include "Analyses/URTTbar/interface/TTObjectSelector.h"
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "TMath.h"
#include "TRandom3.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"

using namespace TMath;

TTObjectSelector::TTObjectSelector(int objsel):
    objsel_(objsel),
    sel_jets_(),
    clean_jets_(),
    sel_muons_(),
    veto_muons_(),
    vmu_sel_("veto_muons"),
    loose_muons_(),
    lmu_sel_("loose_muons"),
    tight_muons_(),
    tmu_sel_("tight_muons"),
    sel_electrons_(),
    veto_electrons_(),
    vel_sel_("veto_electrons"),
    loose_electrons_(),
    lel_sel_("loose_electrons"),
    tight_electrons_(),
    tel_sel_("tight_electrons") {
        URParser &parser = URParser::instance();

        //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
        parser.addCfgParameter<int>("event", "use_trig"   , "", 1);
        parser.addCfgParameter<int>("event", "use_filters", "", 1);
        parser.addCfgParameter<int>("event", "smear_met", "", 1);
        parser.addCfgParameter<int>("event", "allow_loose", "", 0);
        parser.addCfgParameter<float>("event", "MTCut", "", -1);

        parser.addCfgParameter<int>  ("jets", "n_min", "minimum number of jets");
        parser.addCfgParameter<int>  ("jets", "n_max", "maximum number of jets", -1);
        parser.addCfgParameter<int>  ("jets", "applyJER", "");
        parser.addCfgParameter<float>("jets", "ptmin", "minimum pt");
        parser.addCfgParameter<float>("jets", "etamax", "maximum eta");

        parser.parseArguments();
        use_trg_     = (parser.getCfgPar<int>("event", "use_trig"   ) == 1);
        use_filters_ = (parser.getCfgPar<int>("event", "use_filters") == 1);
        smear_met_   = (parser.getCfgPar<int>("event", "smear_met") == 1);
        allow_loose_ = (parser.getCfgPar<int>("event", "allow_loose") == 1);
        cut_mt_      = parser.getCfgPar<float>("event", "MTCut");

        cut_nminjets_ = parser.getCfgPar<int>  ("jets", "n_min" );
        cut_nmaxjets_ = parser.getCfgPar<int>  ("jets", "n_max" );
        apply_jer_ = (parser.getCfgPar<int>("jets", "applyJER" ) > 0.5);
        Logger::log().debug() << "Running JER: " << apply_jer_ << std::endl;
        cut_jet_ptmin_ = parser.getCfgPar<float>("jets", "ptmin" );
        cut_jet_etamax_ = parser.getCfgPar<float>("jets", "etamax");  
    }

void TTObjectSelector::reset() {
    pass_lepton_=0;
    el_trg_ = false;
    mu_trg_ = false;
    evt_type_ = EvtType::NOTSET;
    sel_jets_.clear();
    clean_jets_.clear();
    sel_muons_.clear();
    veto_muons_.clear();
    loose_muons_.clear();
    tight_muons_.clear();
    sel_electrons_.clear();
    veto_electrons_.clear();
    loose_electrons_.clear();
    tight_electrons_.clear();
}

bool TTObjectSelector::select_muons(URStreamer &event, systematics::SysShifts shift) {
    const vector<Muon>& muons = event.muons();
    for(vector<Muon>::const_iterator muon = muons.begin(); muon != muons.end(); ++muon)	{
        IDMuon mu(*muon, event.rho().value());
        bool isveto  = vmu_sel_.pass(mu);
        bool isloose = lmu_sel_.pass(mu);
        bool istight = tmu_sel_.pass(mu);
        if(isveto || isloose || istight) {
            sel_muons_.push_back(mu);
            if(isveto ) veto_muons_.push_back(&(sel_muons_.back()));
            if(isloose) loose_muons_.push_back(&(sel_muons_.back()));
            if(istight) tight_muons_.push_back(&(sel_muons_.back()));
        }
    }
    return (loose_muons_.size() > 0 || tight_muons_.size() > 0);
}

bool TTObjectSelector::select_electrons(URStreamer &event, systematics::SysShifts shift) {
    const vector<Electron>& electrons = event.electrons();
    for(vector<Electron>::const_iterator electron = electrons.begin(); electron != electrons.end(); ++electron)	{
        if(event.rho().value() < 0.) cout << "RHO: " << event.rho().value() << endl;
        IDElectron el(*electron, event.rho().value());
        bool isveto  = vel_sel_.pass(el);
        bool isloose = lel_sel_.pass(el);
        bool istight = tel_sel_.pass(el);
        if(isveto || isloose || istight) {
            sel_electrons_.push_back(el);
            if(isveto ) veto_electrons_.push_back(&(sel_electrons_.back()));
            if(isloose) loose_electrons_.push_back(&(sel_electrons_.back()));
            if(istight) tight_electrons_.push_back(&(sel_electrons_.back()));
        }
    }
    return (loose_electrons_.size() > 0 || tight_electrons_.size() > 0);
}

bool TTObjectSelector::select_jetmet(URStreamer &event, systematics::SysShifts shift) {
    //SET JETS
    const vector<Jet>& jets = event.jets();
    for(vector<Jet>::const_iterator jetit = jets.begin(); jetit != jets.end(); ++jetit)
    {
        IDJet jet(*jetit);
        //set shifts and systematics
        if(shift == systematics::SysShifts::JER_UP) jet.scale(jet.JERUp());
        else if(shift == systematics::SysShifts::JER_DW) jet.scale(jet.JERDown());
        else if(event.run == 1 && apply_jer_) jet.scale(jet.JER());

        if(shift == systematics::SysShifts::JES_DW)      jet.scale(1-jet.JESUnc());
        else if(shift == systematics::SysShifts::JES_UP) jet.scale(1+jet.JESUnc());

        if(jet.Pt() < cut_jet_ptmin_ || Abs(jet.Eta()) > cut_jet_etamax_) {continue;}
        if(!jet.Clean(veto_muons_, veto_electrons_)) {continue;}
        //if(!jet.LooseID() || !jet.Clean(veto_muons_, veto_electrons_)) {continue;}

        sel_jets_.push_back(jet);
        clean_jets_.push_back(&(sel_jets_.back()));
    }
    if( clean_jets_.size() ==0) return false;

    //SET MET
    const vector<Met>& mets = event.METs();
    if(mets.size() == 1) {
        met_ = mets[0];
        if(event.run == 1 && apply_jer_ && smear_met_)   met_.setvect(met_.pxsmear(), met_.pysmear());
        if(shift == systematics::SysShifts::JER_UP)      met_.shiftvect(   met_.pxuncJER(),    met_.pyuncJER());
        else if(shift == systematics::SysShifts::JER_DW) met_.shiftvect(-1*met_.pxuncJER(), -1*met_.pyuncJER());
        else if(shift == systematics::SysShifts::JES_DW) met_.shiftvect(   met_.pxuncJES(),    met_.pyuncJES());
        else if(shift == systematics::SysShifts::JES_UP) met_.shiftvect(-1*met_.pxuncJES(), -1*met_.pyuncJES());
        else if(shift == systematics::SysShifts::MET_DW) met_.shiftvect(met_.pxunc(), met_.pyunc());
        else if(shift == systematics::SysShifts::MET_UP) met_.shiftvect(-1*met_.pxunc(), -1*met_.pyunc());
    }
    else {
        Logger::log().error() << "met collection has more than one entry! which I should pick?" << std::endl;
    }
    return true;
}

bool TTObjectSelector::pass_through(URStreamer &event, systematics::SysShifts shift) {
    reset();
    select_muons(event, shift);
    select_electrons(event, shift);
    select_jetmet(event, shift);
    return true;
}

bool TTObjectSelector::pass_trig(URStreamer &event, systematics::SysShifts shift) {
    //Trigger!
    bool isMC = (event.run == 1);
    if(!use_trg_ && isMC) {
        el_trg_ = true;
        mu_trg_ = true;
        return true;
    }
    el_trg_ = (event.trigger().HLT_Ele27_WPTight_Gsf() == 1);
    //mu_trg_ = (event.trigger().HLT_IsoMu24() == 1 || event.trigger().HLT_IsoTkMu24() == 1);	
    //mu_trg_ = (event.trigger().HLT_IsoMu27() == 1);	
    mu_trg_ = (event.trigger().HLT_IsoMu24() == 1);	

    // cout << event.trigger().HLT_Ele27_WPLoose_Gsf() << " " <<  event.trigger().HLT_IsoMu22() << " " << event.trigger().HLT_IsoTkMu22() << endl;
    // cout << event.trigger().HLT_Ele32_eta2p1_WPTight_Gsf() << " " << event.trigger().HLT_IsoMu24() << " " << event.trigger().HLT_IsoTkMu22() << endl;
    // throw 42;
    return (el_trg_ || mu_trg_);
}

bool TTObjectSelector::pass_filter(URStreamer &event, systematics::SysShifts shift) {
    //Filters
    bool filter_answer = true;
    auto filters = event.filter();
    filter_answer &= (filters.Flag_goodVertices() == 1);
    filter_answer &= (filters.Flag_globalSuperTightHalo2016Filter() == 1);
    filter_answer &= (filters.Flag_HBHENoiseFilter() == 1); 
    filter_answer &= (filters.Flag_HBHENoiseIsoFilter() == 1); 
    filter_answer &= (filters.Flag_EcalDeadCellTriggerPrimitiveFilter() == 1);
    filter_answer &= (filters.Flag_BadPFMuonFilter() == 1);
    filter_answer &= (filters.Flag_BadChargedCandidateFilter() == 1);	
    //filter_answer &= (filters.Flag_eeBadScFilter() == 1);
    filter_answer &= (filters.Flag_ecalBadCalibFilter() == 1);
    return filter_answer;
}

bool TTObjectSelector::pass_vertex(URStreamer &event, systematics::SysShifts shift) {
    auto vtxs = event.vertexs();
    if(vtxs.size() == 0) return false; //should never happen, but also the Titanic should have never sank...
    return !vtxs[0].isFake() && vtxs[0].ndof() > 4. && fabs(vtxs[0].z()) < 24. && vtxs[0].rho() < 2.;
}

bool TTObjectSelector::select(URStreamer &event, systematics::SysShifts shift, bool sync) {
    reset();

    if(!pass_trig(event, shift)) return false;
    if(tracker_) tracker_->track("trigger");

    if(!pass_vertex(event, shift)) return false;
    if(tracker_) tracker_->track("good vtx");

    bool filter_answer = pass_filter(event, shift);
    if(use_filters_ && !filter_answer) return false;
    if(tracker_) tracker_->track("Evt filters");

    bool has_muons = select_muons(event, shift);
    bool has_electrons = select_electrons(event, shift);

    if(!(has_muons || has_electrons)) return false;
    if(tracker_) tracker_->track("has leptons");
    //cout << el_trg_ || mu_trg_ << endl;

    //1 lepton and no veto ones
    if(tight_muons_.size()) {
        if(veto_electrons_.size()) return false;
        if(nvetos(veto_muons_, tight_muons_[0])) return false;
        evt_type_ = TIGHTMU;
        pass_lepton_ = 1;
    }
    else if(tight_electrons_.size()) {
        if(veto_muons_.size()) return false;
        if(nvetos(veto_electrons_, tight_electrons_[0])) return false;
        evt_type_ = TIGHTEL;
        pass_lepton_ = -1;
    }
    else if(loose_muons_.size()) {
        if(veto_electrons_.size()) return false;
        if(loose_electrons_.size()) return false;
        if(nvetos(veto_muons_, loose_muons_[0])) return false;
        vector<IDJet*> new_jets;
        for(auto jet : clean_jets_) {
            if(jet->Clean(loose_muons_, loose_electrons_)) new_jets.push_back(jet);
        }
        clean_jets_ = new_jets;
        evt_type_ = LOOSEMU;
        pass_lepton_ = 1;
    }
    else if(loose_electrons_.size()) {
        if(veto_muons_.size()) return false;
        if(loose_muons_.size()) return false;
        if(nvetos(veto_electrons_, loose_electrons_[0])) return false;
        vector<IDJet*> new_jets;
        for(auto jet : clean_jets_) {
            if(jet->Clean(loose_muons_, loose_electrons_)) new_jets.push_back(jet);
        }
        clean_jets_ = new_jets;
        evt_type_ = LOOSEEL;
        pass_lepton_ = -1;
    }
    else {
        return false;
    }
    if(!allow_loose_ && (evt_type_ == LOOSEMU || evt_type_ == LOOSEEL)) return false;
    if(tracker_) tracker_->track("loose/tight lepton");

    //right lepton
    bool tight_lepton = (evt_type_ == TIGHTMU || evt_type_ == TIGHTEL);
    bool mutype = (evt_type_ == TIGHTMU || evt_type_ == LOOSEMU);
    if(objsel_ == -1 && mutype ) return false;
    if(objsel_ == 1  && !mutype) return false;
    if(tracker_ && tight_lepton) tracker_->track("right lepton");

    //right trigger 
    if(mutype && !mu_trg_) return false;
    if(!mutype && !el_trg_) return false;
    if(tracker_ && tight_lepton) tracker_->track("right trigger");

    if(sync && tight_lepton) {
        cout << "SYNC " << (evt_type_  == TIGHTMU) << " " << event.run << ":" << event.lumi << ":" << event.evt << endl; 
    }

    select_jetmet(event, shift);
    if(clean_jets_.size() < cut_nminjets_) return false;
    if(cut_nmaxjets_ >= 0 && clean_jets_.size() > cut_nmaxjets_) return false;
    if(tracker_ && tight_lepton) tracker_->track("has N jets");

    TLorentzVector* l = lepton();
    float mt = sqrt(pow(l->Pt() + met()->Pt(), 2) - pow(l->Px() + met()->Px(), 2) - pow(l->Py() + met()->Py(), 2));
    if(mt < cut_mt_) return false;
    if(tracker_ && tight_lepton) tracker_->track("TTObjectSelector::MT CUT");

    return true;
}

TLorentzVector* TTObjectSelector::lepton() {
    return (pass_lepton_ == 1) ? dynamic_cast<TLorentzVector*>(muon()) : dynamic_cast<TLorentzVector*>(electron());
}

int TTObjectSelector::lepton_charge() {
    return (pass_lepton_ == 1) ? muon()->charge() : electron()->charge();
}
