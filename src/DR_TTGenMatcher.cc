#include "Analyses/URTTbar/interface/DR_TTGenMatcher.h"
#include "Analyses/URTTbar/interface/GenObject.h"
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"

DR_TTGenMatcher::DR_TTGenMatcher(){
    URParser &parser = URParser::instance();
    parser.addCfgParameter<float>("DR_gen_matching", "drval", "minimum DR", 0.3);

    parser.parseArguments();

    max_dr_ = parser.getCfgPar<float>("DR_gen_matching", "drval");
}

Permutation DR_TTGenMatcher::dr_match(GenTTBar& gen_hyp, std::vector<IDJet*> &jets, TLorentzVector *lep, IDMet *met, int lepcharge){

    Permutation ret;
    Permutation empty_perm;

    ret.BHad( gen_DR_match(gen_hyp.had_b(), jets) );
    ret.BLep( gen_DR_match(gen_hyp.lep_b(), jets) );

    if( gen_hyp.had_W()->first->Pt() > gen_hyp.had_W()->second->Pt() ){
        ret.WJa( gen_DR_match(gen_hyp.had_W()->first , jets) );
        ret.WJb( gen_DR_match(gen_hyp.had_W()->second, jets) );
    }
    else {
        ret.WJa( gen_DR_match(gen_hyp.had_W()->second, jets) );
        ret.WJb( gen_DR_match(gen_hyp.had_W()->first, jets) );
    }

//    const GenObject* lepton = gen_hyp.lepton();
    ret.L( lep );
    ret.MET( met );
    ret.LepCharge( lepcharge );



//    if(fabs(lepton->pdgId()) == ura::PDGID::e){
//    IDElectron* matched = gen_DR_match(lepton, electrons);
//        ret.L(matched);
//    if(matched) ret.LepCharge(matched->charge());
//    }
//    else if(fabs(lepton->pdgId()) == ura::PDGID::mu){
//    IDMuon* matched = gen_DR_match(lepton, muons);
//        ret.L(matched);
//    if(matched) ret.LepCharge(matched->charge());
//    }
//    else { //tau
//        IDMuon* mu = gen_DR_match(lepton, muons);
//        IDElectron* el = gen_DR_match(lepton, electrons);
//        if(mu && el) {
//            vector<TLorentzVector*> leps = {mu, el};
//            TLorentzVector * lep = best_DR_match(lepton, leps);
//            ret.L(lep);
//            if(lep == mu) ret.LepCharge(mu->charge());
//            else ret.LepCharge(el->charge());
//        }
//        else if(mu) { 
//            ret.L(mu);
//            ret.LepCharge(mu->charge());
//        }
//        else if(el) {
//            ret.L(el);
//            ret.LepCharge(el->charge());
//        }       
//    }

////    if( !(ret.BHad() && ret.BLep() && ret.WJa() && ret.WJb()) ) return empty_perm; // require all objects (some should be merged)
//    if( !(ret.BHad() && ret.BLep()) ) return empty_perm; // require both b's
//    if( !(ret.WJa() || ret.WJb()) ) return empty_perm; // require at least one wjet
//    if( ret.BHad() == ret.BLep() ) return empty_perm; // both b's can't be merged
//    if( (ret.BHad() == ret.WJa() && ret.BHad() == ret.WJb()) || (ret.BLep() == ret.WJa() && ret.BLep() == ret.WJb()) ) return empty_perm; // three jets can't be merged
    return ret;

}
