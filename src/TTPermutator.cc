#include "Analyses/URTTbar/interface/TTPermutator.h"
#include "TMath.h"
#include <algorithm>
#include "URAnalysis/AnalysisFW/interface/Logger.h"

using namespace TMath;
using namespace std;

TTPermutator::TTPermutator(std::string cfgname):
    URSelector(cfgname),
    permutations_(),
    jets_(),
    capped_jets_() {
        configure();
    }

void TTPermutator::configure() {
    URSelector::configure();
    if(!is_configured_) {    
        addCfgParameter<size_t>("max_jets", "maximum number of jets to be considered.", 999999);
        addCfgParameter<size_t>("max_bjets", "maximum number of b jets to be considered. (0=disabled)", 0);
        addCfgParameter<float>("hardb_pt", "min pt of the hardest b-jet", 0.);
        addCfgParameter<float>("softb_pt", "min pt of the softer b-jet", 0.);
        addCfgParameter<float>("hardw_pt", "min pt of the hardest w-jet", 0.);
        addCfgParameter<float>("softw_pt", "min pt of the softer w-jet", 0.);
        addCfgParameter<std::string>("tightb", "tightest BTag working point to be applied");
        addCfgParameter<std::string>("looseb", "loosest BTag working point to be applied");

        registerConfiguration();

        cut_max_jets_ = getCfgParameter<size_t>("max_jets" );
        bjet_idx_limit_ = getCfgParameter<size_t>("max_bjets");
        cut_bjetpt_hard_ = getCfgParameter<float>( "hardb_pt" );
        cut_bjetpt_soft_ = getCfgParameter<float>( "softb_pt" );
        cut_wjetpt_hard_ = getCfgParameter<float>( "hardw_pt" );
        cut_wjetpt_soft_ = getCfgParameter<float>( "softw_pt" );
        cut_tight_b_ = IDJet::tag(getCfgParameter<std::string>("tightb"));
        cut_loose_b_ = IDJet::tag(getCfgParameter<std::string>("looseb"));
        is_configured_ = true;
    }
}

bool TTPermutator::preselection(vector<IDJet*> jets, TLorentzVector* lepton, IDMet* met, int lc, double rho, bool track) {
    reset(); //clear everything
    jets_ = jets;
    lepton_ = lepton;
    met_ = met;
    lcharge_ = lc;
    rho_ = rho;
    //keeping only the n leading jets. 
    sort(jets_.begin(), jets_.end(), [](IDJet* A, IDJet* B){return(A->Pt() > B->Pt());});
    int reducedsize = Min(jets_.size(), cut_max_jets_);
    capped_jets_.resize(reducedsize);
    copy(jets_.begin(), jets_.begin()+reducedsize, capped_jets_.begin());
    //check b-tagging conditions
    if(IDJet::id_type(cut_tight_b_) == IDJet::IDType::CSV){
        sort(capped_jets_.begin(), capped_jets_.end(), [](IDJet* A, IDJet* B){return(A->btagCSVV2() > B->btagCSVV2());});
    }
    else if(IDJet::id_type(cut_tight_b_) == IDJet::IDType::MVA){
        sort(capped_jets_.begin(), capped_jets_.end(), [](IDJet* A, IDJet* B){return(A->btagCMVA() > B->btagCMVA());});
    }
    else if(IDJet::id_type(cut_tight_b_) == IDJet::IDType::DEEPFLAVOUR){
        sort(capped_jets_.begin(), capped_jets_.end(), [](IDJet* A, IDJet* B){return(A->btagDeepB() > B->btagDeepB());});
    }
    else if(IDJet::id_type(cut_tight_b_) != IDJet::IDType::NOTSET){
        Logger::log().error() << "Don't know how to sort bjets in Permutations!" << endl;
        throw 42;
    }

    if(!capped_jets_[0]->BTagId(cut_tight_b_)) return false;
    if(tracker_ && track) tracker_->track("tight b cut");
    if(!capped_jets_[1]->BTagId(cut_loose_b_)) return false;
    if(tracker_ && track) tracker_->track("loose b cut");

    return true;
}

void TTPermutator::permutate() {
    if(capped_jets_.size() == 3) capped_jets_.push_back(NULL);
    size_t bjet_cap = (bjet_idx_limit_ > 0) ? bjet_idx_limit_  : capped_jets_.size();

    for(size_t ib1=0; ib1 < bjet_cap ; ++ib1) {
        IDJet* bjet1 = capped_jets_[ib1];
        if(!bjet1) continue; //in 3 jet case it could be NULL, skip

        for(size_t ib2=0; ib2 < bjet_cap ; ++ib2) {
            IDJet* bjet2 = capped_jets_[ib2];
            if(ib1 == ib2) continue; //same jet removal
            if(!bjet2) continue; //in 3 jet case it could be NULL, skip

            if(!(bjet1->BTagId(cut_loose_b_) && bjet2->BTagId(cut_loose_b_))) continue;
            if(!(bjet1->BTagId(cut_tight_b_) || bjet2->BTagId(cut_tight_b_))) continue;
            if(bjet1->Pt() < cut_bjetpt_hard_ && bjet2->Pt() < cut_bjetpt_hard_) continue;
            if(bjet1->Pt() < cut_bjetpt_soft_ || bjet2->Pt() < cut_bjetpt_soft_) continue;

            for(size_t iw1=0; iw1 < capped_jets_.size() ; ++iw1) {
                IDJet* wjet1 = capped_jets_[iw1]; 
                if(ib1 == iw1 || ib2 == iw1) continue; //same jet removal
                if(!wjet1) continue; //in 3 jet case it could be NULL, skip

                for(size_t iw2=0 ; iw2 < capped_jets_.size() ; ++iw2) {
                    IDJet* wjet2 = capped_jets_[iw2]; 
                    if(ib1 == iw2 || ib2 == iw2 || iw1 == iw2) continue;//same jet removal
                    //in case of three jets it SHOULD be NULL
                    if(wjet2) {
                        //disambiguation W jets
                        if(wjet2->Pt() > wjet1->Pt()) continue;
                        if(wjet1->Pt() < cut_wjetpt_hard_ || wjet2->Pt() < cut_wjetpt_soft_) continue;
                    }
                    //std::cout << "    pass ptcut" << std::endl;

                    permutations_.emplace_back(wjet1, wjet2, bjet1, bjet2, lepton_, met_, lcharge_, rho_); //added this to pass lepton charge and rho
                }
            }
        }
    }
    if(!capped_jets_[3]) capped_jets_.pop_back();
    has_run_ = true;
}

/// Joseph added for perm disc
void TTPermutator::permutate_3J(IDJet* wj1, IDJet* wj2, IDJet* bj1, IDJet* bj2, TLorentzVector* lep, IDMet* met, int lc){

    permutations_3J_.emplace_back(wj1, wj2, bj1, bj2, lep, met, lc); //added this to pass lepton charge
    permutations_3J_.emplace_back(wj1, wj2, bj2, bj1, lep, met, lc); //added this to pass lepton charge
    has_run_ = true;

}
///
