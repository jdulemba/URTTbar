#ifndef GENOBJECT_H
#define GENOBJECT_H

#include "Analyses/URTTbar/interface/URStreamer.h"
#include "URAnalysis/AnalysisFW/interface/PDGID.h"
#include "Analyses/URTTbar/interface/LHEParticle.h"
#include <iostream>

using namespace TMath;

class GenObject : public TLorentzVector
{
    private:
        int pdgId_;
        int status_;

    public:
        GenObject(const Genparticle& gp) : TLorentzVector(gp), pdgId_(gp.pdgId()), status_(gp.status())
    {
    }
        /*GenObject(const Pst& gp) : TLorentzVector(gp), pdgId_(gp.pdgId()), status_(gp.status())
          {
          }*/
        GenObject(const Genjet& jet) : TLorentzVector(jet), pdgId_(0), status_(0)
    {
    }
        template<class T>
            GenObject(const T& t): 
                TLorentzVector(t),
                pdgId_(t.pdgId()),
                status_(t.status()) {}
        GenObject() : TLorentzVector(), pdgId_(0), status_(0)
    {
    }

        int pdgId() const {return pdgId_;}
        int status() const {return status_;}

        friend std::ostream & operator<<(std::ostream &os, const GenObject &obj);
};

class GenW { //for some reason the gen-level W it is not tracked
    public:
        enum DecayType {INVALID=0, LEPTONIC=1, HADRONIC=2};
        GenObject *first, *second;
        DecayType type = INVALID;

        GenW():
            first(0), 
            second(0),
            type(DecayType::INVALID) {}

        GenW(GenObject *first_, GenObject *second_=0) {
            if(ura::PDGID::e <= abs(first_->pdgId()) && abs(first_->pdgId()) <= ura::PDGID::nu_tau) {
                type = LEPTONIC;
                if(abs(first_->pdgId()) % 2 == 1) { //charged lepton first
                    first = first_;
                    second = second_;
                }
                else {
                    first = second_;
                    second = first_;
                }
            }
            else {
                type = HADRONIC;
                if(first_->Pt() > second_->Pt()){
                    first =first_;
                    second = second_;
                }
                else {
                    first = second_;
                    second = first_;
                }
            }
        }

        bool is_complete() {
            return first && second;
        }

        GenObject *up()   {return (first->pdgId() % 2 == 0) ? first : second;}
        GenObject *down() {return (first->pdgId() % 2 == 0) ? second : first;}

        friend std::ostream & operator<<(std::ostream &os, const GenW& w);
};

class GenTop: public GenObject {
    public:
        GenObject *b;
        GenW W;

        GenTop():
            GenObject(),
            b(0),
            W() {}

        GenTop(const GenObject& t, GenObject *b_, GenW & W_):
            GenObject(t),
            b(b_),
            W(W_) {}

        bool is_complete() {
            return b && W.is_complete();
        }

        friend std::ostream & operator<<(std::ostream &os, const GenTop& t);
};

class GenTTBar: public TLorentzVector {
    public:
        enum DecayType {INVALID=0, FULLHAD=1, SEMILEP=2, FULLLEP=3};
        GenTop top, tbar;
        DecayType type;

        GenTTBar():
            TLorentzVector(),
            top(),
            tbar(),
            type(DecayType::INVALID) {}

        GenTTBar(const GenTop &t, const GenTop &topbar):
            TLorentzVector(t+topbar),
            top(t),
            tbar(topbar) {
                if(top.W.type == GenW::DecayType::INVALID || tbar.W.type == GenW::DecayType::INVALID)
                    type = DecayType::INVALID;
                else if(top.W.type == GenW::DecayType::HADRONIC && tbar.W.type == GenW::DecayType::HADRONIC)
                    type = DecayType::FULLHAD;
                else if(top.W.type == GenW::DecayType::LEPTONIC && tbar.W.type == GenW::DecayType::LEPTONIC)
                    type = DecayType::FULLLEP;
                else 
                    type = DecayType::SEMILEP;
            }

        static GenTTBar from_collections( vector<GenObject*>& wpartons, vector<GenObject*>& charged_leps,
                vector<GenObject*>& neutral_leps, GenObject* b, GenObject* bbar,
                GenObject* top, GenObject* tbar);

        bool is_complete() {
            return top.is_complete() && tbar.is_complete();
        }

        //// Joseph added
        bool is_bhad_in_acceptance(double pt, double eta) { // check if bhad parton falls w/in pt and eta acceptance
            if( had_b()->Pt() < pt || Abs(had_b()->Eta()) > eta ) {return(false);}
            return(true);
        }

        bool is_blep_in_acceptance(double pt, double eta) { // check if blep parton falls w/in pt and eta acceptance
            if( lep_b()->Pt() < pt || Abs(lep_b()->Eta()) > eta ) {return(false);}
            return(true);
        }

        bool is_wja_in_acceptance(double pt, double eta) { // check if wja parton falls w/in pt and eta acceptance
            if( had_W()->first->Pt() < pt || Abs(had_W()->first->Eta()) > eta ) {return(false);}
            return(true);
        }

        bool is_wjb_in_acceptance(double pt, double eta) { // check if wjb parton falls w/in pt and eta acceptance
            if( had_W()->second->Pt() < pt || Abs(had_W()->second->Eta()) > eta ) {return(false);}
            return(true);
        }

        bool two_partons_in_acceptance(double pt, double eta) { // check if only 2 partons fall w/in pt and eta acceptance
            if( is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && !is_wja_in_acceptance(pt, eta) && !is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( is_bhad_in_acceptance(pt, eta) && !is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && !is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( is_bhad_in_acceptance(pt, eta) && !is_blep_in_acceptance(pt, eta) && !is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( !is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && !is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( !is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && !is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( !is_bhad_in_acceptance(pt, eta) && !is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            return(false);
        }

        bool three_partons_in_acceptance(double pt, double eta) { // check if only 3 partons fall w/in pt and eta acceptance
            if( !is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( is_bhad_in_acceptance(pt, eta) && !is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && !is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            if( is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && !is_wjb_in_acceptance(pt, eta) ) {return(true);}
            return(false);
        }

        bool all_partons_in_acceptance(double pt, double eta){ // check if bhad, blep, wja, and wjb in pt and eta acceptance
            if( is_bhad_in_acceptance(pt, eta) && is_blep_in_acceptance(pt, eta) && is_wja_in_acceptance(pt, eta) && is_wjb_in_acceptance(pt, eta) ) {return(true);}
            return(false);
        }

            // check mergings between two partons
        bool merged_lepblep_partons(double dr_){ // lepton and blep merged
            if( lepton()->DeltaR(*lep_b()) < dr_ ) return true;
            return false;
        }

        bool merged_lepbhad_partons(double dr_){ // lepton and bhad merged
            if( lepton()->DeltaR(*had_b()) < dr_ ) return true;
            return false;
        }

        bool merged_lepwja_partons(double dr_){ // lepton and wja merged
            if( lepton()->DeltaR(*had_W()->first) < dr_ ) return true;
            return false;
        }

        bool merged_lepwjb_partons(double dr_){ // lepton and wjb merged
            if( lepton()->DeltaR(*had_W()->second) < dr_ ) return true;
            return false;
        }

        bool merged_blepbhad_partons(double dr_){ // blep and bhad merged
            if( lep_b()->DeltaR(*had_b()) < dr_ ) return true;
            return false;
        }

        bool merged_blepwja_partons(double dr_){ // blep and wja merged
            if( lep_b()->DeltaR(*had_W()->first) < dr_ ) return true;
            return false;
        }

        bool merged_blepwjb_partons(double dr_){ // blep and wjb merged
            if( lep_b()->DeltaR(*had_W()->second) < dr_ ) return true;
            return false;
        }

        bool merged_bhadwja_partons(double dr_){ // bhad and wja merged
            if( had_b()->DeltaR(*had_W()->first) < dr_ ) return true;
            return false;
        }

        bool merged_bhadwjb_partons(double dr_){ // bhad and wjb merged
            if( had_b()->DeltaR(*had_W()->second) < dr_ ) return true;
            return false;
        }

        bool merged_wjawjb_partons(double dr_){ // wja and wjb merged
            if( had_W()->first->DeltaR(*had_W()->second) < dr_ ) return true;
            return false;
        }
            //
            // only two partons merged
        bool only_merged_lepblep(double dr_){ // only lep and blep merged
            if( merged_lepblep_partons(dr_) && !(merged_lepbhad_partons(dr_) || merged_lepwja_partons(dr_) || 
                merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) || merged_blepwja_partons(dr_) ||
                merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_lepbhad(double dr_){ // only lep and bhad merged
            if( merged_lepbhad_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepwja_partons(dr_) || 
                merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) || merged_blepwja_partons(dr_) ||
                merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_lepwja(double dr_){ // only lep and wja merged
            if( merged_lepwja_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) || merged_blepwja_partons(dr_) ||
                merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_lepwjb(double dr_){ // only lep and wjb merged
            if( merged_lepwjb_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_blepbhad_partons(dr_) || merged_blepwja_partons(dr_) ||
                merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_blepbhad(double dr_){ // only blep and bhad merged
            if( merged_blepbhad_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_lepwjb_partons(dr_) || merged_blepwja_partons(dr_) ||
                merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_blepwja(double dr_){ // only blep and wja merged
            if( merged_blepwja_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) ||
                merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_blepwjb(double dr_){ // only blep and wjb merged
            if( merged_blepwjb_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) ||
                merged_blepwja_partons(dr_) || merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_bhadwja(double dr_){ // only bhad and wja merged
            if( merged_bhadwja_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) ||
                merged_blepwja_partons(dr_) || merged_blepwjb_partons(dr_) || merged_bhadwjb_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_bhadwjb(double dr_){ // only bhad and wjb merged
            if( merged_bhadwjb_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) ||
                merged_blepwja_partons(dr_) || merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) ||
                merged_wjawjb_partons(dr_)) ) return true; 
            return false;
        }

        bool only_merged_wjawjb(double dr_){ // only wja and wjb merged
            if( merged_wjawjb_partons(dr_) && !(merged_lepblep_partons(dr_) || merged_lepbhad_partons(dr_) || 
                merged_lepwja_partons(dr_) || merged_lepwjb_partons(dr_) || merged_blepbhad_partons(dr_) ||
                merged_blepwja_partons(dr_) || merged_blepwjb_partons(dr_) || merged_bhadwja_partons(dr_) ||
                merged_bhadwjb_partons(dr_)) ) return true; 
            return false;
        }
            //



            // check categories of mergings
        bool resolved_had_partons(double dr_){ // hadronic partons resolved at DR = dr_
            if( !(merged_bhadwja_partons(dr_) || merged_bhadwjb_partons(dr_) || merged_wjawjb_partons(dr_)) ) return true;
            return false;
        }

        bool merged_had_partons(double dr_){ // all hadronic partons merged at DR = dr_
            if( merged_bhadwja_partons(dr_) && merged_bhadwjb_partons(dr_) && merged_wjawjb_partons(dr_) ) return true;
            return false;
        }

        bool only_merged_bhadw_partons(double dr_){ // only bhad and one W parton merged at DR = dr_
            if( only_merged_bhadwja(dr_) || only_merged_bhadwjb(dr_) ) return true;
            return false;
        }

        bool partial_hadronic_merged(double dr_){ // w partons or one wjet and bhad merged
            if( only_merged_bhadw_partons(dr_) || only_merged_wjawjb(dr_) ) return true;
            return false;
        }
        ////

        GenTop* lep_top() {
            if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLHAD) return 0;
            return (top.W.type == GenW::DecayType::LEPTONIC) ? &top : &tbar;
        }

        GenTop* had_top() {
            if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLLEP) return 0;
            return (top.W.type == GenW::DecayType::HADRONIC) ? &top : &tbar;
        }

        GenObject* lepton(){
            if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLHAD) return 0;
            return (top.W.type == GenW::DecayType::LEPTONIC) ? top.W.first : tbar.W.first;
        }

        GenObject* lep_b() {
            if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLHAD) return 0;
            return (top.W.type == GenW::DecayType::LEPTONIC) ? top.b : tbar.b;
        }

        GenObject* had_b() {
            if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLLEP) return 0;
            return (top.W.type == GenW::DecayType::HADRONIC) ? top.b : tbar.b;
        }

        GenW* had_W() {
            if(type == GenTTBar::DecayType::INVALID || type == GenTTBar::DecayType::FULLLEP) return 0;
            return (top.W.type == GenW::DecayType::HADRONIC) ? &top.W : &tbar.W;
        }

        friend std::ostream & operator<<(std::ostream &os, const GenTTBar& tt);
};

#endif
