#ifndef IDJET_H
#define IDJET_H
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "Analyses/URTTbar/interface/IDMuon.h"
#include "Analyses/URTTbar/interface/IDElectron.h"
#include <TMath.h>
#include <unordered_map>
#include <string>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/BTagCalibrationStandalone.h"

class IDJet : public Jet, public MCMatchable
{
    private:
        double rndm_;
        TLorentzVector uncorr_;
    public:
        enum BTag {
            NONE, 
            DEEPCSVLOOSE, DEEPCSVMEDIUM, DEEPCSVTIGHT, 
            DEEPJETLOOSE, DEEPJETMEDIUM, DEEPJETTIGHT, 
            DEEPCSVCTAGLOOSE, DEEPCSVCTAGMEDIUM, DEEPCSVCTAGTIGHT,   
            DEEPJETCTAGLOOSE, DEEPJETCTAGMEDIUM, DEEPJETCTAGTIGHT,   
        };
        enum IDType {NOTSET, DEEPCSV, DEEPJET, DEEPCSVCTAG, DEEPJETCTAG};

        static const std::unordered_map<std::string, IDJet::BTag> tag_names;
        static IDJet::BTag tag(std::string label);

        IDJet(const Jet el, double rndm):
            Jet(el),
            MCMatchable(),
            rndm_(rndm),
            uncorr_(el)
    {
    }

        IDJet(const Jet el):
            Jet(el),
            MCMatchable(),
            rndm_(-1),
            uncorr_(el)
    {
    }

        int flavor() const {return (match()) ? match()->pdgId() : partonFlavour();}
        void update_energy(float val) {SetPtEtaPhiE(val*TMath::Sin(Theta()), Eta(), Phi(), val);}
        void resetp4() {SetPtEtaPhiE(uncorr_.Pt(), uncorr_.Eta(), uncorr_.Phi(), uncorr_.E());}
        TLorentzVector original_p4() {return uncorr_;}

        double rndm() const {return rndm_;}

        static std::string tag2string(BTag id);
        static IDType id_type(BTag id);
        static std::string id_string(BTag id);
        static BTagEntry::OperatingPoint tag_tightness(BTag id);

        bool BTagId(BTag wp) const;
        bool CTagId(BTag wp) const;
        bool TagId(BTag wp) const;

        float DeepCSVCvsLtag() const;
        float DeepCSVCvsBtag() const;
        float DeepJetCvsLtag() const;
        float DeepJetCvsBtag() const;

        float DeepCSVbDisc() const {return DeepCSVProbB() + DeepCSVProbBB();}
        float DeepJetbDisc() const {return DeepFlavourProbB() + DeepFlavourProbBB() + DeepFlavourProbLepB();}

        //bool LooseID()	{
        //    // Jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018#Recommendations_for_the_13_TeV_d
        //    //to be filled in new tree version
        //    if(TMath::Abs(Eta()) <= 2.6) {
        //        if(neutralHadronEnergyFraction() >= 0.90){return false;}
        //        if(neutralEmEnergyFraction() >= 0.90){return false;}
        //        if((chargedMultiplicity()+neutralMultiplicity()) <= 1) {return false;}
        //        if(chargedHadronEnergyFraction() <= 0.){return false;}
        //        if(chargedMultiplicity() <= 0.){return false;}
        //        return true;
        //    }
        //    else if(TMath::Abs(Eta()) <= 2.7) {
        //        if(neutralHadronEnergyFraction() >= 0.90) return false; 
        //        if(neutralEmEnergyFraction() >= 0.99) return false;
        //        if(chargedMultiplicity() <= 0) return false;
        //        return true;
        //    }
        //    else if(TMath::Abs(Eta()) <= 3.0) {
        //        if(neutralEmEnergyFraction() >= 0.99 && neutralHadronEnergyFraction() <= 0.02) return false;
        //        if(neutralMultiplicity() <= 2) return false;
        //        return true;
        //    }
        //    else {
        //        if(neutralEmEnergyFraction() >= 0.90) return false;
        //        if(neutralHadronEnergyFraction() >= 0.02) return false; 
        //        if(neutralMultiplicity() <= 10) return false;
        //        return true;
        //    }
        //}

        bool TightID()	{
            // Jet ID https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018#Recommendations_for_the_13_TeV_d
            //to be filled in new tree version
            if(TMath::Abs(Eta()) <= 2.6) {
                if(neutralHadronEnergyFraction() >= 0.90){return false;}
                if(neutralEmEnergyFraction() >= 0.90){return false;}
                if((chargedMultiplicity()+neutralMultiplicity()) <= 1) {return false;}
                if(chargedHadronEnergyFraction() <= 0.){return false;}
                if(chargedMultiplicity() <= 0.){return false;}
                return true;
            }
            else if(TMath::Abs(Eta()) <= 2.7) {
                if(neutralHadronEnergyFraction() >= 0.90) return false; 
                if(neutralEmEnergyFraction() >= 0.99) return false;
                if(chargedMultiplicity() <= 0) return false;
                return true;
            }
            else if(TMath::Abs(Eta()) <= 3.0) {
                if(neutralEmEnergyFraction() >= 0.99 && neutralHadronEnergyFraction() <= 0.02) return false;
                if(neutralMultiplicity() <= 2) return false;
                return true;
            }
            else {
                if(neutralEmEnergyFraction() >= 0.90) return false;
                if(neutralHadronEnergyFraction() >= 0.02) return false; 
                if(neutralMultiplicity() <= 10) return false;
                return true;
            }
        }

        bool Clean(const vector<IDMuon*>& muons, const vector<IDElectron*>& electrons, double distpar = 0.4) {
            for(const IDMuon* mu : muons) {
                if(DeltaR(*mu) < distpar) {
                    return false;
                }
            }
            for(const IDElectron* el : electrons)	{
                if(DeltaR(*el) < distpar)	{
                    return false;
                }
            }
            return true;
        }

        void scale(double f) { SetPtEtaPhiE(f*Pt(), Eta(), Phi(), f*E()); }

};

#endif
