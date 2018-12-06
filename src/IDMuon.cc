#include "Analyses/URTTbar/interface/IDMuon.h"
#include "TMath.h"
#include <iostream>
#include "URAnalysis/AnalysisFW/interface/Logger.h"

using namespace TMath;

const std::map<std::string, IDMuon::IDS> IDMuon::id_names = {
    {"FAIL", IDMuon::IDS::FAIL},
    {"TIGHT_12"  , IDMuon::IDS::TIGHT_12},
    {"LOOSE_12"  , IDMuon::IDS::LOOSE_12},
    {"TIGHT_12Db", IDMuon::IDS::TIGHT_12Db},
    {"LOOSE_12Db", IDMuon::IDS::LOOSE_12Db},
    {"TIGHT_15", IDMuon::IDS::TIGHT_15},
    {"LOOSE_15", IDMuon::IDS::LOOSE_15},
    {"TIGHT_15Db", IDMuon::IDS::TIGHT_15Db},
    {"TIGHT_NOISO", IDMuon::IDS::TIGHT_NOISO},
    {"LOOSE_15Db", IDMuon::IDS::LOOSE_15Db},
    {"ANTILOOSE_15Db", IDMuon::IDS::ANTILOOSE_15Db} 
};

IDMuon::IDS IDMuon::id(const std::string label) {
    try {
        return IDMuon::id_names.at(label);
    }
    catch (std::out_of_range e){
        Logger::log().error() << "Requested IDMuon ID label: " << label << " does not exist!" << std::endl;
        throw 49;
    }
}

IDMuon::IDMuon(const Muons mu, double rho): 
    Muons(mu), 
    MCMatchable(),
    rho_(rho)
{
}

//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
double IDMuon::PFIsoDb()
{
    return pfRelIso04_all(); // same as previous definition but divided by Pt now
    //return (pfChargedIso04() + TMath::Max(pfNeutralIso04() + pfPhotonIso04() - 0.5*pfPUIso04(), 0.));
}

bool IDMuon::isTight() {
    //if(!isGlobal()) return(false);
    //if(!isPF()) return(false);
    //if(chi2()/ndof() > 10.) return(false);  
    //if(validHits() <= 0) return(false);
    //if(numMatchedStations() <= 1) return(false);
    //if(TMath::Abs(dxy()) >= 0.2) return(false);
    //if(TMath::Abs(dz()) >= 0.5) return(false);
    //if(pixelHits() <= 0) return(false);
    //if(trackerLayers() <= 5) return(false);
    //return true;

    if( tightId() ) return(true);
    return false;
}

bool IDMuon::isLoose() {
    //if(!isPFcand()) return(false);
    //if(!isGlobal() && !isTracker()) return(false);
    return true;
}

//double IDMuon::CorPFIsolation2015()
//{
//    double eta = Abs(Eta());
//    double effarea = 0.;
//    if(eta < 0.8){ effarea = 0.0913;}
//    else if(eta < 1.3){ effarea = 0.0765;}
//    else if(eta < 2.0){ effarea = 0.0546;}
//    else if(eta < 2.2){ effarea = 0.0728;}
//    else if(eta < 2.5){ effarea = 0.1177;}
//
//    if(rho_ < 0.){Logger::log().error() << "Store the value of rho in the electrons to use this isolation" << endl;}
//    effarea *= Max(rho_, 0.);
//    return(chargedIso() + Max(neutralIso() + photonIso() - effarea, 0.))/Pt();
//}

bool IDMuon::ID(IDS idtyp)
{
    if(idtyp == FAIL) return false;
    else if(idtyp == TIGHT_12 || idtyp == TIGHT_12Db)
    {
        if(TMath::Abs(Eta()) > 2.4) return(false);
        ////if(!isPF()) return(false);
        ////if(!isGlobal()) return(false);
        ////if(pixelHits() <= 0) return(false);
        ////if(trackerLayers() <= 5) return(false);
        ////if(validHits() <= 0) return(false);
        ////if(numMatchedStations() <= 1) return(false);
        //if(TMath::Abs(dB()) > 0.2) return(false);
        ////if(TMath::Abs(dz()) > 0.5) return(false);
        ////if(chi2()/ndof() > 10.) return(false);

        if( !isTight() ) return(false); /// is this correct????

        if(USEISO && idtyp == TIGHT_12Db && PFIsoDb() > 0.12) return(false);
        //if(USEISO && idtyp == TIGHT_12 && (trackiso())/Pt() > 0.05) return(false);
        return(true);
    }
    else if(idtyp == LOOSE_12 || idtyp == LOOSE_12Db)
    {
        if(TMath::Abs(Eta()) > 2.4) return(false);
        //if(!isPF()) return(false);
        //if(!isGlobal() && !isTracker()) return(false);
        if( !isLoose() ) return(false);        
        if(USEISO && idtyp == LOOSE_12Db && PFIsoDb() > 0.2) return(false);
        //if(USEISO && idtyp == LOOSE_12 && (trackiso())/Pt() > 0.1) return(false);
        return(true);
    }
    //else if(idtyp == TIGHT_15) {
    //    return isTight() && ((trackiso())/Pt() < 0.1); //(PFIsoDb() < 0.15);
    //}
    //else if(idtyp == LOOSE_15) {
    //    return isLoose() && (trackiso())/Pt() < 0.1; //(PFIsoDb() < 0.25);
    //}
    else if(idtyp == IDMuon::IDS::TIGHT_15Db) {
        return isTight() && PFIsoDb() <  0.15;
    }
    else if(idtyp == TIGHT_NOISO) {
        return isTight();
    }
    else if(idtyp == IDMuon::IDS::LOOSE_15Db) {
        return isLoose() && PFIsoDb() < 0.25;
    }
    else if(idtyp == IDMuon::IDS::ANTILOOSE_15Db) {
        double relIso = PFIsoDb();
        return isTight() && 0.15 <= relIso && relIso < 0.43;
    }
    return(false);
}

bool IDMuon::USEISO = true;
