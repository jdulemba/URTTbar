#ifndef IDMET_H
#define IDMET_H
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include <TMath.h>

class IDMet : public Met, public TLorentzVector
{
    public:
        IDMet() {}
        IDMet(const Met met):
            Met(met), TLorentzVector(met.px(), met.py(), 0., TMath::Sqrt(met.px()*met.px() + met.py()*met.py()) )
    {
    }

        void setvect(double px, double py) {SetPxPyPzE(px, py, 0., TMath::Sqrt(px*px+py*py));}
        void shiftvect(double px, double py) {setvect(Px()+px, Py()+py);}
        double pxunctot() {return(TMath::Sqrt(pxunc()*pxunc() +pxuncJES()*pxuncJES() +pxuncJER()*pxuncJER()));}
        double pyunctot() {return(TMath::Sqrt(pyunc()*pyunc() +pyuncJES()*pyuncJES() +pyuncJER()*pyuncJER()));}
};
#endif
