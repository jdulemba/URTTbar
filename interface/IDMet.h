#ifndef IDMET_H
#define IDMET_H
#include "Analyses/URTTbar/interface/MCMatchable.h"
#include "Analyses/URTTbar/interface/URStreamer.h"
#include <TMath.h>

class IDMet : public Met//, public TLorentzVector
{
public:
	IDMet() {}
	IDMet(const Met met):
		Met(met)//, TLorentzVector(met.Pt(), 0., met.Phi(), 0. )
		{
		}

    double Px() const {return Pt()*TMath::Cos(Phi());}
    double Py() const {return Pt()*TMath::Sin(Phi());}
	//void setvect(double px, double py) {SetPxPyPzE(px, py, 0., TMath::Sqrt(px*px+py*py));}
	//void shiftvect(double px, double py) {setvect(Px()+px, Py()+py);}
	//double pxunctot() {return(TMath::Sqrt(pxunc()*pxunc() +pxuncJES()*pxuncJES() +pxuncJER()*pxuncJER()));}
	//double pyunctot() {return(TMath::Sqrt(pyunc()*pyunc() +pyuncJES()*pyuncJES() +pyuncJER()*pyuncJER()));}
};
#endif
