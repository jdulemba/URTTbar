#ifndef PERMUTATION_H
#define PERMUTATION_H
#include <vector>
#include <limits>
#include <TLorentzVector.h>
#include "Analyses/URTTbar/interface/URStreamer.h"
#include "Analyses/URTTbar/interface/IDJet.h"
#include <iostream>
#include <stdexcept>
#include <sstream>
using namespace std;

class IDMet;

std::ostream & operator<<(std::ostream &os, const TLorentzVector& p);

class Permutation
{
	private:

    //// Joseph added for perm discriminant
        double perm_discriminant_ = numeric_limits<double>::max();
        vector<double> perm_discriminant_vec_;
    //


		double prob_ = numeric_limits<double>::max();
		double nu_chisq_          = numeric_limits<double>::max();
		double nu_discriminant_   = numeric_limits<double>::max();
		double btag_discriminant_ = numeric_limits<double>::max();
		double mass_discriminant_ = numeric_limits<double>::max();
		double qgtag_discriminant_ = numeric_limits<double>::max();
		double jratio_discriminant_ = numeric_limits<double>::max();
		IDJet* wja_ = 0;
		IDJet* wjb_ = 0;
		IDJet* bjh_ = 0;
		IDJet* bjl_ = 0;

		TLorentzVector* lep_ = 0;
        int lepcharge_=0;
        double rho_=0;
		IDMet* met_ = 0;
		TLorentzVector nu_;
		bool kinfit_ = false;
        size_t perm_jets_ = 0; //number of jets used for permutations

	public:
		Permutation() {}
//		Permutation(IDJet* wja, IDJet* wjb, IDJet* bjh, IDJet* bjl, TLorentzVector* lep, IDMet* met, int lcharge=0);
		Permutation(IDJet* wja, IDJet* wjb, IDJet* bjh, IDJet* bjl, TLorentzVector* lep, IDMet* met, int lcharge=0, double rho=-1);
        int LepCharge() {return lepcharge_;}
        double RhoVal() {return rho_;}
        void LepCharge(int c) {lepcharge_ = c;}
		void Reset();
		bool IsWHadComplete() const {return(wja_ != 0 && wjb_ != 0);}
		bool IsTHadComplete() const {return(IsWHadComplete() && bjh_ != 0);}
		bool IsTLepComplete() const {return(bjl_ != 0 && lep_ != 0 && met_ != 0);}
		bool IsComplete() const {return (IsTHadComplete() && IsTLepComplete());}
		int NumBJets() const {return((bjl_ != 0 ? 1 : 0) + (bjh_ != 0 ? 1 : 0));}
		int NumWJets() const {return((wja_ != 0 ? 1 : 0) + (wjb_ != 0 ? 1 : 0));}
		int NumTTBarJets() const {return(NumBJets() + NumWJets());}
		IDJet* WJa() const {return(wja_);}
		IDJet* WJb() const {return(wjb_);}
		IDJet* BHad() const {return(bjh_);}
		IDJet* BLep() const {return(bjl_);}
		void SwapWJets() {
			IDJet* tmp = wja_;
			wja_ = wjb_;
			wjb_ = tmp;
		}
        int permutating_jets() const {return perm_jets_;}
        void permutating_jets(int njets) { perm_jets_=njets;}
		TLorentzVector* L() const {return(lep_);}
		IDMet* MET() const {return(met_);}
		void SetMET(IDMet* met) {met_ = met;}
		void WJa(IDJet* wja){wja_=wja;}
		void WJb(IDJet* wjb){wjb_=wjb;}
		void BHad(IDJet* bjh){bjh_=bjh;}
		void BLep(IDJet* bjl){bjl_=bjl;}
		void L(TLorentzVector* lep){lep_=lep;}
		void MET(IDMet* met){met_=met;}
        
  friend std::ostream & operator<<(std::ostream &os, const Permutation& p);

		TLorentzVector Nu() const {return(nu_);}
		void Nu(TLorentzVector v) {nu_ = v;}
		double HT() const {
			double ret=0;
			if(wja_) ret += wja_->Pt();
			if(wjb_) ret += wjb_->Pt();
			if(bjh_) ret += bjh_->Pt();
			if(bjl_) ret += bjl_->Pt();
			return ret;
		}
		const TLorentzVector* NuPtr() const {return(&nu_);}
		TLorentzVector WHad() const {
			if(!WJa()) {
				stringstream err;
				err << "WJa is not defined in the permutation: " << *this;
				throw std::runtime_error(err.str());
			}
			return (WJb()) ? (*WJa() + *WJb()) : *WJa();
		}
		TLorentzVector WLep() const {return((*L() + Nu()));}
		TLorentzVector THad() const {return((WHad() + *BHad()));}
		TLorentzVector TLep() const {return((WLep() + *BLep()));}
        TLorentzVector LVect() const {return THad()+TLep();}


    //// Joseph added for perm discriminant
        double PermDiscr() const {return perm_discriminant_;}
        void PermDiscr(double val) {perm_discriminant_ = val;}

        vector<double> PermDiscr_Vec() const { return perm_discriminant_vec_;}
        void PermDiscr_Vec(vector<double> val) {perm_discriminant_vec_ = val;}

    ////

		double Prob()      const {return prob_             ;}
		double NuChisq() 	 const {return nu_chisq_         ;}
		double NuDiscr() 	 const {return nu_discriminant_  ;}
		double BDiscr()  	 const {return btag_discriminant_;}
		double MassDiscr() const {return mass_discriminant_;}
		double QGDiscr()   const {return qgtag_discriminant_;}
		double JRatioDiscr() const {return jratio_discriminant_;}
		void Prob(     double val) {prob_              = val;}
		void NuChisq(  double val) {nu_chisq_          = val;}
		void NuDiscr(  double val) {nu_discriminant_   = val;}
		void BDiscr(   double val) {btag_discriminant_ = val;}
		void MassDiscr(double val) {mass_discriminant_ = val;}
		void QGDiscr(  double val) {qgtag_discriminant_= val;} 
		void JRatioDiscr(double val) {jratio_discriminant_= val;}

    ////// Joseph added for 3 jet events 
        bool IsEmpty() const{
            if( WJa() == 0 && WJb() == 0 && BHad() == 0 && BLep() == 0 ) {return(true);}
            return(false);
        }

        int n_perm_matches() const{ // number of perm objects that have matches
            int matches = 0;
            if( BHad() ) matches++;
            if( BLep() ) matches++;
            if( WJa() ) matches++;
            if( WJb() ) matches++;
            return matches;
        }

        int unique_matches() const{ // number of perm objects that have unique matches
            int unique_matches = 0;
            if( BHad() ) unique_matches++;
            if( BLep() && ( BLep() != BHad() ) ) unique_matches++;
            if( WJa() && ( WJa() != BHad() ) && ( WJa() != BLep() ) ) unique_matches++;
            if( WJb() && ( WJb() != BHad() ) && ( WJb() != BLep() ) && ( WJb() != WJa() ) ) unique_matches++;
            return unique_matches;
        }

        // merged event categories
        bool Merged_BHadBLep() const{ // only BHad and BLep merged
            if( BHad() && (BHad() == BLep()) && (BHad() != WJa()) && (BHad() != WJb()) && (WJa() != WJb()) ){return(true);}
            return(false);
        }

        bool Merged_BHadWJa() const{ // only BHad and WJa merged
            if( BHad() && (BHad() == WJa()) && (BHad() != WJb()) && (BHad() != BLep()) && (BLep() != WJb() ) ) {return(true);}
            return(false);
        }

        bool Merged_BHadWJb() const{ // only BHad and WJb merged
            if( BHad() && (BHad() == WJb()) && (BHad() != WJa()) && (BHad() != BLep()) && (BLep() != WJa() ) ) {return(true);}
            return(false);
        }

        bool Merged_BLepWJa() const{ // only BLep and WJa merged
            if( BLep() && (BLep() == WJa()) && (BLep() != WJb()) && (BLep() != BHad()) && (BHad() != WJb() ) ) {return(true);}
            return(false);
        }

        bool Merged_BLepWJb() const{ // only BLep and WJb merged
            if( BLep() && (BLep() == WJb()) && (BLep() != WJa()) && (BLep() != BHad()) && (BHad() != WJa() ) ) {return(true);}
            return(false);
        }

        bool Merged_WJets() const{ // only WJa and WJb merged
            if( WJa() && (WJa() == WJb()) && (WJa() != BHad()) && (WJa() != BLep()) && (BHad() != BLep() ) ) {return(true);}
            return(false);
        }

        bool Merged_Event() const{ // event has a merged jet present
            if( Merged_BHadBLep() || Merged_BHadWJa() || Merged_BHadWJb() || Merged_BLepWJa() || Merged_BLepWJb() || Merged_WJets() ) {return(true);}
            return(false);
        }
        // end of merged event categories

        // compare assignments across perms
        bool AreBsFlipped(const Permutation& other) const{ // bjets have flipped positions
            return( BHad() == other.BLep() && BLep() == other.BHad() );
        }

		bool AreBsSame(const Permutation& other) const{ //bjets are selected correct and in the right position
			return(IsBHadCorrect(other) && IsBLepCorrect(other));
//			return( BLep() == other.BLep() && BHad() == other.BHad() );
		}
    /////

		bool IsValid() const
		{
			if(WJa() != 0 && WJa() == WJb()) {return(false);}
			if(BHad() != 0 && (BHad() == WJa() || BHad() == WJb())) {return(false);}
			if(BLep() != 0 && (BLep() == BHad() || BLep() == WJa() || BLep() == WJb())) {return(false);}
			return(true);
		}

		bool AreBsCorrect(const Permutation& other) const //bjets are selected correct, but not necessarily at the right position!!!!!!!!!!
		{
			return((BLep() == other.BLep() && BHad() == other.BHad()) || (BHad() == other.BLep() && BLep() == other.BHad()));
		}
		bool AreJetsCorrect(const Permutation& other) const
		{
			if(BLep() != other.BLep() && BLep() != other.BHad() && BLep() != other.WJa() && BLep() != other.WJb()){return(false);}
			if(BHad() != other.BLep() && BHad() != other.BHad() && BHad() != other.WJa() && BHad() != other.WJb()){return(false);}
			if(WJa() != other.BLep() && WJa() != other.BHad() && WJa() != other.WJa() && WJa() != other.WJb()){return(false);}
			if(WJb() != other.BLep() && WJb() != other.BHad() && WJb() != other.WJa() && WJb() != other.WJb()){return(false);}
			return(true);
		}
		int overlap(const Permutation& other) const
		{
			int olap=0;
			if(BLep() != other.BLep() && BLep() != other.BHad() && BLep() != other.WJa() && BLep() != other.WJb()){olap++;}
			if(BHad() != other.BLep() && BHad() != other.BHad() && BHad() != other.WJa() && BHad() != other.WJb()){olap++;}
			if(WJa() != other.BLep() && WJa() != other.BHad() && WJa() != other.WJa() && WJa() != other.WJb()){olap++;}
			if(WJb() != other.BLep() && WJb() != other.BHad() && WJb() != other.WJa() && WJb() != other.WJb()){olap++;}
			return olap;
		}

		bool AreHadJetsCorrect(const Permutation& other) const
		{
			if(BHad() != other.BHad() && BHad() != other.WJa() && BHad() != other.WJb()){return(false);}
			if(WJa() != other.BHad() && WJa() != other.WJa() && WJa() != other.WJb()){return(false);}
			if(WJb() != other.BHad() && WJb() != other.WJa() && WJb() != other.WJb()){return(false);}
			return(true);
		}

		bool IsBLepCorrect(const Permutation& other) const
		{
			return(BLep() == other.BLep());
		}
		bool IsBHadCorrect(const Permutation& other) const
		{
			return(BHad() == other.BHad());
		}
		bool IsWHadCorrect(const Permutation& other) const
		{
			return((WJa() == other.WJa() && WJb() == other.WJb()) || (WJa() == other.WJb() && WJb() == other.WJa()));
		}
		bool IsTHadCorrect(const Permutation& other) const
		{
			return(IsBHadCorrect(other) && IsWHadCorrect(other));
		}
		bool IsCorrect(const Permutation& other) const
		{
			return(IsTLepCorrect(other) && IsTHadCorrect(other));
		}
		bool IsJetIn(IDJet* jet)
		{
			return(jet == WJa() || jet == WJb() || jet == BHad() || jet == BLep());
		}

		bool IsLooselyCorrect(const Permutation& other) const
		{
			return(IsTLepCorrect(other) && AreHadJetsCorrect(other));
		}
		bool IsTLepCorrect(const Permutation& other) const
		{
			return(IsBLepCorrect(other) && (L() == other.L()) );
		}		
};

bool operator<(const Permutation& A, const Permutation& B);
bool operator>(const Permutation& A, const Permutation& B);

#endif
