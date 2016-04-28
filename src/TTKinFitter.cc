#include "TTKinFitter.h"
#include "TMatrixF.h"
#include "TLorentzVector.h"
#include <cmath>
#include "URParser.h"
#include "TMath.h"
#include "Logger.h"

// Kinematic fit using lagrangian multipliers method
// follows notation described in http://www.phys.ufl.edu/~avery/fitting/kinfit_talk1.pdf
// and AN 2005/025
// Assumes massless jets

TTKinFitter::TTKinFitter():
  jet_scaler_(JetScaler::instance()) {
    URParser &parser = URParser::instance();
    //parser.addCfgParameter(const std::string group, const std::string parameterName, const std::string description, T def_value);
    parser.addCfgParameter<int>  ("kinfit", "maxiter", "maximum number of iterations");
    parser.addCfgParameter<float>("kinfit", "precision", "early termination/convergence condition");
    
    parser.parseArguments();

    maxiter_   = parser.getCfgPar<int  >("kinfit", "maxiter");
    precision_ = parser.getCfgPar<float>("kinfit", "precision");
}

inline TMatrix TTKinFitter::mass_contraints_derivative(double e1, double e2, double eb, double f12, double f1b, double f2b) {
  //computes mass constraints derivative matrix for given point
  TMatrix retval(2, 3);

  float wdelta = 2*e1*e2*f12 - W_mass2_;
  float w_sign = 0;
  if(wdelta > 0) w_sign = 1;
  else if(wdelta < 0) w_sign = -1;

  float tdelta = 2*e1*e2*f12 + 2*e1*eb*f1b + 2*e2*eb*f2b - top_mass2_;
  float t_sign = 0;
  if(tdelta > 0) t_sign = 1;
  else if(tdelta < 0) t_sign = -1;

  //W mass constraint
  retval(0, 0) = w_sign*2*e2*f12;
  retval(0, 1) = w_sign*2*e1*f12;
  retval(0, 2) = 0;

  //top mass constraint
  retval(1, 0) = t_sign*(2*e2*f12 + 2*eb*f1b);
  retval(1, 1) = t_sign*(2*e1*f12 + 2*eb*f2b);
  retval(1, 2) = t_sign*(2*e1*f1b + 2*e2*f2b);

  return retval;
}

inline TMatrix TTKinFitter::mass_constraint(double e1, double e2, double eb, double f12, double f1b, double f2b) {
  // cout << "masses " << sqrt(2*e1*e2*f12) << ", " << sqrt(2*e1*e2*f12 + 2*e1*e2*f1b + 2*e2*eb*f2b) << endl;
  // cout << "energies " << e1 << " " << e2 << " " << eb << endl;
  TMatrix retval(2, 1);
  retval(0, 0) = fabs(2*e1*e2*f12 - W_mass2_);
  retval(1, 0) = fabs(2*e1*e2*f12 + 2*e1*eb*f1b + 2*e2*eb*f2b - top_mass2_);

  // cout << "constraints " << retval(0,0) << " " << retval(1,0) << endl;
  // cout << "-----------------------------------------------" << endl << endl;
  return retval;
}

TTKinFitter::KFResult TTKinFitter::fit(Permutation &perm) {
  //compute opening angles (considered constant) and initial energies
  double f12 = 1 - perm.WJa()->Vect().Unit().Dot(perm.WJb()->Vect().Unit());
  double f1b = 1 - perm.WJa()->Vect().Unit().Dot(perm.BHad()->Vect().Unit());
  double f2b = 1 - perm.WJb()->Vect().Unit().Dot(perm.BHad()->Vect().Unit());

  IDJet* jets[3] = {perm.WJa(), perm.WJb(), perm.BHad()}; 
  TMatrix nrgs_in(3, 1);
  TMatrix cov(3, 3); //covariance matrix
  TMatrix cov_inv(3, 3);
  for(int i=0; i<3; i++) {
    nrgs_in(i, 0) = jets[i]->E();
    float res = pow(jet_scaler_.resolution(*jets[i]), 2);
    cov(i, i) = res;
    cov_inv(i, i) = 1./res;
  }

  TMatrix nrgs_cur = nrgs_in; //current energies  
  
  int nit=0;
  bool converged=false;  
  bool failed=false;
  for(; nit<maxiter_; nit++) {
    //mass constraint    
    TMatrix D = mass_contraints_derivative(
      nrgs_cur(0,0), nrgs_cur(1,0), nrgs_cur(2,0), 
      f12, f1b, f2b
      );
    TMatrix d = mass_constraint(
      nrgs_cur(0,0), nrgs_cur(1,0), nrgs_cur(2,0), 
      f12, f1b, f2b
      );     

    //support matrices
    TMatrix DT(TMatrix::kTransposed, D);
    TMatrix VD = D*cov*DT;
    TMatrix VD_inv = VD;
    VD.InvertFast();
    TMatrix delta = nrgs_cur - nrgs_in;
    TMatrix constraint = D*delta - d;

    //results
    TMatrix delta_new = cov*DT*VD*constraint;
    TMatrix nrgs_new = nrgs_in+delta_new; 

    //check for faliure
    for(int i=0; i<3; i++) failed =  failed || TMath::IsNaN(delta_new(i, 0));
    if(failed) break;

    //check for convergence
    TMatrix ddelta = nrgs_new - nrgs_cur;
    converged=true;
    for(int i=0; i<3; i++) converged = converged && (fabs(ddelta(i, 0))/nrgs_cur(i, 0) <= precision_);

    //prepare for new iteration
    nrgs_cur = nrgs_new;

    if(converged) break;
  }
  
  TTKinFitter::KFResult retval;
  retval.niter = nit;
  if(converged) retval.status = CONVERGED;
  else if(failed) {
    retval.status = FAILED;
    return retval;
  }
  else retval.status = ITERLIMIT;

  //compute chi2
  TMatrix delta = nrgs_cur - nrgs_in;
  TMatrix deltaT(TMatrix::kTransposed, delta);
  retval.chi2 = (deltaT*cov_inv*delta)(0,0);

  //update permutation
  perm.WJa( )->update_energy(nrgs_cur(0,0)); 
  perm.WJb( )->update_energy(nrgs_cur(1,0)); 
  perm.BHad()->update_energy(nrgs_cur(2,0));

  return retval;
}
