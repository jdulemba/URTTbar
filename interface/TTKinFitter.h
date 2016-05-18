#include "Analyses/URTTbar/interface/Permutation.h"
#include "Analyses/URTTbar/interface/JetScaler.h"
#include <cmath>

class TTKinFitter {
public:
  enum Status {CONVERGED, ITERLIMIT, FAILED};
  struct KFResult {
    Status status;
    float chi2 = -1;
    int niter;    
  };

  TTKinFitter();
  KFResult fit(Permutation &perm);  

private:
  inline TMatrix mass_contraints_derivative(double e1, double e2, double eb, double f12, double f1b, double f2b);
  inline TMatrix mass_constraint(double e1, double e2, double eb, double f12, double f1b, double f2b);

  const float W_mass_ = 80.385;
  const float top_mass_ = 173;
  const float W_mass2_   = std::pow(W_mass_, 2);
  const float top_mass2_ = std::pow(top_mass_, 2);
  int maxiter_;
  float precision_; //convergence condition
  JetScaler &jet_scaler_;  
};


