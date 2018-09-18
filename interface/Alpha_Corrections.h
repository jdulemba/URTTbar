#ifndef ALPHA_CORRECTIONS_h
#define ALPHA_CORRECTIONS_h

#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include <string>
#include <fstream>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
#include "URAnalysis/AnalysisFW/interface/json.hpp"

using json = nlohmann::json;

class Permutation;

class Alpha_Corrections {

    private:

        double alpha_;
        double median_alpha_;
        json jsonFile_;

    public:

        Alpha_Corrections(); //configuration
        ~Alpha_Corrections();

        double alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range);
        TLorentzVector Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range);

        double median_alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range);
        TLorentzVector median_Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range);

        double alpha_thad_cthstar(Permutation &perm, string fit_degree, string fit_var, string fit_range);

};
#endif
