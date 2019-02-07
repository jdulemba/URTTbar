#ifndef ALPHA_CORRECTIONS_h
#define ALPHA_CORRECTIONS_h

#include "Analyses/URTTbar/interface/Permutation.h"
#include <map>
#include <string>
#include <fstream>
#include "URAnalysis/AnalysisFW/interface/Logger.h"
//#include "URAnalysis/AnalysisFW/interface/json.hpp"
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include "URAnalysis/AnalysisFW/interface/URParser.h"
#include "URAnalysis/AnalysisFW/interface/DataFile.h"


//using json = nlohmann::json;
using namespace std;

class Permutation;

class Alpha_Corrections {

    private:

        //double alpha_;
        //double median_alpha_;
        ////json jsonFile_;

        //TFile alpha_file_;

        template <class T>
            std::shared_ptr<T> get_from(TFile &file, std::string path, std::string newname) {
                T* original = (T*) file.Get( path.c_str() );
                if(!original) {
                    Logger::log().warning() << "Could not get " << path << " from the file " << file.GetName() << "!" << std::endl;
                    std::shared_ptr<T> ptr;
                    return ptr;
                }
                std::shared_ptr<T> ptr((T*) original->Clone(newname.c_str()));
                return ptr;
            }

        double get_1D_corr( std::shared_ptr<TH1> h, double mthad) const;
        double get_2D_corr( std::shared_ptr<TH2> h, double mthad, double mtt) const;

        // create hists of alpha values
        std::shared_ptr<TH2D> THad_E_Mtt_1d_;
        std::shared_ptr<TH2D> THad_E_Mtt_2d_;
        std::shared_ptr<TH2D> THad_P_Mtt_1d_;
        std::shared_ptr<TH2D> THad_P_Mtt_2d_;
        std::shared_ptr<TH2D> THad_M_Mtt_1d_;
        std::shared_ptr<TH2D> THad_M_Mtt_2d_;
        std::shared_ptr<TH1D> THad_E_All_1d_;
        std::shared_ptr<TH1D> THad_E_All_2d_;
        std::shared_ptr<TH1D> THad_P_All_1d_;
        std::shared_ptr<TH1D> THad_P_All_2d_;
        std::shared_ptr<TH1D> THad_M_All_1d_;
        std::shared_ptr<TH1D> THad_M_All_2d_;


    public:

        Alpha_Corrections(); //configuration
        ~Alpha_Corrections();

        double alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range);
        TLorentzVector Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range);

        //double median_alpha(Permutation &perm, string fit_degree, string fit_var, string fit_range);
        //TLorentzVector median_Alpha_THad(Permutation &perm, string fit_degree, string fit_var, string fit_range);

        double alpha_thad_cthstar(Permutation &perm, string fit_degree, string fit_var, string fit_range);

};
#endif
