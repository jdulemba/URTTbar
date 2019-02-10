#ifndef QCD_WEIGHTPRODUCER_h
#define QCD_WEIGHTPRODUCER_h

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

class QCD_WeightProducer {

    private:

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

        double get_par( std::shared_ptr<TH1> h, string variation ) const;
        
        double SIGMOID( double p0, double p1, double p2, double X ) const;
        
            // create hists for qcd fit parameter values
        std::shared_ptr<TH1D> frac_b_CSV_p0_;

        std::shared_ptr<TH1D> frac_b_nonCSV_p0_;

        std::shared_ptr<TH1D> eff_b_p0_;
        std::shared_ptr<TH1D> eff_b_p1_;
        std::shared_ptr<TH1D> eff_b_p2_;

        std::shared_ptr<TH1D> eff_p_p0_;
        std::shared_ptr<TH1D> eff_p_p1_;
        std::shared_ptr<TH1D> eff_p_p2_;

        std::shared_ptr<TH1D> r_data_nonCSV_p0_;
        std::shared_ptr<TH1D> r_data_nonCSV_p1_;
        std::shared_ptr<TH1D> r_data_nonCSV_p2_;


    public:

        QCD_WeightProducer(); //configuration
        ~QCD_WeightProducer();

        double qcd_weight(double lep_kvar, string fit_variable, string fit_variation);

};
#endif
