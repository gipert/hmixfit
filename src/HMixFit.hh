// Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef _HMIXFIT_H
#define _HMIXFIT_H

// STL
#include <vector>
#include <map>
#include <string>

// ROOT
#include "TH1.h"
#include "TFormula.h"

// BAT
#include "BAT/BCModel.h"

#include "json.hpp"
using json = nlohmann::json;

NLOHMANN_JSON_SERIALIZE_ENUM(BCEngineMCMC::Precision, {
    {BCEngineMCMC::kQuick,    "kQuick"},
    {BCEngineMCMC::kLow,      "kLow"},
    {BCEngineMCMC::kMedium,   "kMedium"},
    {BCEngineMCMC::kHigh,     "kHigh"},
    {BCEngineMCMC::kVeryHigh, "kVeryHigh"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCLog::LogLevel, {
    {BCLog::debug,   "debug"},
    {BCLog::detail,  "detail"},
    {BCLog::summary, "summary"},
    {BCLog::warning, "warning"},
    {BCLog::error,   "error"},
    {BCLog::nothing, "nothing"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCIntegrationMethod, {
    {BCIntegrate::kIntMonteCarlo, "kIntMonteCarlo"},
    {BCIntegrate::kIntGrid,       "kIntGrid"},
    {BCIntegrate::kIntLaplace,    "kIntLaplace"},
    {BCIntegrate::kIntCuba,       "kIntCuba"},
    {BCIntegrate::kIntDefault,    "kIntDefault"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCCubaMethod, {
    {BCIntegrate::kCubaDivonne, "kCubaDivonne"},
    {BCIntegrate::kCubaVegas,   "kCubaVegas"},
    {BCIntegrate::kCubaSuave,   "kCubaSuave"},
    {BCIntegrate::kCubaCuhre,   "kCubaCuhre"},
    {BCIntegrate::kCubaDefault, "kCubaDefault"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCOptimizationMethod, {
    {BCIntegrate::kOptEmpty,      "kOptEmpty"},
    {BCIntegrate::kOptSimAnn,     "kOptSimAnn"},
    {BCIntegrate::kOptMetropolis, "kOptMetropolis"},
    {BCIntegrate::kOptMinuit,     "kOptMinuit"},
    {BCIntegrate::kOptDefault,    "kOptDefault"},
})

struct dataset {
    TH1* data;                               // histogram holding data
    TH1* data_orig;                          // original input data histogram (no rebin or other stuff)
    std::vector<std::pair<int,int>> brange;  // histogram range (bin index!)
    std::map<int, TH1*> comp;                // catalog of fit components
    std::map<int, TH1*> comp_orig;           // catalog of non-rebinned fit components
    std::map<int,int> number_simulated;      // catolog of number of simulated MC events
    double livetime;

};

class HMixFit : public BCModel {

    public:

    // delete dangerous constructors
    HMixFit           ()                   = delete;
    HMixFit           (HMixFit const&) = delete;
    HMixFit& operator=(HMixFit const&) = delete;

    // custom constructor
    HMixFit(json metadata);
    ~HMixFit();

    // methods from BCModel to be overloaded
    double LogLikelihood(const std::vector<double>& parameters);
    void CalculateObservables(const std::vector<double> & parameters);

    void SetIntegrationProperties(json config);
    void PrintOptimizationSummary();
    void PrintShortMarginalizationSummary();
    void SaveHistogramsROOT(std::string filename);
    void SaveHistogramsCSV(std::string folder);
    void WriteResultsTree(std::string filename);
    double GetFastPValue(const std::vector<double>& parameters, long niter);
    double Integrate(bool enable_offset);

    std::vector<dataset> data;
    json config;
    
    private:

    std::map<std::string,TFormula> obs_tformulas;
    double _likelihood_offset = 0.; // for easier integration

    void DumpData();
    TH1* GetFitComponent(std::string filename, std::string objectname, int &nprim,TH1* tf1_hist_format = nullptr, double weight = 1);
};

#endif
