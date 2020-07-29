// MIT License
//
// Copyright (c) 2020 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#ifndef _GERDA_FITTER_H
#define _GERDA_FITTER_H

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
    std::vector<std::pair<int,int>> brange;  // histogram range (bin index!)
    std::map<int, TH1*> comp;                // catalog of fit components
};

class GerdaFitter : public BCModel {

    public:

    // delete dangerous constructors
    GerdaFitter           ()                   = delete;
    GerdaFitter           (GerdaFitter const&) = delete;
    GerdaFitter& operator=(GerdaFitter const&) = delete;

    // custom constructor
    GerdaFitter(json metadata);
    ~GerdaFitter();

    // methods from BCModel to be overloaded
    double LogLikelihood(const std::vector<double>& parameters);
    void CalculateObservables(const std::vector<double> & parameters);

    void SetIntegrationProperties(json config);
    void PrintOptimizationSummary();
    void PrintShortMarginalizationSummary();
    void SaveHistograms(std::string filename);
    void WriteResultsTree(std::string filename);
    double GetFastPValue(const std::vector<double>& parameters, long niter);
    double Integrate(bool enable_offset);

    std::vector<dataset> data;
    json config;

    private:

    std::map<std::string,TFormula> obs_tformulas;
    double _likelihood_offset = 0.; // for easier integration

    void DumpData();
    TH1* GetFitComponent(std::string filename, std::string objectname, TH1* data, int rebin_x = 1, int rebin_y = 1, std::vector<double> change_points = {});
};

#endif
