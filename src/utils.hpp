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

// STL
#include <vector>
#include <string>
#include <exception>

// ROOT
#include "TH2.h"
#include "TF1.h"
#include "TString.h"

#ifndef _GERDA_FITTER_UTILS
#define _GERDA_FITTER_UTILS

namespace utils {

    std::string SafeROOTName(const std::string orig) {
        TString torig(orig);

        for (auto& c : {'.', '-', '/', ':', '|', '+'}) torig.ReplaceAll(c, '_');
        for (auto& c : {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}) {
            if (torig[0] == c) torig[0] = 'N';
        }

        return std::string(torig.Data());
    }

    TF1 ParseTFormula(std::string prefix, std::string expr, double rangelow, double rangeup) {
        // if there's a list of parameters after
        if (expr.find(':') != std::string::npos) {
            auto formula = expr.substr(0, expr.find_first_of(':'));
            auto parlist = expr.substr(expr.find_first_of(':')+1, std::string::npos);
            std::vector<double> parlist_vec;
            TF1 _tformula(
                (prefix + "_prior_tf").c_str(), formula.c_str(),
                rangelow, rangeup
            );

            if (!_tformula.IsValid()) throw std::runtime_error("ParseTFormula(): invalid prior TFormula given");

            // eventually set parameters, if any
            while (!parlist.empty()) {
                auto val = std::stod(parlist.substr(0, parlist.find_first_of(',')));
                parlist_vec.push_back(val);

                if (parlist.find(',') == std::string::npos) parlist = "";
                else parlist.erase(0, parlist.find_first_of(',')+1);
            }
            if ((int)parlist_vec.size() != _tformula.GetNpar()) {
                throw std::runtime_error("ParseTFormula(): number of values specified in '"
                        + prefix + "' TFormula does not match number of TFormula parameters");
            }

            for (size_t j = 0; j < parlist_vec.size(); ++j) _tformula.SetParameter(j, parlist_vec[j]);

            return _tformula;
        }
        // it's just the tformula
        else {
            TF1 _tformula(
                (prefix + "_prior_tf").c_str(), expr.c_str(),
                rangelow, rangeup
            );
            if (_tformula.GetNpar() > 0) {
                throw std::runtime_error("ParseTFormula(): TFormula specified with prefix '"
                        + prefix + "' is parametric but no parameters were specified");
            }
            return _tformula;
        }
    }

    template<typename BasicJsonType>
    std::vector<std::pair<double,double>> CheckAndStoreRanges(BasicJsonType& range) {
        bool _bad_range = false;
        std::vector<std::pair<double,double>> _v_range;
        if (range.is_array()) {
            if (range[0].is_array()) {
                for (auto _r : range) {
                    if (_r[0].is_number() and _r[1].is_number()) {
                        _v_range.push_back({_r[0],_r[1]});
                    }
                    else _bad_range = true;
                }
            }
            else if (range[0].is_number() and range[1].is_number()) {
                _v_range.push_back({range[0],range[1]});
            }
            else _bad_range = true;
        }
        else _bad_range = true;

        if (_bad_range) throw std::runtime_error("CheckAndStoreRanges(): Range is ill-defined");

        return _v_range;
    }

    /*
     * convert strings of the type "0:1:10,23,34,56:4:68" into vector of change points
     */
    std::vector<double> ParseBinChangePoints(std::string input) {
        if (input.empty()) throw std::runtime_error("ParseBinChangePoints(): empty input");

        std::vector<double> change_points;

        unsigned long pos;
        do {
            pos = input.find_first_of(',');
            std::string el = input.substr(0, pos);
            if (el.find(':') != std::string::npos) {
                // is of the type 'start:step:stop'
                auto elpos = el.find_first_of(':');
                std::string sstart = el.substr(0, elpos);
                el.erase(0, elpos+1);
                elpos = el.find_first_of(':');
                if (elpos == std::string::npos) {
                    throw std::invalid_argument("ParseBinChangePoints(): invalid range element format '" + el + "'");
                }
                std::string sstep = el.substr(0, elpos);
                el.erase(0, elpos+1);
                std::string sstop = el;

                double start, step, stop;
                // is it a valid number?
                try {
                    start = std::stod(sstart);
                    step = std::stod(sstep);
                    stop = std::stod(sstop);
                }
                catch (const std::invalid_argument& e) {
                    throw std::invalid_argument("ParseBinChangePoints(): stod: conversion of range element '" + el + "' failed");
                }

                // ok, now some other sanity checks
                if (stop <= start or step <= 0 or step > (stop - start)) {
                    throw std::invalid_argument("ParseBinChangePoints(): range element '" + el + "' does not make sense");
                }

                // finally we can push the change points
                for (double i = start; i < stop; i += step) change_points.push_back(i);
                change_points.push_back(stop);
            }
            else {
                // is it a valid number?
                try {
                    change_points.push_back(std::stod(el));
                }
                catch (const std::invalid_argument& e) {
                    throw std::invalid_argument("ParseBinChangePoints(): stod: conversion of range element '" + el + "' failed");
                }
            }
            input.erase(0, pos+1);
        }
        while (pos != std::string::npos);

        // sort, just to be sure
        std::sort(change_points.begin(), change_points.end());
        // and remove duplicates
        auto _last = std::unique(change_points.begin(), change_points.end());
        change_points.erase(_last, change_points.end());

        return change_points;
    }

    std::vector<std::pair<int,int>> TranslateAxisRangeToBinRange(
        TH1* /*h*/,
        std::vector<std::pair<double,double>> /*x_range*/,
        std::vector<std::pair<double,double>> /*y_range*/
    ) {
        BCLog::OutSummary("TranslateAxisRangeToBinRange(): implement me!");
        std::vector<std::pair<int,int>> _b_range;
        return _b_range;
    }

    double IntegrateHistogram1D(TH1* h, std::vector<std::pair<double,double>> range) {
        double integral = 0.;
        for (auto _r : range) {
            int _b_min = h->FindBin(_r.first);
            int _b_max = h->FindBin(_r.second);
            BCLog::OutDebug(" -> IntegrateHistogram1D ["
                + std::to_string(_b_min) + "," + std::to_string(_b_max)
                + "] n-bins : " + std::to_string(_b_max-_b_min+1)
            );
            integral += h->Integral(_b_min, _b_max);
        } 
        return integral;
    }

    double IntegrateHistogram2D(
        TH2* /*h*/,
        std::vector<std::pair<double,double>> /*x_range*/,
        std::vector<std::pair<double,double>> /*y_range*/
    )
    { 
        BCLog::OutSummary("IntegrateHistogram2D(): implement me!");
        return 0;
    }

    double IntegrateHistogram(
        TH1* h,
        std::vector<std::pair<double,double>> x_range,
        std::vector<std::pair<double,double>> y_range = {}
    ) {
        if (h->GetDimension() > 2) {
            throw std::runtime_error("IntegrateHistogram(): not implemeted for TH3.");
        }
        else if (h->GetDimension() == 2) {
            return IntegrateHistogram2D(dynamic_cast<TH2*>(h), x_range, y_range);
        }
        else if (h->GetDimension() == 1) {
            return IntegrateHistogram1D(h, x_range);
        }
        else {
            throw std::runtime_error("IntegrateHistogram(): something went unexpectedly wrong");
        }
        return 0.;
    }

}

#endif
