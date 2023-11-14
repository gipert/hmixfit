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

// STL
#include <vector>
#include <string>
#include <exception>

// BAT
#include "BAT/BCLog.h"

// ROOT
#include "TH2.h"
#include "TF1.h"
#include "TString.h"

#ifndef _HMIXFIT_UTILS
#define _HMIXFIT_UTILS

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
        throw std::runtime_error("TranslateAxisRangeToBinRange(): implement me!");
        std::vector<std::pair<int,int>> _b_range;
        return _b_range;
    }

    double IntegrateHistogram1D(TH1* h, std::vector<std::pair<double,double>> range) {
        double integral = 0.;
        for (auto _r : range) {
            int _b_min = h->FindBin(_r.first);
            int _b_max = h->FindBin(_r.second);
            BCLog::OutDebug("IntegrateHistogram1D(): about to integrate "
                + std::string(h->GetName()) + " in ["
                + std::to_string(_b_min) + ", " + std::to_string(_b_max)
                + "] (" + std::to_string(_b_max-_b_min+1) + " bins)"
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
        throw std::runtime_error("IntegrateHistogram2D(): implement me!");
        return 0;
    }

    double IntegrateHistogram(
        TH1* h,
        std::vector<std::pair<double,double>> x_range,
        std::vector<std::pair<double,double>> y_range = {}
    ) {
        if (!h) throw std::invalid_argument("IntegrateHistogram(): input histogram is null");

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

    TH1* ReformatHistogram(TH1* original, TH1* htemplate) {

        if (!original or !htemplate) {
            throw std::invalid_argument("ReformatHistogram(): invalid (nullptr) input");
        }

        if (htemplate->GetXaxis()->IsVariableBinSize()) {
            BCLog::OutDebug("ReformatHistogram(): template histogram is variably binned, no consistency checks will be performed! (FIXME)");
        }
        else {
            auto nx = original->GetNbinsX();
            auto ny = original->GetNbinsY();
            auto nz = original->GetNbinsZ();
            auto ax = original->GetXaxis()->GetXmin();
            auto bx = original->GetXaxis()->GetXmax();
            auto ay = original->GetYaxis()->GetXmin();
            auto by = original->GetYaxis()->GetXmax();
            auto az = original->GetZaxis()->GetXmin();
            auto bz = original->GetZaxis()->GetXmax();

            auto nx_ = htemplate->GetNbinsX();
            auto ny_ = htemplate->GetNbinsY();
            auto nz_ = htemplate->GetNbinsZ();
            auto ax_ = htemplate->GetXaxis()->GetXmin();
            auto bx_ = htemplate->GetXaxis()->GetXmax();
            auto ay_ = htemplate->GetYaxis()->GetXmin();
            auto by_ = htemplate->GetYaxis()->GetXmax();
            auto az_ = htemplate->GetZaxis()->GetXmin();
            auto bz_ = htemplate->GetZaxis()->GetXmax();

            BCLog::OutDebug(Form("ReformatHistogram(): %s (original) = h(%f:%d:%f, %f:%d:%f, %f:%d:%f)",
                original->GetName(), ax, nx, bx, ay, ny, by, az, nz, bz));

            BCLog::OutDebug(Form("ReformatHistogram(): %s (htemplate) = h(%f:%d:%f, %f:%d:%f, %f:%d:%f)",
                htemplate->GetName(), ax_, nx_, bx_, ay_, ny_, by_, az_, nz_, bz_));

            if (nx == nx_ and ny == ny_ and nz == nz_ and
                ax == ax_ and ay == ay_ and az == az_ and
                bx == bx_ and by == by_ and bz == bz_) {

                BCLog::OutDebug("ReformatHistogram(): original and template histogram seem already consistent, exiting.");
                auto th_new = dynamic_cast<TH1*>(original->Clone());
                th_new->SetDirectory(nullptr);
                return th_new;
            }
            else {
                BCLog::OutDebug("ReformatHistogram(): template histogram looks different from original, proceeding. No consistency checks will be performed! (FIXME)");
                // FIXME: check whether the original histogram can be legitimately cast into new one
                // if (fmod(ax_, (bx -ax)/nx) != 0) {
                // }
            }
        }

        // clone htemplate
        auto th_new = dynamic_cast<TH1*>(htemplate->Clone());
        th_new->Reset();
        // detach it from original
        th_new->SetDirectory(nullptr);
        th_new->SetName(original->GetName());
        th_new->SetTitle(original->GetTitle());

        // fill new histogram
        for (int b = 0; b <= original->GetNcells(); ++b) {
            // first we get the bin index for each axis
            int bx = 0, by = 0, bz = 0;
            original->GetBinXYZ(b, bx, by, bz);
            // wee need the coordinates of the old bin
            auto x = original->GetXaxis()->GetBinCenter(bx);
            auto y = original->GetYaxis()->GetBinCenter(by);
            auto z = original->GetZaxis()->GetBinCenter(bz);

            // now find the corresponding bin in the new histogram and fill it
            auto newbin = th_new->FindBin(x, y, z);
            th_new->SetBinContent(newbin, th_new->GetBinContent(newbin) + original->GetBinContent(b));
        }

        return th_new;
    }

    /*
     * smallest_poisson_interval(prob_coverage, poisson_mean)
     *
     * Computes the smallest interval boundaries covering `prob_coverage` area of a
     * discrete Poisson distribution with mean `poisson_mean` Returns a
     * `std::pair<float,float>` holding the lower and upper range
     */
    std::pair<float, float> smallest_poisson_interval(double cov, double mu) {

        // sanity check
        if (cov > 1 or cov < 0 or mu < 0) throw std::runtime_error("smallest_poisson_interval(): bad input");

        // initialize lower and upper edges to something
        std::pair<float, float> res = {mu, mu};

        if (mu > 50) { // gaussian approximation os OK
            res = {
                std::round(mu + TMath::NormQuantile((1-cov)/2)*sqrt(mu))-0.5,
                std::round(mu - TMath::NormQuantile((1-cov)/2)*sqrt(mu))+0.5
            };
        }
        else { // do the computation
            // start from the mode, which is the integer part of the mean
            int mode = std::floor(mu);
            int l = mode, u = mode; // let's start from here
            double prob = TMath::PoissonI(mode, mu); // probability covered by the interval

            // check if we're undercovering
            while (prob < cov) {
                // compute probabilities of points just ouside interval
                double prob_u = TMath::PoissonI(u+1, mu);
                double prob_l = TMath::PoissonI(l > 0 ? l-1 : 0, mu);

                // we expand on the right if:
                //  - the lower edge is already at zero
                //  - the prob of the right point is higher than the left
                if (l == 0 or prob_u > prob_l) {
                    u++; // expand interval
                    prob += prob_u; // update coverage
                }
                // otherwhise we expand on the left
                else if (prob_u < prob_l) {
                    l--;
                    prob += prob_l;
                }
                // if prob_u == prob_l we expand on both sides
                else {
                    u++; l--;
                    prob += prob_u + prob_l;
                }
            }
            res = {l == 0 ? 0 : l-0.5, u+0.5};
        }
        return res;
    }

    double normalized_poisson_residual(double mu, double obs) {

        TF1 fu("f", "TMath::Poisson(x,[0])", 0, 15);
        fu.SetNpx(1000);

        auto median = [](double x) {

            double up = x + 1./3;
            double dw = x - std::log(2.) > 0 ? x - std::log(2.) : 0;
            double mprec = 0.01;

            TF1 f_tmp("f_tmp", "TMath::Poisson(x,[0])", 0, 1);
            f_tmp.SetParameter(0, x);
            f_tmp.SetNpx(1000);
            double norm = f_tmp.Integral(0, x < 1 ? 5 : 10*x);

            double I = f_tmp.Integral(0, dw)/norm;
            while (I < 0.5) {
                dw += mprec;
                I = f_tmp.Integral(0, dw)/norm;
                if (dw > up) {
                    std::cout << "normalized_poisson_residual(): "
                              << "WARNING: failed to calculate median" << std::endl;
                    break;
                }
            }

            return dw;
        };

        if (mu < 50) {
            fu.SetParameter(0, mu);

            double prob = fu.Integral(median(mu), obs);
            if (std::abs(prob) >= 0.5) {
	      //std::cout << "normalized_poisson_residual(): "
	      //          << "WARNING: poisson_area(" << obs << ", " << median(mu)
	      //          << ") = " << prob << std::endl;
                return (prob == 1 ? 3.90 : -3.9);
            }
            else return TMath::NormQuantile(0.5+prob);
        }
        else {
            return (obs-mu)/std::sqrt(mu);
        }
    }
}

#endif
