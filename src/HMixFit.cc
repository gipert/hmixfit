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

#include <cassert>

#include "HMixFit.hh"
#include "utils.hpp"

// BAT
#include "BAT/BCMath.h"
#include "BAT/BCAux.h"
#include "BAT/BCTF1Prior.h"
#include "BAT/BCTH1Prior.h"

// ROOT
#include "TF1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TParameter.h"
#include "TObjString.h"
#include <limits>

HMixFit::HMixFit(json outconfig) : config(outconfig) {

    TH1::AddDirectory(false);

    // open log file

    auto outdir = config["output-dir"].get<std::string>();
    auto data_path =config.value("data-path","");
    auto mc_path =config.value("pdf-path","");
    auto hist_name =config.value("hist-dir","raw");
    auto use_priors = config.value("use-priors",true);
    auto livetime =config["livetime"].get<double>();
    
    std::cout<<"Livetime = "<<livetime<<std::endl;
    fLivetime = livetime;
    
    std::system(("mkdir -p " + outdir).c_str());
    auto prefix = outdir + "/hmixfit-" + config["id"].get<std::string>() + "-";

    //    int n_bins = std::stoi(config["n-bins"].get<std::string>());
    int n_bins = config.value("n-bins", 1000);
    
    //    BCLog::SetLogLevelFile(BCLog::debug);
    BCLog::OpenLog(prefix + "output.log");
    BCLog::OutSummary("Saving results in " + outdir);
    //BCLog::SetLogLevelScreen(config.value("logging", BCLog::summary));

    this->SetName(config["id"]);

    /*
     * create fit parameters
     */

    // for (auto& el : config["fit"]["parameters"].items()) {
    for (json::iterator el = config["fit"]["parameters"].begin(); el != config["fit"]["parameters"].end(); ++el) {
        // do we already have a parameter with the same name?
        bool already_exists = false;
        for (unsigned int idx = 0; idx < this->GetNParameters(); ++idx) {
            if (this->GetParameters().At(idx).GetName() == el.key()) {
                already_exists = true;
                break;
            }
        }
        if (already_exists) {
            BCLog::OutDetail("model parameter '" + el.key() + "' already exists, skipping");
        }
        else {
            if (el.value().contains("range")) {
                this->GetParameters().Add(
                    el.key(),
                    el.value()["range"][0].get<double>(),
                    el.value()["range"][1].get<double>(),
                    el.value().value("long-name", ""),
                    "(" + el.value().value("units", "") + ")"
                );
                BCLog::OutDetail("adding model parameter '" + el.key() + "' (\""
                    + el.value().value("long-name", "") + "\" [" + el.value().value("units", "")
                    + "]) in range = [" + std::to_string(el.value()["range"][0].get<double>()) + ","
                    + std::to_string(el.value()["range"][1].get<double>()) + "]");
            }
            else if (el.value().contains("fixed")) {
                // fix it?
                if (!el.value()["fixed"].is_number()) {
                    throw std::runtime_error("\"fixed\" value of parameter " + el.key() + " is not a number");
                }
                this->GetParameters().Add(
                    el.key(),
                    el.value()["fixed"].get<double>()-1,
                    el.value()["fixed"].get<double>()+1,
                    el.value().value("long-name", ""),
                    "(" + el.value().value("units", "") + ")"
                );
                BCLog::OutDetail("fixing parameter " + el.key() + " to "
                        + std::to_string(el.value()["fixed"].get<double>()));
                this->GetParameters().Back().Fix(el.value()["fixed"].get<double>());
            }
        }

	// binning of parameters
	this->GetParameters().SetNBins(n_bins);
	
        // assign prior

	if (use_priors==true &&el.value().contains("prior")) {
	  auto& prior_cfg = el.value()["prior"];
	  TH1* prior_hist = nullptr;
	  BCPrior* prior = nullptr;
	  std::string expr;

	  if (prior_cfg.contains("TFormula") and prior_cfg.contains("histogram")) {
	    throw std::runtime_error("please choose either \"TFormula\" or \"histogram\" for parameter " + el.key());
	  }
	  // a ROOT histogram is given
	  if (prior_cfg.contains("histogram")) {
	    expr = prior_cfg["histogram"].get<std::string>();
	    if (expr.find(':') == std::string::npos) {
		  throw std::runtime_error("invalid \"histogram\" format for parameter " + el.key());
	    }
	    auto filename = expr.substr(0, expr.find_first_of(':'));
	    auto objname = expr.substr(expr.find_first_of(':')+1, std::string::npos);
	    TFile _tff(filename.c_str());
	    if (!_tff.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
	    auto obj = _tff.Get(objname.c_str());
	    if (!obj) throw std::runtime_error("could not find object '" + objname + "' in file " + filename);
	    if (obj->InheritsFrom(TH1::Class())) {
	      prior_hist = dynamic_cast<TH1*>(obj);
	      // the following is a workaround for BAT's bad implementation of the BCTH1Prior constructor
	      prior_hist->SetDirectory(nullptr);
	      _tff.Close();
                }
	    else throw std::runtime_error("object '" + objname + "' in file " + filename + " is not an histogram");
	    BCLog::OutDebug("found histogram prior '" + expr + "' for parameter '" + el.key() + "'");
            }
	  // is a tformula
	  else if (prior_cfg.contains("TFormula")) {
	    expr = prior_cfg["TFormula"].get<std::string>();
	    TF1 _tformula = utils::ParseTFormula(
						 el.key(),
						 expr,
						 el.value()["range"][0].get<double>(),
						 el.value()["range"][1].get<double>()
						 );
                BCLog::OutDebug("found TFormula prior '" + expr + "' for parameter '" + el.key() + "'");

                // I did not manage to use BCTF1Prior from BAT, it segfaults for misterious reasons
                // so I just make a temporary histogram out of the TFormula
                prior_hist = new TH1D((el.key() + "_prior_h").c_str(), "h1_prior", 1000,
                    el.value()["range"][0].get<double>(),
                    el.value()["range"][1].get<double>()
                );
                for (int j = 1; j <= prior_hist->GetNbinsX(); j++) {
                    prior_hist->SetBinContent(j, _tformula.Eval(prior_hist->GetBinCenter(j)));
                }
                BCLog::OutDebug("produced temporary prior 1000-bin histogram for parameter '" + el.key() + "'");
            }
            else throw std::runtime_error("supported prior types: 'histogram', 'TFormula'");

            prior = new BCTH1Prior(prior_hist);
            delete prior_hist;
            BCLog::OutDetail("assigned prior '" + expr + "' to parameter '" + el.key() + "'");

            if (prior != nullptr) {
                if (!prior->IsValid()) {
                    throw std::runtime_error("invalid prior set for parameter " + el.key());
                }
                this->GetParameters().Back().SetPrior(prior);
            }
            else {
                this->GetParameters().Back().SetPriorConstant();
                BCLog::OutDetail("assigned uniform prior to parameter '" + el.key() + "'");
            }
        }
        else {
            this->GetParameters().Back().SetPriorConstant();
            BCLog::OutDetail("assigned uniform prior to parameter '" + el.key() + "'");
        }
    }

    /*
     * read in histograms
     */

    // loop over files with data histograms
    for (auto& el : config["fit"]["theoretical-expectations"].items()) {
        BCLog::OutDebug("opening data file " + el.key());
	std::string file_name = data_path+el.key();
        TFile _tf(file_name.c_str());
        if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + el.key());

        // loop over requested histograms in file
        for (auto& elh : el.value().items()) {
            // get rebin factor before reading all histograms
            std::vector<double> change_points;
            int rebin_x = 1, rebin_y = 1;
            if (elh.value().contains("rebin-factor") and elh.value()["rebin-factor"].is_string()) {
                BCLog::OutDetail("using specified variable-sized rebin for dataset '" + el.key() + "'");
                change_points = utils::ParseBinChangePoints(elh.value()["rebin-factor"].get<std::string>());
            }
            else {
                rebin_x = elh.value().value("rebin-factor-x", 1) * elh.value().value("rebin-factor", 1);
                rebin_y = elh.value().value("rebin-factor-y", 1);
                if (rebin_x != 1)
                    BCLog::OutDetail("using rebin X factor = " + std::to_string(rebin_x) + " for dataset '" + elh.key() + "'");
                if (rebin_y != 1)
                    BCLog::OutDetail("using rebin Y factor = " + std::to_string(rebin_y) + " for dataset '" + elh.key() + "'");
            }

            BCLog::OutDebug("getting data histogram '" + elh.key() + "' from " + el.key());
            auto th_orig = dynamic_cast<TH1*>(_tf.Get(elh.key().c_str()));
            if (!th_orig) throw std::runtime_error("could not find object '" + elh.key() + "' in file " + el.key());
            if (th_orig->GetDimension() != 1 and !change_points.empty()) {
                throw std::runtime_error("cannot apply variable-sized rebin to TH2 or TH3");
            }

            // set explicative histogram names
            auto basename = el.key().substr(
                    el.key().find_last_of('/')+1,
                    el.key().find_last_of('.') - el.key().find_last_of('/') - 1);
            basename += "_" + std::string(th_orig->GetName());

            th_orig->SetName((utils::SafeROOTName(basename) + "_orig").c_str());

            dataset _current_ds;
            _current_ds.data_orig = th_orig;

            // rebin if requested
            TH1* th = nullptr;
            if (change_points.empty()) {
                th = dynamic_cast<TH1*>(th_orig->Clone(utils::SafeROOTName(basename).c_str()));
                th->Rebin(rebin_x);
                if ((rebin_y != 1) and th->InheritsFrom(TH2::Class())) {
                    dynamic_cast<TH2*>(th)->RebinY(rebin_y);
                }
            }
            else {
                // do variable rebin
                th = dynamic_cast<TH1*>(
                    th_orig->Rebin(
                        change_points.size()-1,
                        utils::SafeROOTName(basename).c_str(),
                        &change_points[0]
                    ));
            }

            _current_ds.data = th;

            // eventually get a global value for the gerda-pdfs path
            auto gerda_pdfs_path = elh.value().value("gerda-pdfs", ".");

            BCLog::OutDebug("getting the requested PDFs");
            // loop over requested components
            for (auto& it : elh.value()["components"]) {

                prefix = it.value("prefix", gerda_pdfs_path);

                /* START INTERMEZZO */
                // utility to sum over the requested parts (with weight) given isotope
                auto sum_parts = [&](std::string i) {
                    std::string true_iso = i;
                    if (i.find('-') != std::string::npos) true_iso = i.substr(0, i.find('-'));

                    TH1* sum = nullptr;
                    double sum_prim=0;
                    if (it["part"].is_object()) {
                        for (auto& p : it["part"].items()) {
                            // get volume name
                            auto path_to_part = prefix + "/" + p.key();
                            if (path_to_part.back() == '/') path_to_part.pop_back();
                            auto part = path_to_part.substr(path_to_part.find_last_of('/')+1);
                            path_to_part.erase(path_to_part.find_last_of('/'));
                            auto volume = path_to_part.substr(path_to_part.find_last_of('/')+1);
                            auto filename = prefix + "/" + p.key() + "/" + true_iso + "/" + "pdf-"
                                + volume + "-" + part + "-" + i + ".root";
                            BCLog::OutDebug("opening file " + filename);
                            BCLog::OutDebug("summing object '" + elh.key() + " with weight "
                                    + std::to_string(p.value().get<double>()));
                            // get histogram
                            int nprim_tmp;
                            auto thh = this->GetFitComponent(filename, elh.key(),nprim_tmp, _current_ds.data_orig);
                            // add it with weight
                            if (!sum) {
                                sum = thh;
                                sum->Scale(p.value().get<double>());
                                sum_prim=nprim_tmp*p.value().get<double>();
                            }
                            else {
                                sum->Add(thh, p.value().get<double>());
                                sum_prim+=nprim_tmp*p.value().get<double>();
                                delete thh;
                            }
                        }
                        return sum;
                    }
                    else if (it["part"].is_string()) {
                        // get volume name
                        auto path_to_part = prefix + "/" + it["part"].get<std::string>();
                        if (path_to_part.back() == '/') path_to_part.pop_back();
                        auto part = path_to_part.substr(path_to_part.find_last_of('/')+1);
                        path_to_part.erase(path_to_part.find_last_of('/'));
                        auto volume = path_to_part.substr(path_to_part.find_last_of('/')+1);
                        auto filename = prefix + "/" + it["part"].get<std::string>() + "/" + true_iso + "/" + "pdf-"
                            + volume + "-" + part + "-" + i + ".root";
                        BCLog::OutDebug("getting object '" + elh.key() + "' in file " + filename);
                        // get histogram
                        int nprim=0;
                        auto thh = this->GetFitComponent(filename, elh.key(),nprim, _current_ds.data_orig);

                        return thh;
                    }
                    else throw std::runtime_error("unexpected 'part' value found in [\"fit\"][\""
                            + el.key() + "\"][\"" + elh.key() + "\"]");
                };
                /* END INTERMEZZO */

                // loop over requested isotopes on the relative part
                for (auto& iso : it["components"].items()) {
                    BCLog::OutDetail("adding a PDF for the '" + elh.key() + "' dataset with '"
                        + iso.key() + "' as scaling parameter (debug mode for details)");
                    // throw exception if component name is not in list
                    bool exists = false;
                    int comp_idx = 0;
                    for (unsigned int idx = 0; idx < this->GetNParameters(); ++idx) {
                        if (this->GetParameters().At(idx).GetName() == iso.key()) {
                            comp_idx = idx;
                            exists = true;
                            break;
                        }
                    }
                    if (!exists) throw std::runtime_error(
                        "fit parameter '" + iso.key() + "' not found, is it defined in \"parameters\"?"
                    );

                    // it's a user defined file
                    if (!iso.value().contains("TFormula") and it.contains("root-file")) {
                        BCLog::OutDebug("user-specified ROOT file detected");
                        auto obj_name = iso.value().value("hist-name", elh.key());
                        int nprim;
                        auto thh = this->GetFitComponent(mc_path+it["root-file"].get<std::string>(), obj_name,nprim, _current_ds.data_orig);
                        thh->SetName(utils::SafeROOTName(iso.key()).c_str());
                        _current_ds.comp_orig.insert({comp_idx, thh});
                        _current_ds.number_simulated.insert({comp_idx,nprim});
                        _current_ds.comp.insert({comp_idx, utils::ReformatHistogram(thh, _current_ds.data)});
                    }
                    // it's a user defined file list with weights
                    else if (!iso.value().contains("TFormula") and it.contains("root-files")) {
                        BCLog::OutDebug("user-specified ROOT file list (with weights) detected");
                        auto obj_name = iso.value().value("hist-name", elh.key());

                        TH1* sum = nullptr;
                      
                        if (it["root-files"].is_object()) {
                            for (auto& p : it["root-files"].items()) {
                                // get histogram
                                BCLog::OutDebug("opening file " + p.key());
                                BCLog::OutDebug("summing object '" + obj_name + " with weight "
                                        + std::to_string(p.value().get<double>()));
                                int nprim_tmp;
                                auto thh = this->GetFitComponent(p.key(), obj_name,nprim_tmp, _current_ds.data_orig);
                                // add it with weight
                                if (!sum) {
                                    sum = thh;
                                    sum->Scale(p.value().get<double>());
                                }
                                else {
                                    sum->Add(thh, p.value().get<double>());
                                    delete thh;
                                }
                            }
                        }
                        else throw std::runtime_error("\"root-files\" must be a dictionary");
                        sum->SetName(utils::SafeROOTName(iso.key()).c_str());
                        _current_ds.comp_orig.insert({comp_idx, sum});
                        _current_ds.number_simulated.insert({comp_idx,std::numeric_limits<double>::max()});
                        _current_ds.comp.insert({comp_idx, utils::ReformatHistogram(sum, _current_ds.data)});
                    }
                    // it's a explicit TFormula
                    else if (iso.value().contains("TFormula")) {
                        BCLog::OutDebug("TFormula expression detected");
                        if (!iso.value()["TFormula"].is_string()) {
                            throw std::runtime_error("The \"TFormula\" key must be a string");
                        }
                        auto _tfunc = utils::ParseTFormula(
                            iso.key(),
                            iso.value()["TFormula"].get<std::string>(),
                            _current_ds.data->GetXaxis()->GetXmin(),
                            _current_ds.data->GetXaxis()->GetXmax()
                        );
                        if (_tfunc.GetNdim() != _current_ds.data->GetDimension()) {
                            throw std::runtime_error("the specified PDF TFormula has a different dimensionality than the data");
                        }
                        auto thh = new TH1D(
                            iso.key().c_str(), iso.key().c_str(),
                            _current_ds.data_orig->GetNbinsX(),
                            _current_ds.data_orig->GetXaxis()->GetXmin(),
                            _current_ds.data_orig->GetXaxis()->GetXmax()
                        );
                        for (int b = 1; b < thh->GetNbinsX(); ++b) {
                            thh->SetBinContent(b, _tfunc.Eval(thh->GetBinCenter(b)));
                        }
                        thh->SetName(utils::SafeROOTName(iso.key()).c_str());
                        _current_ds.comp_orig.insert({comp_idx, thh});
                        _current_ds.comp.insert({comp_idx, utils::ReformatHistogram(thh, _current_ds.data)});
                    }
                    else { // look into gerda-pdfs database
                        BCLog::OutDebug("should be a gerda-pdfs PDF");
                        if (iso.value()["isotope"].is_string()) {
                            auto comp = sum_parts(iso.value()["isotope"]);
                            comp->SetName(utils::SafeROOTName(iso.key()).c_str());
                            _current_ds.comp_orig.insert({comp_idx, comp});
                            _current_ds.comp.insert({comp_idx, utils::ReformatHistogram(comp, _current_ds.data)});
                        }
                        else if (iso.value()["isotope"].is_object()) {
                            TH1* comp = nullptr;
                            for (auto& i : iso.value()["isotope"].items()) {
                                BCLog::OutDebug("scaling pdf for " + i.key() + " by a factor "
                                        + std::to_string(i.value().get<double>()));
                                if (!comp) {
                                    comp = sum_parts(i.key());
                                    comp->Scale(i.value().get<double>());
                                }
                                else {
                                    auto _htmp = sum_parts(i.key());
                                    comp->Add(_htmp, i.value().get<double>());
                                    delete _htmp;
                                }

                            }
                            comp->SetName(utils::SafeROOTName(iso.key()).c_str());
                            _current_ds.comp_orig.insert({comp_idx, comp});
                            _current_ds.comp.insert({comp_idx, utils::ReformatHistogram(comp, _current_ds.data)});
                        }
                        else throw std::runtime_error("unexpected entry " + iso.value()["isotope"].dump() + "found in [\"fit\"][\""
                                + el.key() + "\"][\"" + elh.key() + "\"][\"components\"][\"" + iso.key() + "\"][\"isotope\"]");

                        if (!_current_ds.comp[comp_idx] or !_current_ds.comp_orig[comp_idx]) {
                            throw std::runtime_error("invalid pointer found in component list at position " + std::to_string(comp_idx));
                        }
                    }
                }
            }

            // set fit ranges
            // read x-range from config file
            std::vector<std::pair<double,double>> _vxr;
            if (!elh.value().contains("fit-range-x") and !elh.value().contains("fit-range")) {
                _vxr.push_back({
                    _current_ds.data->GetXaxis()->GetBinCenter(1),
                    _current_ds.data->GetXaxis()->GetBinCenter(_current_ds.data->GetNbinsX())
                });
            }
            else { // sanity check x-range
                auto x_range = elh.value().contains("fit-range-x") ? elh.value()["fit-range-x"] : elh.value()["fit-range"];
                if (!x_range.is_array())
                    throw std::runtime_error("wrong \"fit-range-x\" format, array expected");
                if (!(x_range[0].is_array() or x_range[0].is_number()))
                    throw std::runtime_error("wrong \"fit-range-x\" format, array of arrays or numbers expected");

                if (x_range[0].is_array()) {
                    BCLog::OutDebug("\"fit-range-x\" is an array of arrays");
                    for (auto& r : x_range) {
                        if (!r.is_array()) throw std::runtime_error("non-array member detected in \"fit-range-x\"");
                        if (r.size() != 2) throw std::runtime_error("wrong number of entries in item of \"fit-range-x\"");
                        if (!r[0].is_number() or !r[1].is_number()) {
                            throw std::runtime_error("not-a-number value detected in \"fit-range-x\"");
                        }
                        _vxr.push_back({r[0], r[1]});
                    }
                }
                else _vxr.push_back({x_range[0], x_range[1]});
            }

            // set fit range 1D
            if (_current_ds.data->GetDimension() == 1) {
                if (_vxr.empty()) throw std::runtime_error("Something went wrong with the TH1 x-range");
                else {
                    for (auto x : _vxr) {
                        // no under or overflow bins allowed
                        auto _x_ll = _current_ds.data->FindBin(x.first);
                        auto _x_ul = _current_ds.data->FindBin(x.second);
                        while (_x_ll <= 0)                            _x_ll++;
                        while (_x_ul > _current_ds.data->GetNbinsX()) _x_ul--;
                        _current_ds.brange.push_back({_x_ll,_x_ul});
                        BCLog::OutDetail("Adding fit range x [" +
                            std::to_string(_current_ds.brange.back().first) + ", " +
                            std::to_string(_current_ds.brange.back().second) + "] (bins) [" +
                            std::to_string(_current_ds.data->GetBinLowEdge(_x_ll)) + ", " +
                            std::to_string(_current_ds.data->GetBinLowEdge(_x_ul) + _current_ds.data->GetBinWidth(_x_ul))
                            + "] (x-units)"
                        );
                    }
                }
            }
            // set fit range 2D
            else if (_current_ds.data->GetDimension() == 2) {
                // read y-range from config file
                std::vector<std::pair<double,double>> _vyr;
                if (!elh.value().contains("fit-range-y")) {
                    _vyr.push_back({
                        _current_ds.data->GetYaxis()->GetBinCenter(1),
                        _current_ds.data->GetYaxis()->GetBinCenter(_current_ds.data->GetNbinsY())
                    });
                }
                else { // sanity check y-range
                    auto & y_range = elh.value()["fit-range-y"];
                    if (!y_range.is_array())
                        throw std::runtime_error("wrong \"fit-range-y\" format, array expected");
                    if (!(y_range[0].is_array() or y_range[0].is_number()))
                        throw std::runtime_error("wrong \"fit-range-y\" format, array of arrays or numbers expected");

                    if (y_range[0].is_array()) {
                        BCLog::OutDebug("\"fit-range-y\" is an array of arrays");
                        for (auto& r : y_range) {
                            if (!r.is_array()) throw std::runtime_error("non-array member detected in \"fit-range-y\"");
                            if (r.size() != 2) throw std::runtime_error("wrong number of entries in item of \"fit-range-y\"");
                            if (!r[0].is_number() or !r[1].is_number()) {
                                throw std::runtime_error("not-a-number value detected in \"fit-range-y\"");
                            }
                            _vyr.push_back({r[0],r[1]});
                        }
                    }
                    else _vyr.push_back({y_range[0],y_range[1]});
                }
                // last sanity check of ranges
                if (_vxr.empty() or _vyr.empty())
                    throw std::runtime_error("Something went wrong with the TH2 ranges");
                else {
                    // translate ranges in global bin ranges
                    auto _y_bin_width = _current_ds.data->GetYaxis()->GetBinWidth(1);
                    for (auto x : _vxr) {
                        // no under or overflow bins allowed
                        auto _x_ll = _current_ds.data->GetXaxis()->FindBin(x.first);
                        auto _x_ul = _current_ds.data->GetXaxis()->FindBin(x.second);
                        while (_x_ll <= 0)                            _x_ll++;
                        while (_x_ul > _current_ds.data->GetNbinsX()) _x_ul--;
                        for (auto y : _vyr) {
                            // no under or overflow bins allowed
                            auto _y_ll = _current_ds.data->GetYaxis()->FindBin(y.first);
                            auto _y_ul = _current_ds.data->GetYaxis()->FindBin(y.second);
                            while (_y_ll <= 0)                            _y_ll++;
                            while (_y_ul > _current_ds.data->GetNbinsY()) _y_ul--;
                            while (_y_ll <= _y_ul) {
                                _current_ds.brange.push_back({
                                    _current_ds.data->FindBin(_x_ll,_y_ll),
                                    _current_ds.data->FindBin(_x_ul,_y_ll)
                                });
                                BCLog::OutDetail("Adding fit range TH2 global [" +
                                    std::to_string(_current_ds.brange.back().first) + ", " +
                                    std::to_string(_current_ds.brange.back().second) + "] (bins) [" +
                                    std::to_string(_current_ds.data->GetXaxis()->GetBinLowEdge(_x_ll)) + ", " +
                                    std::to_string(_current_ds.data->GetXaxis()->GetBinLowEdge(_x_ul) +
                                    _current_ds.data->GetXaxis()->GetBinWidth(_x_ul)) + " : " +
                                    std::to_string(_current_ds.data->GetYaxis()->GetBinLowEdge(_y_ll)) + ", " +
                                    std::to_string(_current_ds.data->GetYaxis()->GetBinLowEdge(_y_ll) + _y_bin_width) +
                                    "] (x-units:y-units)"
                                );
                                _y_ll += _y_bin_width;
                            }
                        }
                    }
                }
            }
            else throw std::runtime_error("cannot set data ranges : TH3D is not supported yet");

            // now sanity checks on the range
            double _last = - std::numeric_limits<double>::infinity();
            for (auto& r : _current_ds.brange) {
                if (r.first > _last) _last = r.first;
                else throw std::runtime_error("illegal ranges in \"fit-range-x\" detected, did you specify them in ascending order?");
                if (r.second > _last) _last = r.second;
                else throw std::runtime_error("illegal ranges in \"fit-range-x\" detected, did you specify them in ascending order?");
            }

            data.push_back(_current_ds);

            BCLog::OutDebug("data and pdf components got so far:");
            this->DumpData();
        }
    }

    /*
     * define observables
     */

    if (config.contains("observables")) {
        BCLog::OutDebug("Parsing 'observables' section");
        for (auto& el : config["observables"].items()) {
            // sanity checks
            auto _expr = el.value()["TFormula"].get<std::string>();
            TFormula _tformula(el.key().c_str(), _expr.c_str());
            // the following check will always fail with ROOT < 6.12/04
            // The following fixes are essential for this code to work:
            // https://root.cern.ch/doc/v612/release-notes.html#histogram-libraries
            if (!_tformula.IsValid() or _tformula.GetNpar() < 1 or _tformula.GetNdim() > 0) {
                throw std::runtime_error("invalid observables TFormula given");
            }
            // loop over number of parameters in formula
            for (int p = 0; p < _tformula.GetNpar(); ++p) {
                std::string parname = _tformula.GetParName(p);
                bool exists = false;
                // look for BAT internal parameter index
                for (unsigned int idx = 0; idx < this->GetNParameters(); ++idx) {
                    if (this->GetParameters().At(idx).GetName() == parname) {
                        exists = true;
                        // we do this here because we'd rather save BAT's internal parameter index
                        // to access it later in HMixFit::CalculateObservables
                        _tformula.SetParName(p, std::to_string(idx).c_str());
                        break;
                    }
                }
                // throw exception if parameter is undefined
                if (!exists) throw std::runtime_error(
                    "fit parameter '" + parname + "' not found, is it defined in \"parameters\"?"
                );
            }

            // if requested substitute all parameters with :
            // <integral-range> X <parameter> in formula
            if (el.value().contains("multiply-fit-parameter-by-pdf-integral")) {
                BCLog::OutDebug("asked to scale observable '" + el.key() + "' with PDF integral");
                // get range
                std::vector<std::pair<double,double>> _scale_range;
                if (el.value()["multiply-fit-parameter-by-pdf-integral"].contains("range"))
                    _scale_range = utils::CheckAndStoreRanges(el.value()["multiply-fit-parameter-by-pdf-integral"]["range"]);
                else {
                    throw std::runtime_error(
                        "Need range to scale " + el.key() + " with pdf integral."
                    );
                }
                // find dataset 
                long unsigned int _ds_number = 0;
                if (el.value()["multiply-fit-parameter-by-pdf-integral"].contains("dataset")) {
                    std::string _dataset = utils::SafeROOTName(
                        "_" + el.value()["multiply-fit-parameter-by-pdf-integral"]["dataset"].get<std::string>()
                    );
                    for (auto _ds : this->data) {
                        std::string _ds_name = _ds.data->GetName();
                        if (_ds_name.find(_dataset) != std::string::npos) {
                            break;
                        }
                        _ds_number++;
                    }
                    if (_ds_number >= this->data.size()) {
                        throw std::runtime_error(
                            "Corresponding dataset " + _dataset + " not found."
                        );
                    }
                }
                else {
                    throw std::runtime_error(
                        "Integral range for observable " + el.key() + " needs association to dataset."
                    );
                }
                // for each parameter in formula now calculate the integral of
                // the corresponding pdf (the original histogram!) and
                // substitute [par] with ([par]*integral). all parameters are
                // already checked and exist
                BCLog::OutDetail("In observable TFormula scaling parameters : ");
                std::string _expr_ = _tformula.GetExpFormula().Data();
                BCLog::OutDetail(" ┌ original formula : " + _expr_);
                for (int p = 0; p < _tformula.GetNpar(); ++p) {
                    std::string parname = std::string("[") + _tformula.GetParName(p) + "]";
                    int idx = std::stoi(_tformula.GetParName(p));
                    BCLog::OutDebug("considering parameter " + parname + ". Will compute integral of component nr. "
                        + std::to_string(idx) + " in data set nr. " + std::to_string(_ds_number));
                    if (!this->data[_ds_number].comp_orig[idx]) {
                        throw std::runtime_error("histogram for component nr. " + std::to_string(idx)
                            + " in data set nr. " + std::to_string(_ds_number)
                            + " is null. The problem might be in the definition of observable '"
                            + el.key() + "'. Turn on debug mode for details.");
                    }
                    double integral = utils::IntegrateHistogram(this->data[_ds_number].comp_orig[idx], _scale_range);
                    auto _pos = _expr_.find(parname);
                    auto _len = parname.size();
                    _expr_.replace(_pos, _len, Form("(%.5e*%s)", integral, parname.c_str()));
                    auto msg = (p == _tformula.GetNpar()-1 ? " └─ " : " ├─ ")
                        + parname + " -> " + Form("(%.5e*%s)", integral, parname.c_str());
                    BCLog::OutDetail(msg);
                }
                // update TFormula
                _tformula = TFormula(el.key().c_str(), _expr_.c_str());
                // need to reset parameter names because ROOT is stupid
                // FIXME: 'Error in <TFormula::SetParName>: Parameter p1 is not defined.'
                // will be thrown if there are "jumps" in the parameters indices. e.g. is the TFormula
                // is '[0]*x + [7]' then TFormula::GetNpar() will return 8, but not all parameters 0..8
                // are really defined. I don't know how to fix this.
                for (int p = 0; p < _tformula.GetNpar(); ++p) {
                    std::string parname = _tformula.GetParName(p);
                    // delete the p in front of the parameter name
                    parname.erase(0,1);
                    _tformula.SetParName(p, parname.c_str());
                }
            }

            BCLog::OutDetail("adding observable '" + el.key() + "' (\""
                + el.value().value("long-name", "") + "\" [" + el.value().value("units", "")
                + "]) with TFormula = \"" + _tformula.GetExpFormula().Data() + "\" in range = ["
                + std::to_string(el.value()["range"][0].get<double>()) + ","
                + std::to_string(el.value()["range"][1].get<double>()) + "]");

            // register observable in BAT
            this->AddObservable(
                el.key(),
                el.value()["range"][0].get<double>(),
                el.value()["range"][1].get<double>(),
                el.value().value("long-name", ""),
                "(" + el.value().value("units", "") + ")"
            );

            // save TFormula for later use in HMixFit::CalculateObservables
            obs_tformulas.emplace(el.key(), _tformula);
        }
    }
}

HMixFit::~HMixFit() {
    for (auto& h : data) {
        delete h.data;
        for (auto& hh : h.comp) delete hh.second;
    }
}

double HMixFit::LogLikelihood(const std::vector<double>& parameters) {
    double logprob = -1.*_likelihood_offset;
    // loop over datasets
    for (auto& it : data) {
        for (auto& r : it.brange) {
            // ranges are inclusive
            for (int b = r.first; b <= r.second; ++b) {
                // compute theoretical prediction for bin 'b'
                double pred = 0;
                for (auto& h : it.comp) {
                    pred += fLivetime*parameters[h.first]*h.second->GetBinContent(b);
                }
                logprob += BCMath::LogPoisson(it.data->GetBinContent(b), pred);
            }
        }
    }
    return logprob;
}

void HMixFit::CalculateObservables(const std::vector<double>& parameters) {
    // loop over registered observables
    for (unsigned int i = 0; i < this->GetNObservables(); ++i) {
        // get definition
        auto& _tf = obs_tformulas[this->GetObservable(i).GetName()];
        // set TFormula parameters value from current chain state
        for (int p = 0; p < _tf.GetNpar(); ++p) {
            // we saved the BAT internal parameter index in the TFormula parameter name
            _tf.SetParameter(p, parameters[std::stoi(std::string(_tf.GetParName(p)))]);
        }
        // evaluate the observable expression
        this->GetObservable(i) = _tf.Eval(0);
    }
}

TH1* HMixFit::GetFitComponent(std::string filename, std::string objectname,int &num_prim, TH1* tf1_hist_format) {

    // get object
    TFile _tf(filename.c_str());
    if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
    auto obj = _tf.Get(objectname.c_str());
    if (!obj) throw std::runtime_error("HMixFit::GetFitComponent(): could not find object '" +
                                        objectname + "' in file " + filename);

    if (obj->InheritsFrom(TH1::Class())) {

        auto th_orig = dynamic_cast<TH1*>(obj);
        Long64_t _nprim = 1;

        // custom stuff for LEGEND pdfs
        TObjString* _nprim_lgnd = dynamic_cast<TObjString*>(_tf.Get("number_of_primaries"));
        if (_nprim_lgnd) {
            _nprim = std::stoll(std::string(_nprim_lgnd->GetString()));
            delete _nprim_lgnd;
        }

        // custom stuff for GERDA pdfs
        TParameter<Long64_t>* _nprim_gerda = nullptr;
        auto _name_nodir = objectname.substr(objectname.find_last_of('/')+1, std::string::npos);
        if (_name_nodir.substr(0, 3) == "M1_") {
            _nprim_gerda = dynamic_cast<TParameter<Long64_t>*>(_tf.Get("NumberOfPrimariesEdep"));
        }
        else if (_name_nodir.substr(0, 3) == "M2_") {
            _nprim_gerda = dynamic_cast<TParameter<Long64_t>*>(_tf.Get("NumberOfPrimariesCoin"));
        }

        if (_nprim_gerda) {
            _nprim = _nprim_gerda->GetVal();
            delete _nprim_gerda;
        }

        // normalize!
        if (_nprim == 1) {
            BCLog::OutWarning("could not find suitable 'number of primaries' object in '"
                + filename + "', skipping normalization");
        }
        else th_orig->Scale(1./_nprim);
        
        num_prim =_nprim;
        return th_orig;
    }
    else if (obj->InheritsFrom(TF1::Class())) {

        // is there a data format? use it
        TH1* th_new = nullptr;
        if (tf1_hist_format) {
            // clone format
            th_new = dynamic_cast<TH1*>(tf1_hist_format->Clone());
            th_new->Reset();
            th_new->SetName(obj->GetName());
            th_new->SetTitle(obj->GetTitle());
        }
        else {
            throw std::runtime_error("HMixFit::GetFitComponent(): a data format is strictly needed if source is a TF1");
        }

        for (int b = 1; b < th_new->GetNcells(); ++b) {
            // we gotta find x, y, z for every cell
            int bx = 0, by = 0, bz = 0;
            th_new->GetBinXYZ(b, bx, by, bz);
            auto x = th_new->GetXaxis()->GetBinCenter(bx);
            auto y = th_new->GetYaxis()->GetBinCenter(by);
            auto z = th_new->GetZaxis()->GetBinCenter(bz);
            // now we can evaluate the TF1
            th_new->SetBinContent(b, dynamic_cast<TF1*>(obj)->Eval(x, y, z));
        }
        delete obj;
        num_prim=0;
        return th_new;
    }
    else {
        throw std::runtime_error("HMixFit::GetFitComponent(): object '" + objectname +
                                 "' in file " + filename + " isn't of type TH1 or TF1");
    }

    return nullptr;
}

void HMixFit::SetIntegrationProperties(json j) {
    if (j.is_null()) return;

    this->SetIntegrationMethod(j.value("method", BCIntegrate::BCIntegrationMethod::kIntDefault));
    this->SetCubaIntegrationMethod(j.value("cuba-method", BCIntegrate::BCCubaMethod::kCubaDefault));

    auto method = j.value("method", "kIntDefault");
    BCLog::OutDebug("Integration method set: " + method);
    auto cubamethod = j.value("cuba-method", "kCubaDefault");
    BCLog::OutDebug("Cuba integration method set: " + cubamethod);

    if (j["integrator-settings"].is_object()) {
        auto& jj = j["integrator-settings"];
        if (jj[method].is_object()) {
            if (this->GetIntegrationMethod() != BCIntegrate::BCIntegrationMethod::kIntCuba) {
                if (jj[method]["niter-min"].is_number()) {
                    this->SetNIterationsMin(jj[method]["niter-min"].get<long long>());
                }
                if (jj[method]["niter-max"].is_number()) {
                    this->SetNIterationsMax(jj[method]["niter-max"].get<long long>());
                }
                if (jj[method]["rel-precision"].is_number()) {
                    this->SetRelativePrecision(jj[method]["rel-precision"].get<double>());
                }
                if (jj[method]["abs-precision"].is_number()) {
                    this->SetAbsolutePrecision(jj[method]["abs-precision"].get<double>());
                }
            }
            else {
                if (jj[method][cubamethod].is_object()) {
                    if (jj[method][cubamethod]["niter-min"].is_number()) {
                        this->SetNIterationsMin(jj[method][cubamethod]["niter-min"].get<long long>());
                    }
                    if (jj[method][cubamethod]["niter-max"].is_number()) {
                        this->SetNIterationsMax(jj[method][cubamethod]["niter-max"].get<long long>());
                    }
                    if (jj[method][cubamethod]["rel-precision"].is_number()) {
                        this->SetRelativePrecision(jj[method][cubamethod]["rel-precision"].get<double>());
                    }
                    if (jj[method][cubamethod]["abs-precision"].is_number()) {
                        this->SetAbsolutePrecision(jj[method][cubamethod]["abs-precision"].get<double>());
                    }
                }
                auto& jjj = jj[method][cubamethod];

                auto set_base_props = [&jjj](BCCubaOptions::General& m) {
                    if (jjj["ncomp"].is_number()) m.ncomp = jjj["ncomp"].get<int>();
                    if (jjj["flags"].is_number()) m.flags = jjj["flags"].get<int>();
                };

                switch (this->GetCubaIntegrationMethod()) {

                    case BCIntegrate::BCCubaMethod::kCubaVegas : {
                        auto o = this->GetCubaVegasOptions();
                        set_base_props(o);
                        if (jjj["nstart"]   .is_number()) o.nstart    = jjj["nstart"]   .get<int>();
                        if (jjj["nincrease"].is_number()) o.nincrease = jjj["nincrease"].get<int>();
                        if (jjj["nbatch"]   .is_number()) o.nbatch    = jjj["nbatch"]   .get<int>();
                        if (jjj["gridno"]   .is_number()) o.gridno    = jjj["gridno"]   .get<int>();
                        this->SetCubaOptions(o);
                        break;
                    }
                    case BCIntegrate::BCCubaMethod::kCubaSuave : {
                        auto o = this->GetCubaSuaveOptions();
                        set_base_props(o);
                        if (jjj["nmin"]    .is_number()) o.nmin     = jjj["nmin"]    .get<int>();
                        if (jjj["nnew"]    .is_number()) o.nnew     = jjj["nnev"]    .get<int>();
                        if (jjj["flatness"].is_number()) o.flatness = jjj["flatness"].get<double>();
                        this->SetCubaOptions(o);
                        break;
                    }
                    case BCIntegrate::BCCubaMethod::kCubaDivonne : {
                        auto o = this->GetCubaDivonneOptions();
                        set_base_props(o);
                        if (jjj["key1"]        .is_number()) o.key1         = jjj["key3"]        .get<int>();
                        if (jjj["key2"]        .is_number()) o.key2         = jjj["key2"]        .get<int>();
                        if (jjj["key3"]        .is_number()) o.key3         = jjj["key1"]        .get<int>();
                        if (jjj["maxpass"]     .is_number()) o.maxpass      = jjj["maxpass"]     .get<int>();
                        if (jjj["border"]      .is_number()) o.border       = jjj["border"]      .get<double>();
                        if (jjj["maxchisq"]    .is_number()) o.maxchisq     = jjj["maxchisq"]    .get<double>();
                        if (jjj["mindeviation"].is_number()) o.mindeviation = jjj["mindeviation"].get<double>();
                        this->SetCubaOptions(o);
                        break;
                    }
                    case BCIntegrate::BCCubaMethod::kCubaCuhre : {
                        auto o = this->GetCubaCuhreOptions();
                        set_base_props(o);
                        if (jjj["key"].is_number()) o.key = jjj["key"].get<int>();
                        this->SetCubaOptions(o);
                        break;
                    }
                    case BCIntegrate::BCCubaMethod::kCubaDefault : break;
                    case BCIntegrate::BCCubaMethod::NCubaMethods : break;
                }
            }
        }
    }
}

void HMixFit::SaveHistogramsROOT(std::string filename) {
    BCLog::OutSummary("Saving histograms (ROOT) to output file " + filename);

    TFile tf(filename.c_str(), "recreate");
    for (auto& it : data) {

        tf.mkdir(it.data->GetName(), "Dataset directory");
        tf.mkdir(Form("%s/originals", it.data->GetName()), "Original histograms");

        auto dataset_dir = Form("%s:/%s", filename.c_str(), it.data->GetName());
        auto dataset_orig_dir = Form("%s/originals", dataset_dir);

        TH1* sum = nullptr;

        BCLog::OutDebug("writing histograms in " + std::string(dataset_dir));
        tf.cd(dataset_dir);
        it.data->Write("fitted_data");
        for (auto& h : it.comp) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(fLivetime*this->GetBestFitParameters()[h.first]);
            // compute total model without changing the components
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            hcopy->Write(h.second->GetName());
            delete hcopy;
        }
        sum->SetTitle("total_model");
        sum->Write("total_model");

        sum = nullptr;
        BCLog::OutDebug("writing histograms in " + std::string(dataset_orig_dir));
        tf.cd(dataset_orig_dir);
        it.data_orig->Write("fitted_data");
        for (auto& h : it.comp_orig) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(fLivetime*this->GetBestFitParameters()[h.first]);
            // compute total model without changing the components
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            hcopy->Write(h.second->GetName());
            delete hcopy;
        }
        sum->SetTitle("total_model");
        sum->Write("total_model");

        // write fit range
        BCLog::OutDebug("writing fit ranges in " + std::string(dataset_dir));
        tf.cd(dataset_dir);
        if (it.data->GetDimension() == 1) {
            BCLog::OutDebug("Writing fit range 1D dataset to file");
            TParameter<double> range_low("fit_range_lower", it.data->GetBinLowEdge(it.brange[0].first));
            TParameter<double> range_upp("fit_range_upper", it.data->GetBinLowEdge(it.brange.back().second)
                + it.data->GetBinWidth(it.brange.back().second));
            range_low.Write();
            range_upp.Write();
        }
        // for TH2 write y-range
        if (it.data->GetDimension() == 2) {
            BCLog::OutDebug("Writing fit range 2D dataset to file");
            int binx_low, binx_up, biny, binz;
            it.data->GetBinXYZ(it.brange[0].first, binx_low, biny, binz);
            it.data->GetBinXYZ(it.brange.back().second, binx_up, biny, binz);
            auto bw_x = it.data->GetXaxis()->GetBinWidth(1);
            auto bw_y = it.data->GetYaxis()->GetBinWidth(1);

            // full x-range
            TParameter<double> range_low("fit_range_lower_x", it.data->GetXaxis()->GetBinLowEdge(binx_low));
            TParameter<double> range_upp("fit_range_upper_x", it.data->GetXaxis()->GetBinLowEdge(binx_up) + bw_x);
            range_low.Write();
            range_upp.Write();

            // all y-ranges for projection
            int idx = 0;
            int pre_biny = biny, cur_biny = biny;
            for (auto r : it.brange) {
                it.data->GetBinXYZ(r.first, binx_low, cur_biny, binz);
                // detector new y-range and write it to file
                if ((cur_biny - pre_biny) > 1) {
                    TParameter<double> range_low_y(Form("fit_range_lower_y%i",idx),
                        it.data->GetYaxis()->GetBinLowEdge(biny));
                    TParameter<double> range_upp_y(Form("fit_range_upper_y%i",idx++),
                        it.data->GetYaxis()->GetBinLowEdge(pre_biny) + bw_y);
                    range_low_y.Write();
                    range_upp_y.Write();
                    biny = cur_biny;
                }
                pre_biny = cur_biny;
            }
        }

        tf.cd();
    }
}

void HMixFit::SaveHistogramsCSV(std::string folder) {

    BCLog::OutSummary("Saving histograms (CSV) in output folder " + folder);
    std::system(("mkdir -p " + folder + "/originals").c_str());

    for (auto& it : data) {

        if (it.data->GetDimension() != 1) {
            BCLog::OutDebug(Form("HMixFit::SaveHistogramsCSV(): %s is not a 1D histogram, skipping.", it.data->GetName()));
        }

        std::ofstream fout(Form("%s/%s.csv", folder.c_str(), it.data->GetName()));

        std::vector<TH1*> comps;
        TH1* sum = nullptr;

        for (auto& h : it.comp) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(fLivetime*this->GetBestFitParameters()[h.first]);
            // compute total model without changing the components
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            comps.push_back(hcopy);
        }
        sum->SetName("total_model");

        fout << "xunit, xunit_low, fitted_data, total_model";
        for (auto c : comps) fout << ", " << c->GetName();
        fout << ", 1sig_p, 1sig_m, 2sig_p, 2sig_m, 3sig_p, 3sig_m, norm_pois_res";
        fout << '\n';

        for (int b = 0; b <= it.data->GetNbinsX()+1; ++b) {
            fout << it.data->GetBinCenter(b)
                 << ", " << it.data->GetBinLowEdge(b)
                 << ", " << it.data->GetBinContent(b)/it.data->GetBinWidth(b)
                 << ", " << sum->GetBinContent(b)/sum->GetBinWidth(b);
            for (auto c : comps) {
                fout << ", " << c->GetBinContent(b)/c->GetBinWidth(b);
            }

            auto sig1 = utils::smallest_poisson_interval(0.682, sum->GetBinContent(b));
            auto sig2 = utils::smallest_poisson_interval(0.954, sum->GetBinContent(b));
            auto sig3 = utils::smallest_poisson_interval(0.997, sum->GetBinContent(b));

            fout << ", " << sig1.second/sum->GetBinWidth(b)
                 << ", " << sig1.first/sum->GetBinWidth(b)
                 << ", " << sig2.second/sum->GetBinWidth(b)
                 << ", " << sig2.first/sum->GetBinWidth(b)
                 << ", " << sig3.second/sum->GetBinWidth(b)
                 << ", " << sig3.first/sum->GetBinWidth(b);

            fout << ", " << utils::normalized_poisson_residual(sum->GetBinContent(b), it.data->GetBinContent(b));

            fout << '\n';
        }
        for (auto c : comps) delete c;
        comps.clear();
        fout.close();

        std::ofstream fout_orig(Form("%s/originals/%s.csv", folder.c_str(), it.data->GetName()));
        sum = nullptr;

        for (auto& h : it.comp_orig) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(fLivetime*this->GetBestFitParameters()[h.first]);
            // compute total model without changing the components
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            comps.push_back(hcopy);
        }
        sum->SetName("total_model");

        fout_orig << "xunit, xunit_low, fitted_data, total_model";
        for (auto c : comps) fout_orig << ", " << c->GetName();
        fout_orig << ", 1sig_p, 1sig_m, 2sig_p, 2sig_m, 3sig_p, 3sig_m, norm_pois_res";
        fout_orig << '\n';

        // it.data_orig->Rebin(15);
        // for (auto c : comps) c->Rebin(15);
        // sum->Rebin(15);

        for (int b = 1; b <= it.data_orig->GetNbinsX(); ++b) {
            fout_orig << it.data_orig->GetBinCenter(b)
                      << ", " << it.data_orig->GetBinLowEdge(b)
                      << ", " << it.data_orig->GetBinContent(b)
                      << ", " << sum->GetBinContent(b);
            for (auto c : comps) {
                fout_orig << ", " << c->GetBinContent(b);
            }

            auto sig1 = utils::smallest_poisson_interval(0.682, sum->GetBinContent(b));
            auto sig2 = utils::smallest_poisson_interval(0.954, sum->GetBinContent(b));
            auto sig3 = utils::smallest_poisson_interval(0.997, sum->GetBinContent(b));

            fout_orig << ", " << sig1.second
                      << ", " << sig1.first
                      << ", " << sig2.second
                      << ", " << sig2.first
                      << ", " << sig3.second
                      << ", " << sig3.first;

            fout_orig << ", " << utils::normalized_poisson_residual(sum->GetBinContent(b), it.data_orig->GetBinContent(b));
            fout_orig << '\n';
        }
        for (auto c : comps) delete c;
        comps.clear();
        fout_orig.close();
    }
}

void HMixFit::WriteResultsTree(std::string filename) {
    BCLog::OutSummary("Writing Parameters to output Tree in " + filename);
    TFile tf(filename.c_str(), "recreate");
    // define parameters
    std::string par_name;
    double marg_mode;
    double marg_qt16, marg_qt84, marg_qt90;
    double glob_mode, glob_mode_error;
    // build results tree
    TTree tt("fit_par_results", "Results of the fitting procedure");
    tt.Branch("par_name",          &par_name);
    tt.Branch("marg_mode",         &marg_mode,       "marg_mode/D");
    tt.Branch("marg_quantile_16",  &marg_qt16,       "marg_quantile_16/D");
    tt.Branch("marg_quantile_84",  &marg_qt84,       "marg_quantile_84/D");
    tt.Branch("marg_quantile_90",  &marg_qt90,       "marg_quantile_90/D");
    tt.Branch("glob_mode",         &glob_mode,       "glob_mode/D");
    tt.Branch("glob_mode_error",   &glob_mode_error, "glob_mode_error/D");

    for (unsigned int p = 0; p < this->GetNVariables(); p++) {
        BCLog::OutDebug("Writing Parameter " + this->GetVariable(p).GetName() + " to tree");
        bool isfixed = (p < this->GetNParameters() and this->GetParameter(p).Fixed());
        par_name = std::string(this->GetVariable(p).GetName().data());
        if (!isfixed) {
            auto bch_marg = this->GetMarginalized(p);
            marg_mode = bch_marg.GetLocalMode();
            marg_qt16 = bch_marg.GetQuantile(0.16);
            marg_qt84 = bch_marg.GetQuantile(0.84);
            marg_qt90 = bch_marg.GetQuantile(0.90);
            glob_mode = p < this->GetNParameters() ? this->GetBestFitParameters()[p] : 0;
            glob_mode_error = p < this->GetNParameters() ? this->GetBestFitParameterErrors()[p] : 0;
        }
        else {
            marg_mode = this->GetParameter(p).GetFixedValue();
            marg_qt16 = 0;
            marg_qt84 = 0;
            marg_qt90 = 0;
            glob_mode = this->GetParameter(p).GetFixedValue();
            glob_mode_error = 0;
        }
        tt.Fill();
    }
    tt.Write();

    // calculate integral of raw histogram in fit range
    // from the original histogram, *not* the rebinned one!
    for (auto ds : data) {
        std::string comp_name;
        double orig_range, orig_bi;
        double best_range, best_bi;
        double bestErr_range, bestErr_bi;
        double marg_range, marg_bi;
        double qt16_range, qt16_bi;
        double qt84_range, qt84_bi;
        double qt90_range, qt90_bi;

        TTree ttds(Form("counts_%s", ds.data->GetName()), "counts in selected regions for each parameter");
        ttds.Branch("comp_name",               &comp_name);
        ttds.Branch("fit_range_orig",          &orig_range,    "fit_range_orig/D");
        ttds.Branch("fit_range_glob_mode",     &best_range,    "fit_range_glob_mode/D");
        ttds.Branch("fit_range_glob_mode_err", &bestErr_range, "fit_range_glob_mode_err/D");
        ttds.Branch("fit_range_marg_mod",      &marg_range,    "fit_range_marg_mod/D");
        ttds.Branch("fit_range_qt16",          &qt16_range,    "fit_range_qt16/D");
        ttds.Branch("fit_range_qt84",          &qt84_range,    "fit_range_qt84/D");
        ttds.Branch("fit_range_qt90",          &qt90_range,    "fit_range_qt90/D");
        ttds.Branch("bi_range_orig",           &orig_bi,       "bi_range_orig/D");
        ttds.Branch("bi_range_glob_mode",      &best_bi,       "bi_range_glob_mode/D");
        ttds.Branch("bi_range_glob_mode_err",  &bestErr_bi,    "bi_range_glob_mode_err/D");
        ttds.Branch("bi_range_marg_mode",      &marg_bi,       "bi_range_marg_mode/D");
        ttds.Branch("bi_range_qt16",           &qt16_bi,       "bi_range_qt16/D");
        ttds.Branch("bi_range_qt84",           &qt84_bi,       "bi_range_qt84/D");
        ttds.Branch("bi_range_qt90",           &qt90_bi,       "bi_range_qt90/D");

        for (auto c : ds.comp_orig) {
            comp_name = std::string(this->GetVariable(c.first).GetName().data());
            auto ch = c.second;
            bool isfixed   = this->GetParameter(c.first).Fixed();
            double best    = isfixed ? this->GetParameter(c.first).GetFixedValue() : this->GetBestFitParameters()[c.first];
            double bestErr = isfixed ? 0 : this->GetBestFitParameterErrors()[c.first];
            auto bch_marg  = this->GetMarginalized(c.first);
            orig_range = 0;
            for (auto r : ds.brange) {
                // this trick here is needed to find the right bin indices
                auto bl_ = ch->FindBin(ds.comp[c.first]->GetBinLowEdge(r.first));
                auto bu_ = ch->FindBin(ds.comp[c.first]->GetBinLowEdge(r.second) + ds.comp[c.first]->GetBinWidth(r.second));
                orig_range += ch->Integral(bl_, bu_);
            }
            best_range    = orig_range*best;
            bestErr_range = orig_range*bestErr;
            marg_range    = orig_range*(isfixed ? this->GetParameter(c.first).GetFixedValue() : bch_marg.GetLocalMode());
            qt16_range    = isfixed ? 0 : orig_range*bch_marg.GetQuantile(0.16);
            qt84_range    = isfixed ? 0 : orig_range*bch_marg.GetQuantile(0.84);
            qt90_range    = isfixed ? 0 : orig_range*bch_marg.GetQuantile(0.90);
            orig_bi = 0., best_bi = 0., bestErr_bi = 0.;
            if (!ds.data_orig->InheritsFrom(TH2::Class())) {
                std::vector<int> bins = { // BI window
                    ch->FindBin(1930), ch->FindBin(2099),
                    ch->FindBin(2109), ch->FindBin(2114),
                    ch->FindBin(2124), ch->FindBin(2190)
                };
                bool bad_bins = false;
                for (auto b : bins) {
                    bad_bins = ch->IsBinUnderflow(b) or ch->IsBinOverflow(b);
                }
                if (bad_bins) continue;
                orig_bi += ch->Integral(bins[0], bins[1]);
                orig_bi += ch->Integral(bins[2], bins[3]);
                orig_bi += ch->Integral(bins[4], bins[5]);
                best_bi    = orig_bi*best;
                bestErr_bi = orig_bi*bestErr;
                marg_bi    = orig_bi*(isfixed ? this->GetParameter(c.first).GetFixedValue() : bch_marg.GetLocalMode());
                qt16_bi    = isfixed ? 0 : orig_bi*bch_marg.GetQuantile(0.16);
                qt84_bi    = isfixed ? 0 : orig_bi*bch_marg.GetQuantile(0.84);
                qt90_bi    = isfixed ? 0 : orig_bi*bch_marg.GetQuantile(0.90);
            }
            ttds.Fill();
        }
        ttds.Write();
    }
}

void HMixFit::DumpData() {
    for (auto& it : data) {
        std::ostringstream addr;
        addr << " (" << it.data << "), orig (" << it.data_orig << ")";
        BCLog::OutDebug(it.data->GetName() + addr.str());
        addr.str(""); addr.clear();
        size_t i = 0;
        for (auto& comp : it.comp) {
            addr << " (" << comp.second << "), orig (" << it.comp_orig[comp.first] << ")";
            auto msg = (i == it.comp.size()-1 ? " └─ [" : " ├─ [")
                + std::to_string(comp.first) + "] " + std::string(comp.second->GetName())
                + addr.str();
            BCLog::OutDebug(msg);
            addr.str(""); addr.clear();
            i++;
        }
    }
}

double HMixFit::GetFastPValue(const std::vector<double>& parameters, long niter) {

    BCLog::OutSummary("Computing fast (corrected) p-value with " + std::to_string(niter) + " iterations");

    if (parameters.size() != this->GetParameters().Size()) {
        throw std::runtime_error("GetFastPValue: input number of parameters does not match the number of defined model parameters");
    }

    TRandom3 rnd(0);
    std::vector<unsigned> observed_comb;
    std::vector<double>   expected_comb;
    int nbins_comb = 0;

    for (auto& it : data) {

        std::vector<unsigned> observed;
        std::vector<double>   expected;
        int nbins = 0;

        // compute total model
        TH1* sum = nullptr;
        for (auto& h : it.comp) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(fLivetime*parameters[h.first]);
            // compute total model
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            delete hcopy;
        }
        for (auto& r : it.brange) {
            for (int b = r.first; b <= r.second; ++b) {
                observed.push_back(it.data->GetBinContent(b));
                observed_comb.push_back(it.data->GetBinContent(b));
                expected.push_back(sum->GetBinContent(b));
                expected_comb.push_back(sum->GetBinContent(b));
                nbins++;
                nbins_comb++;
            }
        }
        delete sum;

        double p = BCMath::CorrectPValue(
            BCMath::FastPValue(observed, expected, niter, rnd.GetSeed()),
            this->GetNFreeParameters(), nbins
        );

        BCLog::OutSummary(Form("p-value for data set '%s' => %g", it.data->GetName(), p));
    }

    double pvalue = BCMath::CorrectPValue(
        BCMath::FastPValue(observed_comb, expected_comb, niter, rnd.GetSeed()),
        this->GetNFreeParameters(), nbins_comb
    );

    return pvalue;
}

double HMixFit::Integrate(bool enable_offset) {
    // sanity check
    if (GetBestFitParameters().size() != GetNParameters() &&
        GetBestFitParameters().size() != GetNVariables()) {
        BCLog::OutSummary("No best fit parameters available, needed for integration. Please do the minimization first.");
        return 0.;
    }
    double _likelihood_integral = 0.;
    if (enable_offset) {
        // set offset to a reasonable value and integrate, finally reset the offset
        _likelihood_offset = this->LogLikelihood(this->GetBestFitParameters());
        _likelihood_integral = log(BCIntegrate::Integrate());
        BCLog::OutSummary("LogIntegral (LogOffset " + std::to_string(_likelihood_offset) +
            "): " + std::to_string(_likelihood_integral + _likelihood_offset));
        _likelihood_offset = 0.;
    }
    else {
      _likelihood_integral = BCIntegrate::Integrate();
    }
    return _likelihood_integral;
}

void HMixFit::PrintOptimizationSummary() {
    // ┌─┬┐
    // │ ││
    // ├─┼┤
    // └─┴┘
    BCLog::OutSummary(Form("Optimization summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("  Number of parameters:  Npar = %i", GetNParameters()));
    BCLog::OutSummary("  Best fit parameters (global mode):");
    auto logl_best = this->LogLikelihood(this->GetBestFitParameters());
    BCLog::OutSummary("  LogLikelihood value at global mode: " + std::to_string(logl_best));

    auto best = this->GetBestFitParameters();

    std::string line;
    int maxnamelength = this->GetMaximumParameterNameLength(true);
    line = "    ┌"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┬───────────────────┐";
    BCLog::OutSummary(line);
    BCLog::OutSummary(Form("    │ %-*s │ global mode       │", maxnamelength+5, "parameters"));
    line = "    ├"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┼───────────────────┤";
    BCLog::OutSummary(line);
    for (size_t i = 0; i < best.size(); ++i) {
        auto val = Form("%.*g", 2, best[i]);
        auto err = Form("%.*g", 2, this->GetBestFitParameterErrors()[i]);
        line = Form("    │ [%2i] %-*s │ %7s ± %-7s │",
                int(i),
                maxnamelength,
                this->GetVariable(i).GetName().data(),
                val, err
        );
        if (this->GetParameter(i).Fixed()) line += " (fixed)";
        BCLog::OutSummary(line);
    }
    line = "    └"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┴───────────────────┘";
    BCLog::OutSummary(line);
}

void HMixFit::PrintShortMarginalizationSummary() {
    // ┌─┬┐
    // │ ││
    // ├─┼┤
    // └─┴┘
    BCLog::OutSummary(Form("Marginalization summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("  Number of parameters:  Npar = %i", GetNParameters()));
    BCLog::OutSummary("  Parameters values (marginalized mode):");

    std::string line;
    int maxnamelength = this->GetMaximumParameterNameLength(true);

    line = "    ┌"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┬─────────────────────────────┐";
    BCLog::OutSummary(line);
    BCLog::OutSummary(Form("    │ %-*s │ marg. mode ± 1σ             │", maxnamelength+5, "parameters"));
    line = "    ├"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┼─────────────────────────────┤";
    BCLog::OutSummary(line);
    for (size_t i = 0; i < this->GetNVariables(); ++i) {
        bool isfixed = false;
        if (i < this->GetNParameters() and this->GetParameter(i).Fixed()) isfixed = true;
        auto val  = Form("%.*g", 2, isfixed ? this->GetParameter(i).GetFixedValue() : this->GetMarginalized(i).GetLocalMode());
        auto err1 = Form("%.*g", 2, isfixed ? 0 : this->GetMarginalized(i).GetLocalMode() - this->GetMarginalized(i).GetQuantile(0.16));
        auto err2 = Form("%.*g", 2, isfixed ? 0 : this->GetMarginalized(i).GetQuantile(0.84) - this->GetMarginalized(i).GetLocalMode());
        auto qt90 = Form("%.*g", 2, isfixed ? 0 : this->GetMarginalized(i).GetQuantile(0.90));
        if (isfixed) {
            line = Form("    │ [%2i] %-*s │ %7s (fixed)             │",
                int(i),
                maxnamelength,
                this->GetVariable(i).GetName().data(),
                val
            );
        }
        else if (this->GetMarginalized(i).GetLocalMode() > this->GetMarginalized(i).GetQuantile(0.84) or
                 this->GetMarginalized(i).GetQuantile(0.16) > this->GetMarginalized(i).GetLocalMode() ) {
            line = Form("    │ [%2i] %-*s │ < %7s (90%%C.L.)         │",
                int(i),
                maxnamelength,
                this->GetVariable(i).GetName().data(),
                qt90
             );
        }
        else {
            line = Form("    │ [%2i] %-*s │ %7s + %-7s - %-7s │",
                int(i),
                maxnamelength,
                this->GetVariable(i).GetName().data(),
                val, err2, err1
            );
        }
        BCLog::OutSummary(line);
    }
    line = "    └"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┴─────────────────────────────┘";
    BCLog::OutSummary(line);
}
