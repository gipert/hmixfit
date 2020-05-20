/* GerdaFitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#include "GerdaFitter.hh"

// BAT
#include "BAT/BCMath.h"
#include "BAT/BCAux.h"
#include "BAT/BCTF1Prior.h"
#include "BAT/BCTH1Prior.h"

// ROOT
#include "TH2.h"
#include "TF1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TParameter.h"

GerdaFitter::GerdaFitter(json outconfig) : config(outconfig) {
    // open log file
    auto outdir = config["output-dir"].get<std::string>();
    std::system(("mkdir -p " + outdir).c_str());
    auto prefix = outdir + "/gerda-fitter-" + config["id"].get<std::string>() + "-";

    BCLog::SetLogLevelFile(BCLog::debug);
    BCLog::OpenLog(prefix + "output.log");
    BCLog::OutSummary("Saving results in " + outdir);
    BCLog::SetLogLevelScreen(config.value("logging", BCLog::summary));

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

        // assign prior
        if (el.value().contains("prior")) {
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
                TF1 _tformula = this->ParseTFormula(
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
                prior_hist->SetDirectory(nullptr);
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
        TFile _tf(el.key().c_str());
        if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + el.key());

        // loop over requested histograms in file
        for (auto& elh : el.value().items()) {
            // get rebin factor before reading all histograms
            auto rebin_x = elh.value().value("rebin-factor-x", 1) * elh.value().value("rebin-factor", 1);
            auto rebin_y = elh.value().value("rebin-factor-y", 1);
            if (rebin_x != 1)
                BCLog::OutDetail("using rebin X factor = " + std::to_string(rebin_x) + " for dataset '" + elh.key() + "'");
            if (rebin_y != 1)
                BCLog::OutDetail("using rebin Y factor = " + std::to_string(rebin_y) + " for dataset '" + elh.key() + "'");

            BCLog::OutDebug("getting data histogram '" + elh.key() + "' from " + el.key());
            auto th = dynamic_cast<TH1*>(_tf.Get(elh.key().c_str()));
            if (!th) throw std::runtime_error("could not find object '" + elh.key() + "' in file " + el.key());
            th->SetDirectory(nullptr);

            // rebin if requested
            th->Rebin(rebin_x);
            if ((rebin_y != 1) and (_tf.Get(elh.key().c_str())->InheritsFrom(TH2::Class()))) {
                dynamic_cast<TH2*>(_tf.Get(elh.key().c_str()))->RebinY(rebin_y);
            }
            auto basename = el.key().substr(
                    el.key().find_last_of('/')+1,
                    el.key().find_last_of('.') - el.key().find_last_of('/') - 1);
            basename += "_" + std::string(th->GetName());
            th->SetName(this->SafeROOTName(basename).c_str());

            dataset _current_ds;
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
                            auto thh = this->GetFitComponent(filename, elh.key(), _current_ds.data, rebin_x, rebin_y);
                            // add it with weight
                            if (!sum) {
                                sum = thh;
                                sum->SetDirectory(nullptr); // please do not delete it when the TFile goes out of scope
                                sum->Scale(p.value().get<double>());
                            }
                            else sum->Add(thh, p.value().get<double>());
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
                        auto thh = this->GetFitComponent(filename, elh.key(), _current_ds.data, rebin_x, rebin_y);
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
                        auto thh = this->GetFitComponent(it["root-file"].get<std::string>(), obj_name, _current_ds.data, rebin_x, rebin_y);
                        _current_ds.comp.insert({comp_idx, thh});
                    }
                    // it's a explicit TFormula
                    else if (iso.value().contains("TFormula")) {
                        BCLog::OutDebug("TFormula expression detected");
                        if (!iso.value()["TFormula"].is_string()) {
                            throw std::runtime_error("The \"TFormula\" key must be a string");
                        }
                        auto _tfunc = this->ParseTFormula(
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
                            _current_ds.data->GetNbinsX(),
                            _current_ds.data->GetXaxis()->GetXmin(),
                            _current_ds.data->GetXaxis()->GetXmax()
                        );
                        for (int b = 1; b < thh->GetNbinsX(); ++b) {
                            thh->SetBinContent(b, _tfunc.Eval(thh->GetBinCenter(b)));
                        }
                        thh->SetDirectory(nullptr);
                        thh->SetName(this->SafeROOTName(iso.key()).c_str());
                        _current_ds.comp.insert({comp_idx, thh});
                    }
                    else { // look into gerda-pdfs database
                        BCLog::OutDebug("should be a gerda-pdfs PDF");
                        if (iso.value()["isotope"].is_string()) {
                            auto comp = sum_parts(iso.value()["isotope"]);
                            comp->SetName(this->SafeROOTName(iso.key() + "_" + std::string(comp->GetName())).c_str());
                            _current_ds.comp.insert({comp_idx, comp});
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
                                else comp->Add(sum_parts(i.key()), i.value().get<double>());

                            }
                            comp->SetName(this->SafeROOTName(iso.key() + "_" + std::string(comp->GetName())).c_str());
                            _current_ds.comp.insert({comp_idx, comp});
                        }
                        else throw std::runtime_error("unexpected entry " + iso.value()["isotope"].dump() + "found in [\"fit\"][\""
                                + el.key() + "\"][\"" + elh.key() + "\"][\"components\"][\"" + iso.key() + "\"][\"isotope\"]");

                        if (!_current_ds.comp[comp_idx]) {
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
                auto & x_range = elh.value().contains("fit-range-x")
                    ? elh.value()["fit-range-x"] : elh.value()["fit-range"];
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
                        _vxr.push_back({r[0],r[1]});
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
                        // to access it later in GerdaFitter::CalculateObservables
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
                    _scale_range = CheckAndStoreRanges(el.value()["multiply-fit-parameter-by-pdf-integral"]["range"]);
                else {
                    throw std::runtime_error(
                        "Need range to scale " + el.key() + " with pdf integral."
                    );
                }
                // find dataset 
                long unsigned int _ds_number = 0;
                if (el.value()["multiply-fit-parameter-by-pdf-integral"].contains("dataset")) {
                    std::string _dataset = SafeROOTName(
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
                // for each parameter in formula now calculate the integral of the corresponding pdf
                // and substitute [par] with ([par]*integral).
                // all parameters are already checked and exist
                BCLog::OutDetail("In observable TFormula scaling parameters : ");
                std::string _expr_ = _tformula.GetExpFormula().Data();
                BCLog::OutDetail(" ┌ original formula : " + _expr_);
                for (int p = 0; p < _tformula.GetNpar(); ++p) {
                    std::string parname = std::string("[") + _tformula.GetParName(p) + "]";
                    int idx = std::stoi(_tformula.GetParName(p));
                    double integral = this->IntegrateHistogram(this->data[_ds_number].comp[idx], _scale_range);
                    auto _pos = _expr_.find(parname);
                    auto _len = parname.size();
                    _expr_.replace(_pos,_len,Form("(%.5e*%s)",integral,parname.c_str()));
                    auto msg = (p == _tformula.GetNpar()-1 ? " └─ " : " ├─ ")
                        + parname + " -> " + Form("(%.5e*%s)",integral,parname.c_str());
                    BCLog::OutDetail(msg);
                }
                // update TFormula
                _tformula = TFormula(el.key().c_str(), _expr_.c_str());
                // need to reset parameter names because ROOT is stupid
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

            // save TFormula for later use in GerdaFitter::CalculateObservables
            obs_tformulas.emplace(el.key(), _tformula);
        }
    }
}

GerdaFitter::~GerdaFitter() {
    for (auto& h : data) {
        delete h.data;
        for (auto& hh : h.comp) delete hh.second;
    }
}

double GerdaFitter::LogLikelihood(const std::vector<double>& parameters) {
    double logprob = -1.*_likelihood_offset;
    // loop over datasets
    for (auto& it : data) {
        for (auto& r : it.brange) {
            // ranges are inclusive
            for (int b = r.first; b <= r.second; ++b) {
                // compute theoretical prediction for bin 'b'
                double pred = 0;
                for (auto& h : it.comp) pred += parameters[h.first]*h.second->GetBinContent(b);
                logprob += BCMath::LogPoisson(it.data->GetBinContent(b), pred);
            }
        }
    }
    return logprob;
}

void GerdaFitter::CalculateObservables(const std::vector<double>& parameters) {
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

TH1* GerdaFitter::GetFitComponent(std::string filename, std::string objectname, TH1* dataformat, int rebin_x, int rebin_y) {
    TFile _tf(filename.c_str());
    if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
    auto obj = _tf.Get(objectname.c_str());
    if (!obj) throw std::runtime_error("could not find object '" + objectname + "' in file " + filename);
    // please do not delete it when the TFile goes out of scope
    if (obj->InheritsFrom(TH1::Class())) {
        auto _th = dynamic_cast<TH1*>(obj);
        // rebin component
        _th->Rebin(rebin_x);
        if ((rebin_y != 1) and (obj->InheritsFrom(TH2::Class()))) {
            dynamic_cast<TH2*>(obj)->RebinY(rebin_y);
            BCLog::OutDebug("rebinned combonent : [" + std::to_string(rebin_x) + "," + std::to_string(rebin_y) + "]");
        }
        // sanity checks
        if (_th->GetDimension() > 2) throw std::runtime_error("TH3 are not supported yet");
        if (_th->GetNcells() != dataformat->GetNcells() or // ncells : nbins + over- and underflow
            _th->GetDimension() != dataformat->GetDimension() or
            _th->GetXaxis()->GetXmin() != dataformat->GetXaxis()->GetXmin() or
            _th->GetXaxis()->GetXmax() != dataformat->GetXaxis()->GetXmax()) {
            throw std::runtime_error("histogram '" + objectname + "' in file " + filename
                + " and corresponding data histogram do not have the same number of bins and/or same ranges");
        }
        if (_th->GetDimension() == 2) {
            if (_th->GetYaxis()->GetXmin() != dataformat->GetYaxis()->GetXmin() or
                _th->GetYaxis()->GetXmax() != dataformat->GetYaxis()->GetXmax()) {
                throw std::runtime_error("histogram '" + objectname + "' in file " + filename
                + " and corresponding data histogram do not have the same y-ranges");
            }
        }
        _th->SetDirectory(nullptr);
        TParameter<Long64_t>* _nprim = nullptr;
        auto _name_nodir = objectname.substr(objectname.find_last_of('/')+1, std::string::npos);
        if (_name_nodir.substr(0, 3) == "M1_") {
            _nprim = dynamic_cast<TParameter<Long64_t>*>(_tf.Get("NumberOfPrimariesEdep"));
        }
        else if (_name_nodir.substr(0, 3) == "M2_") {
            _nprim = dynamic_cast<TParameter<Long64_t>*>(_tf.Get("NumberOfPrimariesCoin"));
        }
        if (!_nprim) {
            BCLog::OutWarning("could not find suitable 'NumberOfPrimaries' object in '"
                + filename + "', skipping normalization");
        }
        else _th->Scale(1./_nprim->GetVal());

        return _th;
    }
    else if (obj->InheritsFrom(TF1::Class())) {
        auto _th = new TH1D(obj->GetName(), obj->GetTitle(), dataformat->GetNbinsX(),
            dataformat->GetXaxis()->GetXmin(), dataformat->GetXaxis()->GetXmax());
        for (int b = 1; b < _th->GetNbinsX(); ++b) {
            _th->SetBinContent(b, dynamic_cast<TF1*>(obj)->Eval(_th->GetBinCenter(b)));
        }
        _th->SetDirectory(nullptr);
        return _th;
    }
    else {
        throw std::runtime_error("object '" + objectname + "' in file " + filename + " isn't of type TH1 or TF1");
    }

    return nullptr;
}

void GerdaFitter::SetIntegrationProperties(json j) {
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

void GerdaFitter::SaveHistograms(std::string filename) {
    TFile tf(filename.c_str(), "recreate");
    for (auto& it : data) {
        tf.mkdir(it.data->GetName());
        tf.cd(it.data->GetName());
        it.data->Write("fitted_data");

        TH1* sum = nullptr;
        for (auto& h : it.comp) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(this->GetBestFitParameters()[h.first]);
            // compute total model without changing the components
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            hcopy->Write(h.second->GetName());
            delete hcopy;
        }
        sum->SetTitle("total_model");
        sum->Write("total_model");

        // write fit range
        if (it.data->GetDimension() == 1) {
            BCLog::OutDetail("Writing fit range 1D dataset to file");
            TParameter<double> range_low("fit_range_lower", it.data->GetBinLowEdge(it.brange[0].first));
            TParameter<double> range_upp("fit_range_upper", it.data->GetBinLowEdge(it.brange.back().second)
                + it.data->GetBinWidth(it.brange.back().second));
            range_low.Write();
            range_upp.Write();
        }
        // for TH2 write y-range
        if (it.data->GetDimension() == 2) {
            BCLog::OutDetail("Writing fit range 2D dataset to file");
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

void GerdaFitter::WriteResultsTree(std::string filename) {
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

    for (unsigned int p = 0; p < this->GetNParameters(); p++) {
        BCLog::OutDebug("Writing Parameter " + this->GetVariable(p).GetName() + " to tree");
        bool isfixed = this->GetParameter(p).Fixed();
        par_name = std::string(this->GetVariable(p).GetName().data());
        if (!isfixed) {
            auto bch_marg = this->GetMarginalized(p);
            marg_mode = bch_marg.GetLocalMode();
            marg_qt16 = bch_marg.GetQuantile(0.16);
            marg_qt84 = bch_marg.GetQuantile(0.84);
            marg_qt90 = bch_marg.GetQuantile(0.90);
            glob_mode = this->GetBestFitParameters()[p];
            glob_mode_error = this->GetBestFitParameterErrors()[p];
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

        for (auto c : ds.comp) {
            comp_name = std::string(this->GetVariable(c.first).GetName().data());
            auto & ch = c.second;
            bool isfixed   = this->GetParameter(c.first).Fixed();
            double best    = isfixed ? this->GetParameter(c.first).GetFixedValue() : this->GetBestFitParameters()[c.first];
            double bestErr = isfixed ? 0 : this->GetBestFitParameterErrors()[c.first];
            auto bch_marg  = this->GetMarginalized(c.first);
            orig_range = 0.;
            for (auto r : ds.brange) {
                orig_range += ch->Integral(r.first, r.second);
            }
            best_range    = orig_range*best;
            bestErr_range = orig_range*bestErr;
            marg_range    = orig_range*(isfixed ? this->GetParameter(c.first).GetFixedValue() : bch_marg.GetLocalMode());
            qt16_range    = isfixed ? 0 : orig_range*bch_marg.GetQuantile(0.16);
            qt84_range    = isfixed ? 0 : orig_range*bch_marg.GetQuantile(0.84);
            qt90_range    = isfixed ? 0 : orig_range*bch_marg.GetQuantile(0.90);
            orig_bi = 0., best_bi = 0., bestErr_bi = 0.;
            if (!ds.data->InheritsFrom(TH2::Class())) {
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

void GerdaFitter::DumpData() {
    for (auto& it : data) {
        BCLog::OutDebug(it.data->GetName());
        size_t i = 0;
        for (auto& comp : it.comp) {
            auto msg = (i == it.comp.size()-1 ? " └─ [" : " ├─ [")
                + std::to_string(comp.first) + "] " + std::string(comp.second->GetName());
            BCLog::OutDebug(msg);
            i++;
        }
    }
}

double GerdaFitter::GetFastPValue(const std::vector<double>& parameters, long niter) {

    BCLog::OutSummary("Computing fast (corrected) p-value with " + std::to_string(niter) + " iterations");

    if (parameters.size() != this->GetParameters().Size()) {
        throw std::runtime_error("GetFastPValue: input number of parameters does not match the number of defined model parameters");
    }

    std::vector<unsigned> observed;
    std::vector<double>   expected;
    int nbins = 0;

    for (auto& it : data) {
        // compute total model
        TH1* sum = nullptr;
        for (auto& h : it.comp) {
            auto hcopy = dynamic_cast<TH1*>(h.second->Clone());
            hcopy->Scale(parameters[h.first]);
            // compute total model
            if (!sum) sum = dynamic_cast<TH1*>(hcopy->Clone());
            else      sum->Add(hcopy);
            delete hcopy;
        }
        for (auto& r : it.brange) {
            for (int b = r.first; b <= r.second; ++b) {
                observed.push_back(it.data->GetBinContent(b));
                expected.push_back(sum->GetBinContent(b));
                nbins++;
            }
        }
        delete sum;
    }

    TRandom3 rnd(0);
    double pvalue = BCMath::CorrectPValue(
        BCMath::FastPValue(observed, expected, niter, rnd.GetSeed()),
        this->GetNFreeParameters(), nbins
    );

    return pvalue;
}

double GerdaFitter::Integrate(bool enable_offset) {
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

void GerdaFitter::PrintOptimizationSummary() {
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
    int maxnamelength = this->GetMaximumParameterNameLength(false);
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

void GerdaFitter::PrintShortMarginalizationSummary() {
    // ┌─┬┐
    // │ ││
    // ├─┼┤
    // └─┴┘
    BCLog::OutSummary(Form("Marginalization summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("  Number of parameters:  Npar = %i", GetNParameters()));
    BCLog::OutSummary("  Parameters values (marginalized mode):");

    std::string line;
    int maxnamelength = this->GetMaximumParameterNameLength(false);

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

std::string GerdaFitter::SafeROOTName(const std::string orig) {
    TString torig(orig);

    for (auto& c : {'.', '-', '/', ':', '|', '+'}) torig.ReplaceAll(c, '_');
    for (auto& c : {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}) {
        if (torig[0] == c) torig[0] = 'N';
    }

    return std::string(torig.Data());
}

TF1 GerdaFitter::ParseTFormula(std::string prefix, std::string expr, double rangelow, double rangeup) {
    // if there's a list of parameters after
    if (expr.find(':') != std::string::npos) {
        auto formula = expr.substr(0, expr.find_first_of(':'));
        auto parlist = expr.substr(expr.find_first_of(':')+1, std::string::npos);
        std::vector<double> parlist_vec;
        TF1 _tformula(
            (prefix + "_prior_tf").c_str(), formula.c_str(),
            rangelow, rangeup
        );

        if (!_tformula.IsValid()) throw std::runtime_error("invalid prior TFormula given");

        // eventually set parameters, if any
        while (!parlist.empty()) {
            auto val = std::stod(parlist.substr(0, parlist.find_first_of(',')));
            parlist_vec.push_back(val);

            if (parlist.find(',') == std::string::npos) parlist = "";
            else parlist.erase(0, parlist.find_first_of(',')+1);
        }
        if ((int)parlist_vec.size() != _tformula.GetNpar()) {
            throw std::runtime_error("number of values specified in '"
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
            throw std::runtime_error("TFormula specified with prefix '"
                    + prefix + "' is parametric but no parameters were specified");
        }
        return _tformula;
    }
}

template<typename BasicJsonType>
std::vector<std::pair<double,double>> GerdaFitter::CheckAndStoreRanges(BasicJsonType& range) {
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

    if (_bad_range) throw std::runtime_error("Range is ill-defined");

    return _v_range;
}

std::vector<std::pair<int,int>> GerdaFitter::TranslateAxisRangeToBinRange(
    TH1* /*h*/,
    std::vector<std::pair<double,double>> /*x_range*/,
    std::vector<std::pair<double,double>> /*y_range*/
) {
    BCLog::OutSummary("Implement me : TranslateAxisRangeToBinRange");
    std::vector<std::pair<int,int>> _b_range;
    return _b_range;
}

double GerdaFitter::IntegrateHistogram(
    TH1* h,
    std::vector<std::pair<double,double>> x_range,
    std::vector<std::pair<double,double>> y_range
) {
    if (h->GetDimension() > 2) {
        throw std::runtime_error("IntegrateHistogram not implemeted for TH3.");
    }
    else if (h->GetDimension() == 2) {
        return this->IntegrateHistogram2D(dynamic_cast<TH2*>(h), x_range, y_range);
    }
    else if (h->GetDimension() == 1) {
        return this->IntegrateHistogram1D(h, x_range);
    }
    else {
        throw std::runtime_error("Something went wrong in IntegrateHistogram.");
    }
    return 0.;
}

double GerdaFitter::IntegrateHistogram1D(TH1* h, std::vector<std::pair<double,double>> range) {
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

double GerdaFitter::IntegrateHistogram2D(
    TH2* /*h*/,
    std::vector<std::pair<double,double>> /*x_range*/,
    std::vector<std::pair<double,double>> /*y_range*/
)
{ 
    BCLog::OutSummary("Implement me : IntegrateHistogram::TH2");
    return 0;
}
