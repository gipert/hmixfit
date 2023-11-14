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

#include "HMixFit.hh"

// STL
#include <iostream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <chrono>

// BAT
#include "BAT/BCLog.h"

int main(int argc, char** argv) {

    /*
     * get command line args
     */

    std::string progname(argv[0]);

    auto usage = [&]() {
        std::cerr << "USAGE: " << progname << " [-h|--help] json-config\n";
    };

    const char* const short_opts = ":h";
    const option long_opts[] = {
        { "help",  no_argument, nullptr, 'h' },
        { nullptr, no_argument, nullptr, 0   }
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'h': // -h or --help
            case '?': // Unrecognized option
            default:
                usage();
                return 1;
        }
    }

    // extra arguments
    std::vector<std::string> args;
    for(; optind < argc; optind++){
        args.emplace_back(argv[optind]);
    }

    if (args.empty() or args.size() > 1) {usage(); return 1;}

    std::ifstream fconfig(args[0]);
    if (!fconfig.is_open()) {
        BCLog::OutError("config file " + args[0] + " does not exist");
        return 1;
    }
    json config;
    fconfig >> config;

    /*
     * main routine
     */

    HMixFit* model;
    try {
        model = new HMixFit(config);
    }
    catch(std::exception& e) {
        BCLog::OutError(e.what());
        BCLog::OutError("caught exception while initializing model, aborting...");
        return 1;
    }

      BCLog::SetLogLevelScreen(config.value("logging", BCLog::summary));

    // set precision (number of samples in Markov chain)
    model->SetPrecision(config.value("precision", BCEngineMCMC::kMedium));

    // let's be more generous (default 100)
    model->SetInitialPositionAttemptLimit(1000);

    // run MCMC and marginalize posterior w/r/t all parameters
    // An estimation of the global mode is available but without uncertainties
    auto start = std::chrono::system_clock::now();
    model->MarginalizeAll();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);
    BCLog::OutSummary("Time spent: " + std::to_string(elapsed.count()) + "s");
    model->PrintShortMarginalizationSummary();
    BCLog::OutSummary("");

    // run mode finding, by default using Minuit
    auto opt_method = BCIntegrate::kOptMinuit;
    if (config["global-mode-search"].is_object()) {
        opt_method = config["global-mode-search"].value("method", BCIntegrate::kOptMinuit);
    }
    model->FindMode(opt_method, model->GetBestFitParameters());
    model->PrintOptimizationSummary();
    BCLog::OutSummary("");

    // integration
    if (config["integration"].is_object()) {
        if (config["integration"].value("enabled", false)) {
            model->SetIntegrationProperties(config["integration"]);
            model->Integrate((bool)config["integration"].value("use-best-fit-likelihood-offset", false));
            auto logpost = model->LogProbability(model->GetBestFitParameters());
            BCLog::OutSummary("Posterior at global mode: " + std::to_string(logpost));
        }
    }

    // fast p-value computation
    if (config["p-value"].is_object()) {
        if (config["p-value"].value("enabled", false)) {
            start = std::chrono::system_clock::now();
            auto pvalue = model->GetFastPValue(model->GetBestFitParameters(), config["p-value"].value("iterations", 1E07));
            elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);
            BCLog::OutSummary("p-value = " + std::to_string(pvalue));
            BCLog::OutSummary("Time spent: " + std::to_string(elapsed.count()) + "s");
        }
    }

    // OUTPUT
    auto outdir = config["output-dir"].get<std::string>();
    std::system(("mkdir -p " + outdir).c_str());
    auto prefix = outdir + "/hmixfit-" + config["id"].get<std::string>() + "-";

    // draw parameter plot
    model->PrintParameterPlot(prefix + "parameters.pdf");
    model->PrintParameterLatex(prefix + "parameters.tex");
    model->PrintCorrelationPlot(prefix + "par-correlation.pdf");
    model->PrintCorrelationMatrix(prefix + "correlation-matrix.pdf");

    // draw/save all marginalized distributions
    model->WriteMarginalizedDistributions(prefix + "marginalized.root", "recreate");
    model->SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    model->PrintKnowledgeUpdatePlots(prefix + "know-update.pdf");
    model->SaveHistogramsROOT(prefix + "histograms.root");
    model->SaveHistogramsCSV(prefix + "histograms.csv");
    model->WriteResultsTree(prefix + "analysis.root");

    std::ofstream fcfg_copy(prefix + "config.json");
    fcfg_copy << std::setw(4) << config;
    fcfg_copy.close();

    BCLog::OutSummary("Exiting");
    // close log file
    BCLog::CloseLog();

    delete model;

    return 0;
}
