# hmixfit ![CI](https://github.com/gipert/hmixfit/workflows/CI/badge.svg)

A fully JSON-configurable bayesian fitting engine (based on
[BAT](https://github.com/bat/bat) and
[ROOT](https://github.com/root-project/root)) for data in the form of ROOT
histograms. Taylored on GERDA data and Probability Density Functions.

### Compile and install

Requirements
 - [ROOT](https://github.com/root-project/root) ≥ v6.12/04
 - [BAT](https://github.com/bat/bat) ≥ v1.0.0 (with Cuba enabled)

Then just `PREFIX=/path/to/prefix make install`.

Alternatively, a Singularity container can be used:
```console
$ sudo singularity build hmixfit.sif Singularity.def
$ singularity exec hmixfit.sif hmixfit -h
USAGE: hmixfit [-h|--help] json-config
```

### Usage

The `hmixfit` executable acceps a JSON config file as the only argument.
Examples can be found in this repository under `config/`.

The JSON config file begins with some general settings:
```js
{
    "id" : "phIIAfterLAr",       // model name
    "logging" : "summary",       // BAT verbosity level, see manual
    "precision" : "kMedium",     // precision (number and length of Markov chains), see BAT manual
    "output-dir" : "../results", // folder with fit results
    // ...
```
settings about the global mode search algorithm:
```js
    "global-mode-search" : {
        "method" : "kOptMinuit" // see the BAT manual to learn about the other algorithms
    },
```
settings about the numerical integration needed to compute the evidence:
```js
    "integration" : {
        "enabled" : false,               // enable/disable the integration step
        "method" : "kIntCuba",           // see the BAT manual to learn about the other algorithms
        "cuba-method" : "kCubaDivonne",  // see the Cuba manual
        "integrator-settings" : {
            "kIntCuba" : {               // here you can tweak the Cuba integration settings
                "kCubaDivonne" : {       // here for the Divonne algorithm
                    "niter-max" : 1E07,
                    "niter-min" : 0,
                    "flags" : 0
                }
                "kCubaVegas" : {         // here for Vegas...
                    // ...
                }
                // ...
            }
        }
    },
    // ...
```
settings about the p-value determination
```js
    "p-value" : {
        "enabled" : false,   // enable/disable the computation
        "iterations" : 1E07  // play with this number until the p-value is stable
    },
```
and finally the fit configuration section `"fit"`, where everything about the data and
the fit components is specified in a modular fashion:
```js
    // ...
    "fit" : {
        "parameters" : { /* ... */ },  // define fit parameters globally
        "theoretical-expectations" : { /* ... */ }  // import PDFs and associated parameters
    }
}
```
Let's start with the `"parameters"` section, here the fit parameters must be defined:
```js
"parameters" : {
    "alpha-slope-bege" : {  // unique internal name
        "range" : [2E-5, 1E-4],
        "long-name" : "#alpha-model BEGe - slope",
        "units" : "cts",
        "prior" : { "histogram" : "priorfile.root:objname" }  // specify prior via external TH1
    },
    "alpha-offset-bege" : {
        "range" : [0, 1E-1],
        "long-name" : "#alpha-model BEGe - offset",
        "units" : "cts"
        "prior" : { "TFormula" : "gaus:1,10,5" }  // specify prior via TFormula
    },
    "background" : {
        "fixed" : 1234,  // parameters can be fixed to a value (not fit parameters anymore)
        "long-name" : "Background model",
        "units" : "cts"
    },
    // ...
}
```
and then associated to PDFs in the `"theoretical-expectations"` section:
```js
"theoretical-expectations" : { // takes a list of files with data histograms
    "../data/gerda-data-bkgmodel-phaseII-v04.00-lar.root" : {  // takes a list of object names in the file
        "M1_enrBEGe" : {  // this is a 1D histogram
            "gerda-pdfs" : "../data/gerda-pdfs/v2.1",  // set here the path to the gerda-pdfs, if you want
            "fit-range" : [[560, 2014], [2064, 5300]],  // note the possibility to skip regions
            "rebin-factor" : 5,  // fixed-size rebin
            "rebin-factor" : "560:10:700,700:20:900,1000:100:5300",  // support for variable binning!
            "components" : [  // here you must specify a list of PDFs you want to use
                { /* ... */ }, { /* ... */ }, // ...
            ]
        },
        "M1_enrCoax" : { /* ... */ },
        "M2_enrGe" : {  // this is a 2D histogram
            "fit-range-x" : [[560, 2014], [2064, 5300]],
            "fit-range-y" : [700, 5300],
            "rebin-factor-x" : 5,  // or just "rebin-factor" to rebin both axes
            "rebin-factor-y" : 2,
            "components" : [  // here you must specify a list of PDFs you want to use
                { /* ... */ }, { /* ... */ }, // ...
            ]
        },
        // ...
    },
    "../data/gerda-data-bkgmodel-phaseII-v04.00-raw.root" : { /* ... */ }
    // ...
}
```
the keys in the `"theoretical-expectations"` dictionary must be paths to the
files that contain histograms to be fitted (the data). Then for each of these
files the user must specify what histograms (ROOT objects) the program should
try to fit. For every data histogram a list of fit components must be provided
in the `"components"` array. The array is filled with JSON objects that can be
of multiple types.

As instance, one might want to use the GERDA PDFs distributed within
[gerda-mage-sim](https://github.com/mppmu/gerda-mage-sim) using the following
structure:
```js
{
    "gerda-pdfs" : "../data/gerda-pdfs/v2.1"  // the gerda-pdfs path might be set here to override the global one
    "part": "cables/cables_all",
    "components" : {
        "Th228-cables" : {  // this parameter name must be defined in the "parameters" section!
            "isotope" : { "Tl208-larveto" : 0.3539, "Bi212-larveto" : 1 },  // specify a mixture of isotopes
        },
        "Co60-cables" : {
            "isotope": "Co60-run68pca", // no mixture here
            // ...
        },
        // ...
    }
},
{
    "part": {  // you can also specify a mixture of parts!
        "calib/single_s2_8220" : 52183,
        "calib/single_s2_8408" : 25337,
        "calib/single_s2_8570" : 79868,
        "calib/single_s3_8220" : 55438,
        "calib/single_s3_8405" : 43433,
        "calib/single_s3_8570" : 24130
    },
    "components" : { /* ... */ }
}
```
or even provide manually a ROOT histogram:
```js
{
    "root-file" : "../data/gerda-pdfs/v2.0-RC/alphas/analytic/pdf-functions.root",
    "components" : {
        "alpha-offset" : {
            "hist-name" : "flat",
        },
        // ...
    }
},
```
or even a ROOT `TFormula` in the form `"formula:par1,par2,..."`:
```js
{
    "components" : {
        "alpha-slope" : {
            "TFormula" : "gaus:1,34,2",
        },
        // ...
    }
},
```
Last but not least, observables that depend on the model parameters only can be
defined via JSON file with the following syntax:
```js
"parameters" : {
    "2nbb-half-life-bege" : {  // unique internal name
        "TFormula": "1.13380E26/[2nbb-bege]",  // ROOT's TFormula
        "multiply-fit-parameter-by-pdf-integral" : {  // there's the possibility to multiply each parameter
                                                      // above by the pdf integral in a range:
                                                      // [2nbb-bege] -> ([2nbb-bege]*Int)
            "range" : [[10,19], [80,89]],  // range for the integral
            "dataset" : "h_data"  // dataset pdf refers to
        },
        "range" : [2E-5, 1E-4],
        "long-name" : "T_{1/2}^{2#nu} - BEGe",
        "units" : "cts",
    },
    // ...
```
Model parameters must be specified as they were a `TFormula` parameter,
enclosing their name in square brackets.
