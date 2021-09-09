---
title: Config
linktitle: Pysotope config
type: project
order: 5
---

The pysotope package is highly configurable using javascript object notation files for storing isotope system configuration. 
These files contain all information about how to read the data files, which filetype plugin to use, spike-standard compostion etc.
Here an example configuration file is presented with comments on each parameter. Comments are not allowed in the regular json files used for data reduction.

<!--more-->


A version of this file is contained in the pysotope github repository along with configurations for running Cr and Cd isotopes on IsotopX Phoenix and Isoprobe-T mass thermal ionization mass spectrometers (TIMS), and Cr on the ThermoFisher Neptune multicollector inductively coupled plasma mass spectrometer (MC-ICP-MS).



[github.com/tarting/pysotope](https://github.com/tarting/pysotope)


```javascript
{ 
  // Description of intended use of the file, e.g. spike and standard 
  // lot numbers could be annotated here, as well as where and when 
  // spike and standard were calibrated.
  // It's reccomended that creator info/contact is added as well.
  "description": "Data reduction for Cr isotopes, calibrated for 50Cr-54Cr spike pair, with reporting of 53Cr/52Cr ratio relative to NIST SRM 3112a, on the Phoenix and Isoprobe-T TIMS at the Geology Section, Department of Geoscience and Natural Resources, University of Copenhagen",

  // The version parameter is used by the file spec_read to see if the file spec
  // is compatible with the current version of pysotope.
  "version": 2,
  // The element symbol, used for generating lookup text for ratios.
  "element" : "Cr",

  // List of column in the data files, must be formatted as:
  //   "mass_number" : column_index
  // Mass number must be enclosed in double qoutes: "49"
  // Index must be an integer: 0
  // The index starts at 0 for the first column.
  "cycle_columns" : 
    { "49" : 0, 
      "50" : 1, 
      "51" : 2, 
      "52" : 3, 
      "53" : 4, 
      "54" : 5, 
      "56" : 6 
    },
 
  // List of fractions n to report in the final data table, formatted as:
  //    ["denominator", ["numerator 1", "numerator 2", ... , "numerator n"]]
  // Mass number for denom. and numer. must be enclosed in double quotes, and
  // may not contain any text except for the mass number.
  "report_fracs" : 
    ["52", ["50", "53", "54"]],

  // Follows the same format as 'report_fracs', however may only contain
  // three numerators.
  // The denominator and first numerator should refer to the spiked isotopes,
  // And the two second numerators should be the reference isotopes.
  "reduce_fracs" : 
    ["54", ["50", "52", "53"]],
    
  // reporting of relative isotope composition, following the format:
  //   {"label 1" : ["numerator 1", "denominator", factor, include_unfiltered],
  //    "label 2" : ["numerator 2", "denominator", factor, include_unfiltered]}
  // The label may be any relevant text e.g. as shown below.
  // Numerator and denominator must be an isotope reference in double quotes e.g.
  //    "53Cr"
  // The factor is an integer multiplied to the deviation from the standard.
  //   e.g. 1000 for delta value and 10000 for epsilon value.
  // Include_unfiltered can be either true or false (case sensitive) 
  // given true statistics for non-outlier rejected data will be produced in 
  // addition to outlier-rejected statistics. false will only report outlier rejected.
  "rel_report":
    {
        "d53Cr_SRM3112a" : ["53Cr", "52Cr", 1000, true]
    },

  // Dict of used isotopes:
  //   "analyte_mass" : [["ref_mass", "interfering_element"], ...]
  // A list of used isotopes and their interferences.
  // A reference isotope must be referenced for each interfering elements.
  // An empty list [] denotes no interference correction.
  // The natural or reference ratio between the "analyte_mass" and "ref_mass"
  // of the interfering element must be denoted in "nat_ratios"
  "used_isotopes" :
    { 
      "50" : [["51", "V"], ["49", "Ti"]],
      "52" : [],
      "53" : [],
      "54" : [["56", "Fe"]]    
    },

  // A complete list of atomic masses of all isotopes used in the data reduction, 
  // incuding those of interfering elements.
  // Labels must be formatted as below with atomic number before element symbol.
  "masses" :
    { "50Cr" : 49.946049,
      "52Cr" : 51.940512,
      "53Cr" : 52.940653,
      "54Cr" : 53.938885,
      "49Ti" : 48.947871,
      "50Ti" : 49.944792,
      "51V"  : 50.943964,
      "50V"  : 49.947163,
      "54Fe" : 53.939613,
      "56Fe" : 55.934941},

  // Reference ratios for interfering elements.
  // labels must be two isotope labels, with the interfering isotope as the denominator
  // and reference isotope as numerator, these must be separated by a forward slash: /
  "nat_ratios" :
    { "49Ti/50Ti" :   1.0185185185,
      "51V/50V"   : 399.0,
      "56Fe/54Fe" :  15.698587},

  // Spike composition
  // Must contain:
  // - All isotopes of analyte element listed as mass fraction.
  // - A field labeled "amu", the atomic mass of the spike composition
  // - All ratios for the spike which can be generated by the report_fracs and reduce_fracs.
  "spike" :
    { "name" : "RF_UCPH_Cr",
      "50Cr" : 0.5379706252508730,
      "52Cr" : 0.0350570627662604,
      "53Cr" : 0.0108693041097927,
      "54Cr" : 0.4161030078730740
    },

  // Reference standard composition
  // Adheres to the same specifications as the spike compositon, but additionaly contains:
  // - An optional "name" field with a standard name, just for reference.
  // - Any ratio generated by the "rel_report" field.
  "standard" :
    { "name" : "SRM3112a",
      "50Cr" : 0.0434502788351757,
      "52Cr" : 0.8378541590693160,
      "53Cr" : 0.0950600848660935,
      "54Cr" : 0.0236354772294146
    },

  // Parameters for otlier rejection
  // - iqr_limit:    Multiple of the inter quartile range outside of which a 
  //                 point is considered an outlier
  // - max_fraction: Maximum fraction of rows rejected.
  "outlier_rejection":
    { "iqr_limit" : 1.5,
      "max_fraction" : 0.05},

  // Initial guess for inversion function
  // Mass bias factors an lambda value.
  // [alpha_nat, beta_ins, lambda]
  "initial_parameters":
    [ 0.003,-0.06,0.8 ],
  
  // Info for parsing dates reported in data files
  // See strftime.org or python 3 datetime for explanation.
  "date":{
    "field" : ["PARAMETERS", "Analysis Start Time"],
    "report_format" : "%Y-%m-%d %H:%M:%S"},
  
  // Filetype plugin
  // This is the import path to a python module able to read the specified
  // datatype. These can be installed by putting a copy of the plugin in the 
  // plugin directory under the pysotope code folder. 
  // The filetype plugin must expose a single function:
  //    read_file(file_path: str, file_spec: dict) -> dict
  // Where dict contains at a key "CYCLES" with a list of each analytical cycle
  // needed for data reduction.
  "filetype_plugin":
        "pysotope.plugins.filetype_xls",

  // Data file specification.
  // Specific for the provided xls reader, as it outputs all sheets as a single
  // csv file. Configured to read data files from the IonVantage software suite.
  "file_spec":
        {
        // This is a parameter lookup definition
        // I.e. returns a list of key value pairs
        "PARAMETERS": [ // Label in output dict
            "params", // param lookup identifier
            {
                // Search string data collection starts at the first
                // match if search_string.
                // In this case it is the name of the excel sheet.
                "start_string": "TUNING PARAMETERS", 
                // The algorightm searches for this value in the csv
                // and returns the value in the cell left of the match.
                // The params file has a colon (:) at the end of these labels
                // it is skipped in the definition.
                "labels": [ // list of lookup labels
                    "Method Name",
                    "Analysis Start Time",
                    "Cycles per Block",
                    "Number of Blocks",
                    "No. of Sequences",
                    "Sample/Vial Number",
                    "Sample ID",
                    "User name",
                    "Serial Number",
                    "XL Doc Name",
                    "Temperature (°C)",
                    "Temp Error (°C)",
                    "Current (Amps)"]
            }
        ],
        // Table lookup
        "CYCLES": [ // Label for output dict
            "table", // table lookup identifier
            { // Dictionary defining the data lookup
                // First row lookup string
                "start_string": "Cycle,Time", 
                // After last column lookup, skips any 
                // matches prior to start string
                "end_string": "\nBLOCK Data", 
                "first_col": 2, // First column to include (starts at 0)
                "n_columns": 7, // Number of columns to include
                // Number of rows to skip (value of 0 includes the row 
                // matced by the start_string)
                "skip_rows": 1 
            },
            {
                "search_string": ",Function:", // Search string to match for label row
                "first_col": 2, // Number of columns to include
                "n_columns": 7 // Number of columns to include
            }
        ],
        // Another table lookup
        "MONITORING": [
            "table",
            {
                "start_string": "FAC Temp",
                "end_string": "TUNING PARAMETERS",
                "first_col": 7,
                "n_columns": 10,
                "skip_rows": 1
            },
            {
                "search_string": "FAC Temp",
                "first_col": 7,
                "n_columns": 10
            }
        ]
    }
}
```
