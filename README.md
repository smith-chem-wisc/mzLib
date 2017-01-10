#mzLib

A library for mass spectrometry projects.

 [![Build status](https://ci.appveyor.com/api/projects/status/d6jjrjfk8ou3waky/branch/master?svg=true)](https://ci.appveyor.com/project/stefanks/mzlib/branch/master)
 [![Build Status](https://travis-ci.org/smith-chem-wisc/mzLib.svg?branch=master)](https://travis-ci.org/smith-chem-wisc/mzLib)
 [![Coverage Status](https://coveralls.io/repos/github/smith-chem-wisc/mzLib/badge.svg?branch=master)](https://coveralls.io/github/smith-chem-wisc/mzLib?branch=master)
 [![Coverity Scan Build Status](https://scan.coverity.com/projects/10000/badge.svg)](https://scan.coverity.com/projects/mzlib)
 


Releases are here: https://www.nuget.org/packages/mzLib/

## Usage

### Populate the periodic table

To automatically populate the periodic table with current NIST values, use
```csharp
            Loaders.LoadElements("elements.dat");
```

### Chemical Formulas

Elements (and isotopes) can be combined to create chemical formulas:
```csharp
            ChemicalFormula formula = new ChemicalFormula("C2H3NO");
```

### Isotopic Distribution

Isotopic distribution of a chemical compound can be calculated from its chemical formula:
```csharp
            IsotopicDistribution dist = new IsotopicDistribution(formula);
```

### Proteomics

Peptides can be created for a sequence of amino acids

### Mass Spectrometry

Spectra can be read from raw files. 

## License
Code heavily borrowed from https://github.com/dbaileychess/CSMSL and distrubuted under the appropriate license, LGPL.
