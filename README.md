# mzLib

A library for mass spectrometry projects.

 [![Build status](https://ci.appveyor.com/api/projects/status/d6jjrjfk8ou3waky/branch/master?svg=true)](https://ci.appveyor.com/project/stefanks/mzlib/branch/master)
 [![Build Status](https://travis-ci.org/smith-chem-wisc/mzLib.svg?branch=master)](https://travis-ci.org/smith-chem-wisc/mzLib)
 [![Coverage Status](https://coveralls.io/repos/github/smith-chem-wisc/mzLib/badge.svg?branch=master)](https://coveralls.io/github/smith-chem-wisc/mzLib?branch=master)
 [![Coverity Scan Build Status](https://scan.coverity.com/projects/10000/badge.svg)](https://scan.coverity.com/projects/mzlib)
 


Releases are here: https://www.nuget.org/packages/mzLib/

# Usage

## Loading Databases From Online Sources
```
Loaders.LoadElements("elements.dat"); // void, loads elements into static PeriodicTable class 
IEnumerable<ModificationWithLocation> unimodMods = Loaders.LoadUnimod("unimod.dat");
IEnumerable<ModificationWithLocation> uniprotMods = Loaders.LoadUniprot("uniprot.dat");
```

## Reading Protein Database Files
To read .fasta, .xml, or .xml.gz files, use 
```
var ok = ProteinDbLoader.LoadProteinDb("proteins.xml", generateDecoys, allKnownModifications, IsContaminant, out unknownModifications);
```
The parameters are:
* ```bool generateDecoys``` True if wish to generate decoy proteins.
* ```IDictionary<string, IList<Modification>> allKnownModifications``` Dictionary of modifications with keys that correspond to modifications in the xml file.
* ```bool IsContaminant``` True if it is a contaminant database
* ```out Dictionary<string, Modification> unknownModifications``` An auxiliary output of modifications that were in the xml file but are not known. 

## Reading Modification Files
To load modifications from ptmlist formatted files, use
```
var ptms = PtmListLoader.ReadMods("ptms.txt")
```
# License
Code heavily borrowed from https://github.com/dbaileychess/CSMSL and distrubuted under the appropriate license, LGPL.
