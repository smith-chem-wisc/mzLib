# mzLib

A library for mass spectrometry projects.

[![GitHub Action Builds](https://github.com/smith-chem-wisc/mzLib/actions/workflows/dotnet.yml/badge.svg)](https://github.com/smith-chem-wisc/mzLib/actions/workflows/dotnet.yml)
[![codecov](https://codecov.io/gh/smith-chem-wisc/mzLib/branch/master/graph/badge.svg)](https://codecov.io/gh/smith-chem-wisc/mzLib)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/10000/badge.svg)](https://scan.coverity.com/projects/mzlib)
[![NuGet Badge](https://buildstats.info/nuget/mzLib)](https://www.nuget.org/packages/mzLib/)

![image](https://user-images.githubusercontent.com/16841846/113908189-df7a6e80-979b-11eb-9a2d-571a53e167ac.png)

NuGet packages are released here: https://www.nuget.org/packages/mzLib/

GitHub release tags are recoreded here: https://github.com/smith-chem-wisc/mzLib/releases

# Usage
## Reading Spectra Files
To read Thermo or mzML files, use
```
ThermoStaticData staticThermo = ThermoStaticData.LoadAllStaticData(@"spectra.raw");
ThermoDynamicData dynamicThermo = ThermoDynamicData.InitiateDynamicConnection(@"spectra.raw")
Mzml mzmlFile = Mzml.LoadAllStaticData(@"spectra.mzML");
```
Both filetypes implement the same interface that has all of the necessary functionality to interact with spectra files:
```
IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> thermoFile = new ThermoRawFile(@"spectra.RAW");
IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> mzmlFile = new Mzml(@"spectra.mzML");
```
## Loading Databases From Online Sources
```
Loaders.LoadElements("elements.dat"); // void, loads elements into static PeriodicTable class 
IEnumerable<ModificationWithLocation> unimodMods = Loaders.LoadUnimod("unimod.dat");
IEnumerable<ModificationWithLocation> uniprotMods = Loaders.LoadUniprot("uniprot.dat");
```
## Reading Protein Database Files
To read .fasta, .xml, or .xml.gz files, use 
```
List<Protein> proteins = ProteinDbLoader.LoadProteinDb("proteins.xml", generateDecoys, allKnownModifications, IsContaminant, out unknownModifications);
```
The parameters are:
* ```bool generateDecoys``` True if wish to generate decoy proteins.
* ```IDictionary<string, IList<Modification>> allKnownModifications``` Dictionary of modifications with keys that correspond to modifications in the xml file.
* ```bool IsContaminant``` True if it is a contaminant database
* ```out Dictionary<string, Modification> unknownModifications``` An auxiliary output of modifications that were in the xml file but are not known.

## Reading Modification Files
To load modifications from ptmlist formatted files use
```
IEnumerable<ModificationWithLocation> ptms = PtmListLoader.ReadMods("ptms.txt")
```
# License
Code heavily borrowed from https://github.com/dbaileychess/CSMSL and distrubuted under the appropriate license, LGPL.

