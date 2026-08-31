# mzLib

A library for mass spectrometry projects.

[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/mzLib)](https://github.com/smith-chem-wisc/mzLib/releases/latest)
[![Build and Test](https://github.com/smith-chem-wisc/mzLib/actions/workflows/dotnet.yml/badge.svg?branch=master)](https://github.com/smith-chem-wisc/mzLib/actions/workflows/dotnet.yml)
[![codecov](https://codecov.io/gh/smith-chem-wisc/mzLib/branch/master/graph/badge.svg)](https://codecov.io/gh/smith-chem-wisc/mzLib)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/10000/badge.svg)](https://scan.coverity.com/projects/mzlib)

[![Github All Releases](https://img.shields.io/nuget/dt/mzlib)]()

![image](https://user-images.githubusercontent.com/16841846/113908189-df7a6e80-979b-11eb-9a2d-571a53e167ac.png)

NuGet packages are released here: https://www.nuget.org/packages/mzLib/

GitHub release tags are recoreded here: https://github.com/smith-chem-wisc/mzLib/releases

# Documentation

**The [mzLib wiki](https://github.com/smith-chem-wisc/mzLib/wiki) is the usage documentation, and it is
kept current.** Among other things it covers:

- [Reading mass spectrometry files](https://github.com/smith-chem-wisc/mzLib/wiki/File-Reading:-Mass-Spec)
- [Reading search-result formats](https://github.com/smith-chem-wisc/mzLib/wiki/File-Reading:-Result-Formats)
- [Reading sequence databases](https://github.com/smith-chem-wisc/mzLib/wiki/File-Reading:-Sequence-Databases)
- [Chemistry](https://github.com/smith-chem-wisc/mzLib/wiki/Chemistry) and [Mass Spectrometry](https://github.com/smith-chem-wisc/mzLib/wiki/Mass-Spectrometry)
- Omics: [base foundation](https://github.com/smith-chem-wisc/mzLib/wiki/Omics:-Base-Foundation), [digestion](https://github.com/smith-chem-wisc/mzLib/wiki/Omics:-Digestion), [modifications](https://github.com/smith-chem-wisc/mzLib/wiki/Omics:-Modifications), [fragmentation](https://github.com/smith-chem-wisc/mzLib/wiki/Omics:-Fragmentation), [decoy generation](https://github.com/smith-chem-wisc/mzLib/wiki/Omics:-Decoy-Generation), [proteomics](https://github.com/smith-chem-wisc/mzLib/wiki/Omics:-Proteomics)

The package also ships XML documentation, so mzLib types and members have tooltips in your IDE.

Usage examples used to live in this file. They were removed rather than rewritten because they had
drifted out of date - one called a `ProteinDatabaseLoader` type that no longer exists - and a second
copy of the documentation is a second thing to keep current. The wiki is the single place now.

# Contributing

Tests live in `mzLib/Test`. Run them with
`dotnet test mzLib/Test/Test.csproj --filter "Category!=ExternalService"`; the `ExternalService`
category covers tests that reach live services such as UniProt, and is deliberately kept out of the
required build so that an outage elsewhere cannot fail a pull request.

# License
Code heavily borrowed from https://github.com/dbaileychess/CSMSL and distributed under the appropriate license, LGPL.

mzLib additionally redistributes third-party data files and native libraries whose licences are
separate from its own; see [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).
