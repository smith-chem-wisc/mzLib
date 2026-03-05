# AGENTS.md - mzLib Development Guide

A C# class library for mass spectrometry proteomics and transcriptomics analysis.

## Build Commands

```bash
# Build entire solution
dotnet build mzLib/mzLib.sln

# Build in Release mode
dotnet build mzLib/mzLib.sln -c Release

# Build specific project
dotnet build mzLib/Proteomics/Proteomics.csproj
```

## Test Commands (NUnit)

```bash
# Run all tests
dotnet test mzLib/mzLib.sln

# Run all tests in main test project
dotnet test mzLib/Test/Test.csproj

# Run all tests in FlashLFQ test project
dotnet test mzLib/TestFlashLFQ/TestFlashLFQ.csproj

# Run tests by namespace/class (use --filter)
dotnet test mzLib/Test/Test.csproj --filter "FullyQualifiedName~TestPeptideWithSetMods"
dotnet test mzLib/Test/Test.csproj --filter "FullyQualifiedName~TestProteinDigestion"

# Run a single test by name
dotnet test mzLib/Test/Test.csproj --filter "Name=TestPeptideWithSetModsEquals"
dotnet test mzLib/Test/Test.csproj --filter "Name=TestDigestIndices"

# Run tests matching a pattern
dotnet test mzLib/Test/Test.csproj --filter "Name~Fragment"

# Combine filters (AND logic)
dotnet test mzLib/Test/Test.csproj --filter "FullyQualifiedName~TestProtein&Name~Digest"

# Run with verbose output
dotnet test mzLib/Test/Test.csproj --filter "Name=TestPeptideWithSetModsEquals" -v n
```

## Code Style

### Formatting
- Use 4 spaces for indentation (no tabs)
- Opening braces on same line for methods and control structures
- One class per file, filename matches class name
- Use `var` when type is obvious from right-hand side

### Naming Conventions
- **Classes/Interfaces**: PascalCase (`PeptideWithSetModifications`, `IBioPolymer`)
- **Methods/Properties**: PascalCase (`GetDigestProducts()`, `BaseSequence`)
- **Private fields**: camelCase with underscore prefix (`_baseSequence`)
- **Parameters/locals**: camelCase (`peptideSequence`, `fragmentIons`)
- **Constants**: PascalCase (`MaxMissedCleavages`)

### Imports
- System namespaces first, then third-party, then project namespaces
- Remove unused usings
- No global usings in this codebase

### Error Handling
- Use exceptions for exceptional conditions, not control flow
- Validate public method parameters with ArgumentNullException/ArgumentException
- Use nullable reference types where appropriate

---

## Domain Model - Critical Information

### Class Hierarchy Overview

```
IBioPolymer (parent biopolymer - undigested)
├── Protein (Proteomics namespace)
└── NucleicAcid (Transcriptomics namespace)

IBioPolymerWithSetMods (digested fragment with modifications)
├── PeptideWithSetModifications (from Protein)
└── OligoWithSetMods (from NucleicAcid)
```

### Key Interfaces

**IBioPolymer** - Represents full undigested biopolymer (protein or nucleic acid)
- `BaseSequence`: Full amino acid or nucleotide sequence
- `Accession`: Unique identifier
- `OneBasedPossibleLocalizedModifications`: Potential modification sites

**IBioPolymerWithSetMods** - Represents digested fragment with fixed modifications
- `BaseSequence`: Unmodified sequence of the fragment
- `FullSequence`: Sequence with modification annotations (e.g., `PEPTK[Acetyl]IDE`)
- `Parent`: Reference to source IBioPolymer
- `AllModsOneIsNterminus`: Dictionary of applied modifications (see indexing below)

---

## CRITICAL: Indexing Conventions

### Sequence Indexing (Zero-Based)
`BaseSequence` is a string with standard zero-based indexing:
```csharp
string seq = "PEPTIDE";
char first = seq[0];  // 'P'
char last = seq[6];   // 'E'
```

### Residue Position (One-Based)
Residue positions in the parent protein/nucleic acid are ONE-BASED:
```csharp
peptide.OneBasedStartResidue  // First residue position in parent (1-indexed)
peptide.OneBasedEndResidue    // Last residue position in parent (1-indexed)
```

### Modification Dictionary Indexing (Special)
`AllModsOneIsNterminus` uses a SPECIAL indexing scheme:

| Index | Position |
|-------|----------|
| 1 | N-terminus (5' for RNA) |
| 2 | Residue 1 (first amino acid/nucleotide) |
| 3 | Residue 2 |
| ... | ... |
| Length+1 | Residue at position Length |
| Length+2 | C-terminus (3' for RNA) |

**Formula**: For residue at zero-based index `i`, modification index = `i + 2`

```csharp
// Example: 7-residue peptide "PEPTIDE"
// Index 1 = N-terminus
// Index 2 = P (residue 1)
// Index 3 = E (residue 2)
// Index 8 = E (residue 7, the last one)
// Index 9 = C-terminus

// To get mod at zero-based position i:
if (AllModsOneIsNterminus.TryGetValue(i + 2, out Modification mod)) { ... }
```

---

## Fragmentation Classes

### Product (Theoretical Fragment)
Represents a theoretical fragment ion:
- `ProductType`: Ion type (b, y, c, z for peptides; a, b, c, d, w, x, y, z for oligos)
- `FragmentNumber`: Position in sequence where fragmentation occurs
- `NeutralMass`: Calculated mass without charge
- `ResiduePosition`: One-based position in parent sequence

### MatchedFragmentIon (Experimental Match)
Links theoretical Product to observed peak:
- `NeutralTheoreticalProduct`: The matched Product
- `Charge`: Observed charge state
- `Mz`: Observed mass-to-charge ratio
- `Intensity`: Peak intensity

---

## Equality Semantics (Critical for Parsimony)

### PeptideWithSetModifications Equality
Two peptides are equal ONLY if ALL of these match:
1. `FullSequence` (includes modifications)
2. `Parent` protein reference
3. `DigestionParams` (protease, missed cleavages, etc.)
4. `OneBasedStartResidue` and `OneBasedEndResidue`

**Important**: Identical sequences from different proteases are NOT equal:
```csharp
// These are NOT equal even with same sequence:
var peptide1 = protein.Digest(trypsinParams).First();
var peptide2 = protein.Digest(chymotrypsinParams).First();
// peptide1.Equals(peptide2) == false (different DigestionParams)
```

### Cross-Type Equality
`PeptideWithSetModifications` and `OligoWithSetMods` are NEVER equal, even with identical sequences. Type safety is enforced.

---

## Common Patterns

### Digesting a Protein
```csharp
var protein = new Protein("MPEPTIDEK", "P12345");
var digestionParams = new DigestionParams(
    protease: "trypsin",
    maxMissedCleavages: 2,
    minPeptideLength: 7,
    maxPeptideLength: 30);
var peptides = protein.Digest(digestionParams, fixedMods, variableMods);
```

### Generating Fragments
```csharp
var products = new List<Product>();
peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
```

### Matching Fragments to Spectrum
```csharp
var matchedIons = MetaMorpheusEngine.MatchFragmentIons(
    spectrum, products, commonParameters);
```

---

## Project Structure

```
mzLib/
├── mzLib.sln                   # Main solution file
│
├── # Core Domain Libraries
├── Proteomics/                 # Protein classes, digestion, modifications
├── Transcriptomics/            # Nucleic acid (RNA/DNA) classes, digestion
├── Omics/                      # Shared interfaces (IBioPolymer, fragmentation)
├── Chemistry/                  # Chemical formulas, elements, isotopes
├── MassSpectrometry/           # Spectrum handling, peak detection, deconvolution
│
├── # File I/O
├── Readers/                    # Unified file reading infrastructure
├── FileReaders/                # Specific file format readers
├── MzML/                       # mzML file format support
├── Mgf/                        # MGF file format support
├── MzIdentML/                  # mzIdentML format support
├── pepXML/                     # pepXML format support
├── ThermoRawFileReader/        # Thermo RAW file support
│
├── # Quantification
├── FlashLFQ/                   # Label-free quantification engine
├── Quantification/             # Quantification utilities
├── Chromatography/             # Chromatographic peak handling
│
├── # Analysis Tools
├── SpectralAveraging/          # Spectrum averaging algorithms
├── BayesianEstimation/         # Bayesian statistical methods
├── Predictions/                # Prediction models (retention time, etc.)
├── PredictionClients/          # Clients for external prediction services
│
├── # Utilities
├── MzLibUtil/                  # General utilities, chemistry helpers
├── UsefulProteomicsDatabases/  # Database downloading and parsing
├── mzPlot/                     # Plotting utilities
│
├── # Development/Testing
├── Test/                       # Main NUnit test project
├── TestFlashLFQ/               # FlashLFQ-specific tests
├── Benchmark/                  # Performance benchmarks
├── Development/                # Development utilities
└── Profiling/                  # Performance profiling tools
```
