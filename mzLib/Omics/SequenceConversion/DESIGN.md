# Universal Sequence Converter - Design Documentation

## Overview

The **SequenceConversion** system provides a universal infrastructure for converting peptide and oligonucleotide sequences between different format representations. It uses a canonical intermediate representation (`CanonicalSequence`) that acts as a format-agnostic bridge, enabling conversion from any supported input format to any supported output format.

## Design Intent

### Core Goals

1. **Format-agnostic representation** - A `CanonicalSequence` serves as the universal intermediate representation (IR) that can be created from any input format and serialized to any output format

2. **Multi-format support** - Extensible parser/serializer architecture supporting:
   - **mzLib format**: `PEP[Oxidation on M]TIDE`, `[Acetyl]PEPTIDE`
   - **Mass shift format**: `PEP[+15.995]TIDE`
   - **Chronologer format**: `-PEPmTIDE_` (single-character encoding for deep learning)
   - Future formats (UNIMOD, ProForma, etc.)

3. **Modification enrichment** - `IModificationLookup` interface enables resolving partial modification information (e.g., just a mass) to fully-resolved mzLib `Modification` objects with UNIMOD IDs, chemical formulas, etc.

4. **Graceful degradation** - Multiple handling modes (`ThrowException`, `ReturnNull`, `RemoveIncompatibleElements`, `UsePrimarySequence`) for dealing with incompatible modifications

5. **Lazy enrichment** - Modifications store only what the source format provides; additional fields can be populated on-demand during serialization

## Architecture

```
┌──────────────────────────────────────────────────────────────────────────┐
│                    SequenceConversionService                              │
│  (Orchestrates parsing/serialization, manages registered formats)        │
└─────────────────────────────┬────────────────────────────────────────────┘
                              │
         ┌────────────────────┼────────────────────┐
         │                    │                    │
         ▼                    ▼                    ▼
   ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
   │  Parsers    │    │ Canonical   │    │ Serializers │
   │             │    │ Sequence    │    │             │
   │ MzLib       │───▶│ (Universal  │───▶│ MzLib       │
   │ MassShift   │    │  IR)        │    │ Chronologer │
   └─────────────┘    └─────────────┘    └─────────────┘
         │                    │                    │
         │                    │                    │
         ▼                    ▼                    ▼
   ┌─────────────┐    ┌─────────────────┐   ┌─────────────┐
   │  Schemas    │    │ Modification    │   │  Schemas    │
   │             │    │ Lookups         │   │             │
   │ MzLib       │    │                 │   │ MzLib       │
   │ MassShift   │    │ Unimod          │   │ Chronologer │
   │ Chronologer │    │                 │   │ MassShift   │
   └─────────────┘    └─────────────────┘   └─────────────┘
```

## Core Components

### Data Types

| Component | Type | Purpose |
|-----------|------|---------|
| **`CanonicalSequence`** | `readonly record struct` | Immutable universal intermediate representation with base sequence + modifications |
| **`CanonicalModification`** | `readonly record struct` | Immutable modification representation with position, mass, IDs, and resolved `Modification` reference |
| **`CanonicalSequenceBuilder`** | `class` | Mutable builder for constructing `CanonicalSequence` during parsing |
| **`ModificationPositionType`** | `enum` | Position type: `NTerminus`, `CTerminus`, `Residue` |

### Interfaces

| Interface | Purpose |
|-----------|---------|
| **`ISequenceParser`** | Parses string → `CanonicalSequence` |
| **`ISequenceSerializer`** | Serializes `CanonicalSequence` → string |
| **`IModificationLookup`** | Resolves/enriches modifications |
| **`ISequenceConversionService`** | Orchestrates conversions, manages registry |

### Base Classes

| Class | Purpose |
|-------|---------|
| **`SequenceFormatSchema`** | Abstract base class defining format syntax (brackets, separators) |
| **`ModificationLookupBase`** | Base class for modification lookups with shared resolution logic |

### Service Layer

| Component | Purpose |
|-----------|---------|
| **`SequenceConversionService`** | Main orchestrator with parser/serializer registry, auto-detection, and conversion API |

### Utility Types

| Component | Purpose |
|-----------|---------|
| **`ConversionWarnings`** | Accumulates warnings/errors during conversion |
| **`SequenceConversionException`** | Typed exception for conversion failures |
| **`SequenceConversionHandlingMode`** | Configurable behavior for incompatible elements |
| **`ConversionFailureReason`** | Enum of failure reasons |

## Implemented Formats

### mzLib Format

**Parser**: `MzLibSequenceParser`  
**Serializer**: `MzLibSequenceSerializer`  
**Schema**: `MzLibSequenceFormatSchema`

**Format characteristics**:
- Uses square brackets `[ ]` for modification annotations
- N-terminal modifications directly precede sequence with **no separator**: `[Acetyl]PEPTIDE`
- C-terminal modifications use hyphen separator: `PEPTIDE-[Amidated]`
- Modification identifiers use `IdWithMotif` format: `Oxidation on M`
- Can include modification type prefix: `Common Fixed:Carbamidomethyl on C`

**Examples**:
```
PEPTIDE                                          # Unmodified
PEP[Oxidation on M]TIDE                         # Residue modification
[Acetyl]PEPTIDE                                  # N-terminal (NO separator)
PEPTIDE-[Amidated]                              # C-terminal (with separator)
[Acetyl]PEP[Oxidation on M]TIDE-[Amidated]     # Multiple modifications
```

### Mass Shift Format

**Parser**: `MassShiftSequenceParser`  
**Serializer**: ❌ Not implemented  
**Schema**: `MassShiftSequenceFormatSchema`

**Format characteristics**:
- Uses square brackets `[ ]` for modification annotations
- Modifications specified as mass shifts with sign: `[+15.995]` or `[-18.011]`
- N-terminal modifications directly precede sequence with **no separator**: `[+42.011]PEPTIDE`
- C-terminal modifications use hyphen separator: `PEPTIDE-[+0.984]`

**Examples**:
```
PEPTIDE                                    # Unmodified
PEP[+15.995]TIDE                          # Oxidation by mass
[+42.011]PEPTIDE                          # N-terminal acetylation
PEPTIDE-[+0.984]                          # C-terminal amidation
PEP[-18.011]TIDE                          # Water loss (negative shift)
```

**Use cases**:
- Exact modification identities unknown
- Working with raw mass spectrometry data
- Converting between systems with different modification databases
- Non-standard modifications not in databases

### Chronologer Format

**Parser**: ❌ Not needed (one-way conversion)  
**Serializer**: `ChronologerSequenceSerializer`  
**Schema**: `ChronologerSequenceFormatSchema`

**Format characteristics**:
- Single-character encoding for deep learning models
- Lowercase letters represent modified residues
- Special tokens for N/C-terminus states
- Maximum sequence length: 50 amino acids (+ 2 for termini)
- Does NOT use bracket-based modification annotations

**Alphabet**:
- **Canonical amino acids**: `ACDEFGHIKLMNPQRSTVWY` (positions 1-20)
- **Modified residues**: `cmdestyabunopqrxz` (positions 21-37)
- **N/C terminus states**: `-^()&*_` (positions 38-44)
- **User-defined slots**: `0123456789` (positions 45-54)

**Modification encoding**:
| Code | Modification | Target | Mass Shift |
|------|--------------|--------|------------|
| `m` | Oxidation | M | +15.995 |
| `c` | Carbamidomethyl | C | +57.021 |
| `d` | Alternative C mod | C | +39.99 |
| `e` | PyroGlu | E/Q | -18.01/-17.02 |
| `s` | Phosphorylation | S | +79.966 |
| `t` | Phosphorylation | T | +79.966 |
| `y` | Phosphorylation | Y | +79.966 |
| `a` | Acetylation | K | +42.011 |
| `b` | Succinylation | K | +100.0 |
| `u` | Ubiquitination | K | +114.0 |
| `n` | Methylation | K | +14.016 |
| `o` | Dimethylation | K | +28.031 |
| `p` | Trimethylation | K | +42.047 |
| `q` | Methylation | R | +14.016 |
| `r` | Dimethylation | R | +28.031 |
| `z` | GlyGly | K | +224.1 |
| `x` | Heavy GlyGly | K | +229.1 |

**N-terminus states**:
| Token | Meaning |
|-------|---------|
| `-` | Free N-terminus |
| `^` | N-terminal acetylation |
| `)` | PyroGlu at N-terminus (from E) |
| `(` | Cyclized CAM-Cys at N-terminus |
| `&` | N-terminal GlyGly |
| `*` | N-terminal heavy GlyGly |

**C-terminus state**:
| Token | Meaning |
|-------|---------|
| `_` | C-terminus |

**Examples**:
```
-PEPTIDE_                  # Unmodified
-PEPmTIDE_                # Oxidized methionine
^PEPTIDE_                 # N-terminal acetylation
-PEPsTIDE_                # Phosphorylation on serine
```

## Modification Resolution

### `IModificationLookup` Interface

Enables resolving and enriching modifications by looking up additional information from modification databases.

**Resolution strategies**:
1. **UNIMOD ID**: `UNIMOD:35` → Oxidation
2. **mzLib ID**: `Oxidation on M`
3. **Mass**: Fuzzy matching within tolerance
4. **Chemical formula**: Exact formula matching
5. **Name patterns**: Flexible name matching

### Implemented Lookups

**`MzLibModificationLookup`**:
- Resolves via mzLib `Mods` class
- Searches protein and/or RNA modifications
- Supports mass-based fallback (0.001 Da tolerance)
- Singleton instances: `Instance`, `ProteinOnly`, `RnaOnly`

**`UnimodModificationLookup`**:
- Status: Implementation exists but not verified in this review
- Should resolve via UNIMOD database references

## Error Handling

### `SequenceConversionHandlingMode` Enum

Configures behavior when encountering incompatible elements:

| Mode | Behavior |
|------|----------|
| `ThrowException` | Throw `SequenceConversionException` immediately |
| `ReturnNull` | Return `null` without throwing |
| `RemoveIncompatibleElements` | Remove incompatible modifications with warnings |
| `UsePrimarySequence` | Fall back to unmodified base sequence |

### `ConversionWarnings` Class

Accumulates issues during conversion:
- **Warnings**: Non-fatal issues (e.g., using original representation)
- **Errors**: Issues that affected conversion
- **Incompatible items**: Specific modifications that couldn't be converted
- **Failure reason**: Fatal error reason if conversion failed

### `SequenceConversionException`

Typed exception with:
- `ConversionFailureReason` enum
- List of incompatible items
- List of warnings

## Usage Examples

### Basic Conversion

```csharp
// Parse mzLib format to canonical
var service = SequenceConversionService.Default;
var canonical = service.Parse("[Acetyl]PEP[Oxidation on M]TIDE", "mzLib");

// Serialize to Chronologer format
var chronologer = service.Serialize(canonical.Value, "Chronologer");
// Result: "^PEPmTIDE_"
```

### Direct Format Conversion

```csharp
var service = SequenceConversionService.Default;
var chronologer = service.Convert(
    "[Acetyl]PEPTIDE", 
    sourceFormat: "mzLib", 
    targetFormat: "Chronologer");
// Result: "^PEPTIDE_"
```

### Auto-Detection

```csharp
var service = SequenceConversionService.Default;
var detectedFormat = service.DetectFormat("PEP[+15.995]TIDE");
// Result: "MassShift"

var canonical = service.ParseAutoDetect("PEP[+15.995]TIDE");
```

### With Error Handling

```csharp
var warnings = new ConversionWarnings();
var canonical = service.Parse(
    "PEP[UnknownMod]TIDE",
    "mzLib",
    warnings,
    SequenceConversionHandlingMode.RemoveIncompatibleElements);

if (warnings.HasIncompatibleItems)
{
    Console.WriteLine($"Removed: {string.Join(", ", warnings.IncompatibleItems)}");
}
```

### Using Modification Lookup

```csharp
var lookup = MzLibModificationLookup.Instance;
var serializer = new MzLibSequenceSerializer(lookup);

// Canonical sequence with only mass information
var canonical = new CanonicalSequence(
    "PEPTIDE",
    ImmutableArray.Create(
        CanonicalModification.AtResidue(2, 'P', "+79.966", mass: 79.966)
    ),
    "MassShift");

// Serializer will use lookup to resolve mass to "Phospho on S"
var mzLibFormat = serializer.Serialize(canonical);
```

## Implementation Status

### Completed ✅

1. **Core Data Types**
   - `CanonicalSequence` with immutable modifications
   - `CanonicalModification` with all metadata fields
   - `CanonicalSequenceBuilder` for fluent construction
   - `ModificationPositionType` enum

2. **Schema Architecture**
   - Abstract `SequenceFormatSchema` base class with inheritance
   - `MzLibSequenceFormatSchema` with empty N-term separator
   - `MassShiftSequenceFormatSchema`
   - `ChronologerSequenceFormatSchema` with encoding constants

3. **Parsers**
   - `MzLibSequenceParser` - bracket annotation parsing
   - `MassShiftSequenceParser` - mass shift parsing

4. **Serializers**
   - `MzLibSequenceSerializer` - mzLib format output
   - `ChronologerSequenceSerializer` - single-char encoding with UNIMOD/mass mappings

5. **Modification Resolution**
   - `IModificationLookup` interface
   - `ModificationLookupBase` abstract base with shared logic
   - `MzLibModificationLookup` implementation

6. **Service Layer**
   - `SequenceConversionService` with full orchestration
   - Format auto-detection
   - Parser/serializer registry

7. **Error Handling**
   - `ConversionWarnings` accumulator
   - `SequenceConversionException` typed exception
   - `SequenceConversionHandlingMode` enum
   - `ConversionFailureReason` enum

8. **Build Status**: ✅ Compiles with 0 errors

### Remaining Work ❌

| Task | Priority | Notes |
|------|----------|-------|
| **Unit Tests** | **High** | Comprehensive test coverage needed for all components |
| **MassShiftSequenceSerializer** | Medium | Enable outputting `[+15.995]` format |
| **UnimodModificationLookup verification** | Medium | File exists but implementation not verified |
| **Integration with mzLib domain** | Low | Wire into `PeptideWithSetModifications.FullSequence` after testing |
| **ProForma format support** | Low | Not a current priority |
| **Documentation** | Low | XML docs exist, consider user guide |

## Design Decisions

### Why Canonical Sequence as Intermediate Representation?

- **Separation of concerns**: Parsing and serialization are independent
- **Extensibility**: New formats only require implementing one interface
- **Composability**: Can chain conversions or apply transformations
- **Testability**: Each parser/serializer can be tested in isolation
- **Performance**: Reusable intermediate representation avoids re-parsing

### Why Lazy Enrichment?

- **Flexibility**: Source formats provide different levels of detail
- **Efficiency**: Don't resolve modifications unless needed
- **Compatibility**: Can work with partial information
- **Preservation**: Original representation always maintained

### Why Schema Inheritance?

- **Type safety**: Each format has its own schema class
- **Discoverability**: Format-specific constants co-located with schema
- **Extensibility**: Easy to add format-specific properties
- **Clean API**: Singleton pattern for easy access

### Why Separate Mass Shift Format?

- **Distinct semantics**: Mass shifts represent different information than named modifications
- **Parser complexity**: Different detection and validation logic
- **User clarity**: Clear separation between formats
- **Future flexibility**: Could have different serialization rules

## Future Considerations

### Potential Format Additions

1. **ProForma** - `M[Oxidation]` or `M[U:35]` notation
2. **UNIMOD** - `(UniMod:35)` style notation  
3. **PSI-MOD** - PSI Modification Ontology format
4. **MaxQuant** - MaxQuant-specific notation
5. **Mascot** - Mascot search engine format

### Potential Enhancements

1. **Batch conversion** - Process collections of sequences efficiently
2. **Streaming API** - Handle large files without loading all into memory
3. **Format validation** - Validate sequences against format rules
4. **Format migration** - Tools for migrating between database formats
5. **Ambiguity resolution** - Handle cases where modifications are ambiguous
6. **Localization scoring** - Preserve modification localization confidence

### Integration Points

1. **`PeptideWithSetModifications.FullSequence`** - Use for reading/writing
2. **`OligoWithSetMods`** - RNA/DNA sequence support
3. **File readers** - Parse sequences from identification files
4. **Search engines** - Input/output format conversions
5. **Visualization** - Display sequences in user-preferred formats

## Project Structure

```
mzLib/Omics/SequenceConversion/
├── CanonicalSequence.cs               # Universal IR
├── CanonicalModification.cs           # Modification IR
├── CanonicalSequenceBuilder.cs        # Builder for parsing
├── SequenceFormatSchema.cs            # Abstract schema base
├── SequenceConversionService.cs       # Main orchestrator
├── ISequenceParser.cs                 # Parser interface
├── ISequenceSerializer.cs             # Serializer interface
├── ISequenceConversionService.cs      # Service interface
├── IModificationLookup.cs             # Lookup interface
│
├── Schema/
│   ├── MzLibSequenceFormatSchema.cs   # mzLib format definition
│   ├── MassShiftSequenceFormatSchema.cs
│   └── ChronologerSequenceFormatSchema.cs
│
├── Parsers/
│   ├── MzLibSequenceParser.cs
│   └── MassShiftSequenceParser.cs
│
├── Serializers/
│   ├── MzLibSequenceSerializer.cs
│   └── ChronologerSequenceSerializer.cs
│
├── ModificationLookup/
│   ├── ModificationLookupBase.cs      # Base class
│   ├── MzLibModificationLookup.cs
│   └── UnimodModificationLookup.cs
│
├── Util/
│   ├── ConversionWarnings.cs
│   ├── SequenceConversionException.cs
│   ├── SequenceConversionHandlingMode.cs
│   └── ConversionFailureReason.cs
│
└── Converters/                        # Empty - future use
```

## Contributors

This design was developed as part of the mzLib project to provide universal sequence format conversion capabilities for mass spectrometry proteomics and transcriptomics analysis.

## Version History

- **Current**: Initial implementation with mzLib, MassShift, and Chronologer format support
- **Future**: Testing, integration with domain model, and additional format support
