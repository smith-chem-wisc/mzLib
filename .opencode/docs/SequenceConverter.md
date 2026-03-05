# SequenceConverter Usage Guide

## Overview
`SequenceConverter` standardizes peptide/oligo sequences between modification naming conventions and
provides mass-shift formatting for tools that cannot parse full modification names.

Supported conventions:
- MetaMorpheus
- UniProt
- Unimod

## Quick Start

```csharp
using Omics.Modifications;
using Omics.Modifications.Conversion;

var converter = SequenceConverter.Default;
string uniprot = converter.ConvertFullSequence(
    "PEPT[Common Fixed:Carbamidomethyl on C]IDE",
    ModificationNamingConvention.MetaMorpheus,
    ModificationNamingConvention.UniProt);
```

## Common Use Cases

### 1. Convert a peptide from MetaMorpheus to Unimod

```csharp
string unimod = SequenceConverter.Default.ConvertFullSequence(
    peptide.FullSequence,
    ModificationNamingConvention.MetaMorpheus,
    ModificationNamingConvention.Unimod);
```

### 2. Convert with fallback handling

```csharp
bool ok = SequenceConverter.Default.TryConvertFullSequence(
    peptide,
    ModificationNamingConvention.UniProt,
    SequenceConversionHandlingMode.RemoveIncompatibleMods,
    out var converted,
    out var reason);
```

### 3. Mass-shift notation for display

```csharp
string massShift = SequenceConverter.ToMassShiftNotation(peptide.FullSequence);
```

### 4. Parse mass-shift notation back to modifications

```csharp
var mods = SequenceConverter.FromMassShiftNotation(massShift);
```

### 5. Chronologer-ready sequences

```csharp
string chronologer = peptide.ToChronologerSequence(
    SequenceConversionHandlingMode.RemoveIncompatibleMods);
```

### 6. ProteinDbWriter normalization

```csharp
var options = new ProteinDbWriterConversionOptions
{
    Enabled = true,
    TargetConvention = ModificationNamingConvention.UniProt,
    HandlingMode = SequenceConversionHandlingMode.RemoveIncompatibleMods
};

ProteinDbWriter.WriteXmlDatabase(
    additionalMods,
    proteins,
    outputPath,
    updateTimeStamp: false,
    conversionOptions: options);
```

### 7. Koina normalization (opt-in)

```csharp
var options = new KoinaSequenceConversionOptions
{
    Enabled = true,
    TargetConvention = ModificationNamingConvention.Unimod,
    HandlingMode = SequenceConversionHandlingMode.RemoveIncompatibleMods
};

var model = new Prosit2019iRT(peptides, out var warnings, options);
```

## Naming Conventions

Examples:
- MetaMorpheus: `PEPT[Common Fixed:Carbamidomethyl on C]IDE`
- UniProt: `PEPT[UniProt:Carboxymethyl cysteine on C]IDE`
- Unimod: `PEPT[Unimod:Carbamidomethyl on C]IDE`

## Handling Modes

`SequenceConversionHandlingMode` controls what happens when a mod is not convertible:
- `ThrowException`: throw `SequenceConversionException`.
- `RemoveIncompatibleMods`: drop incompatible mods.
- `UsePrimarySequence`: return the base sequence.
- `UseMassShifts`: return mass-shift notation if possible.
- `KeepOriginalAnnotation`: keep the original annotation.
- `ReturnNull`: return null (Try* methods only).

## Performance Tips

- Reuse `SequenceConverter.Default` instead of new instances.
- Prefer `TryConvert*` for bulk pipelines to avoid exceptions.
- Use `ModificationNamingConvention.Mixed` for unknown source inputs.

## Troubleshooting

- `SequenceConversionException` with `NoEquivalent` means no compatible mod exists in the target.
- `AmbiguousEquivalent` indicates multiple candidates match; refine the input or target convention.
- If mass-shift conversion fails, the mod may have no monoisotopic mass.
