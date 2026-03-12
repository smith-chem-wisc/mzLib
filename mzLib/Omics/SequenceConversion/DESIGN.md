# Universal Sequence Converter – Interaction Wiki

## System at a Glance

```
┌───────────────────────────── SequenceConversionService ─────────────────────────────┐
│  Detects formats ▸ chooses parser ▸ builds CanonicalSequence ▸ serializes ▸ reports │
└──────────────┬────────────────────────────┬────────────────────────────┬────────────┘
               │                            │                            │
        Parsers & Schemas          Canonical IR & Builder        Serializers & Schemas
      (format-specific input)      (shared data contracts)       (format-specific output)
               │                            │                            │
               └───────────► Modification Lookups ◄────────────┘
                                (lazy enrichment)
```

Every conversion flows through a canonical intermediate representation (IR). Parsers only worry about “How do I get to IR?” while serializers focus on “How do I express IR for my format?”. Modification lookups and handling policies participate along the way to enrich or prune data without breaking format responsibilities.

## Conversion Lifecycle

1. **Format detection / selection** – `SequenceConversionService` auto-detects or respects the caller-provided format key.
2. **Parsing** – The format’s parser consumes the raw string and drives a `CanonicalSequenceBuilder` using its paired `SequenceFormatSchema` for lexical rules.
3. **Canonicalization** – The builder emits an immutable `CanonicalSequence` containing the base sequence and a list of `CanonicalModification` values (with whatever metadata the source could provide).
4. **Lookup enrichment (optional)** – Serializers or downstream code supply an `IModificationLookup` (e.g., mzLib or UNIMOD) to resolve missing IDs, formulas, or motifs.
5. **Handling incompatible elements** – `SequenceConversionHandlingMode` and `ConversionWarnings` coordinate whether we throw, return null, strip modifications, or fall back to primaries.
6. **Serialization** – The requested serializer, guided by its schema, recreates the output string (including any lazy-enriched metadata).

## Core Concepts & Their Interactions

### Canonical IR

| Type | Interaction Highlights |
|------|------------------------|
| `CanonicalSequence` | Immutable hub shared by all formats. Holds `BaseSequence`, the originating `FormatKey`, and an immutable array of `CanonicalModification` entries. Serializers read it; parsers write it. |
| `CanonicalModification` | Carries residue index, position type (`ModificationPositionType`), masses, IDs, and an optional resolved `Modification`. Lookups may inject enriched data before serialization. |
| `CanonicalSequenceBuilder` | A parser-owned mutable object. Ensures we can stream characters without repeated allocations, then produce the immutable IR when parsing completes. |

### Format Abstractions

| Component | Role in the flow |
|-----------|------------------|
| `SequenceFormatSchema` | Describes tokens, separators, and character casing rules. Parsers reference it while reading; serializers reference it when writing. Keeping schemas separate allows multiple parsers/serializers to share the same grammar. |
| `ISequenceParser` / `ISequenceSerializer` | Concrete implementations live under `Parsers/` and `Serializers/`. Both operate strictly on `CanonicalSequence` to stay composable. |
| `ISequenceConverter` | Composes a parser + serializer into a single pipeline. Converter names follow the convention `{sourceFormat}-{targetFormat}` (e.g., `mzLib-UNIMOD`). |

### Modification Lookups

`IModificationLookup` sits between raw IR and serialized output. Parsers deliberately keep whatever information they received (names, masses, partial IDs). When a serializer needs a canonical identifier (e.g., `[UNIMOD:35]`), it asks the lookup to resolve the `CanonicalModification`. Because lookups inherit from `ModificationLookupBase`, they all share caching, fuzzy matching, and tolerance handling.

Available implementations:

- **`MzLibModificationLookup`** – Backs onto the mzLib `Mods` catalog (protein/RNA variant singleton configurations).
- **`UnimodModificationLookup`** – Accepts a candidate set of `Modification` objects (often filtered) and resolves against UNIMOD accessions or references.

### Handling Modes & Warnings

`SequenceConversionHandlingMode` governs failure policy for both parsers and serializers:

| Mode | Interaction effect |
|------|--------------------|
| `ThrowException` | Bail immediately with `SequenceConversionException`. Used by workflows that require strict fidelity. |
| `ReturnNull` | Propagate `null` and record the reason inside `ConversionWarnings`. Upstream code can skip the sequence gracefully. |
| `RemoveIncompatibleElements` | Strip offending modifications/residues, emit warnings, and continue. Essential for prediction services that can tolerate partial information. |
| `UsePrimarySequence` | Drop all modifications and keep only the base sequence. Primarily used by tools that explicitly work on unmodified strings. |

`ConversionWarnings` travels alongside conversions and aggregates textual warnings, incompatible element summaries, and the eventual `ConversionFailureReason` (if any). Consumers such as the Koina prediction clients surface these warnings back to the user interface.

## Supported Formats & How They Interplay

| Format | Parser | Serializer | Schema | Typical Consumers |
|--------|--------|------------|--------|-------------------|
| **mzLib** | `MzLibSequenceParser` | `MzLibSequenceSerializer` | `MzLibSequenceFormatSchema` | Native mzLib APIs (`PeptideWithSetModifications`, readers, writers). |
| **Mass Shift** | `MassShiftSequenceParser` | _pending_ | `MassShiftSequenceFormatSchema` | Spectral libraries or pipelines that only know delta masses. |
| **Chronologer** | _n/a (one-way)_ | `ChronologerSequenceSerializer` | `ChronologerSequenceFormatSchema` | Machine-learning inputs (single-character encodings). |

Adding a format means supplying a schema plus at least one parser or serializer; the rest of the stack (service, warnings, lookup plumbing) immediately works. The service can auto-create converters for any registered parser/serializer pair.

## Service Layer Responsibilities

`SequenceConversionService` exposes the public API:

- `Parse(formatKey)` / `ParseAutoDetect` – produce `CanonicalSequence?`
- `Serialize(formatKey)` – write strings from IR
- `Convert(sourceFormat, targetFormat)` – convenience wrapper for parse + serialize
- Converter registry – converter keys follow `{sourceFormat}-{targetFormat}`. The service builds converters on demand whenever a parser and serializer exist for the requested pair.

Internally, the service also instantiates `ConversionWarnings`, copies handling modes into parser/serializer calls, and orchestrates lookup injection so consumers only need to provide a lookup once.

## Interaction With the Rest of mzLib

- **Proteomics/Transcriptomics domain** – `PeptideWithSetModifications` and `OligoWithSetMods` can adopt canonical conversions when reading/writing full sequences, ensuring consistent format handling.
- **Prediction clients (Koina)** – Koina models build on this infrastructure by providing an `ISequenceConverter` per model: parse mzLib strings, restrict allowed UNIMOD IDs via custom lookups, serialize to `[UNIMOD:*]`, and hand sequences to remote models while still exposing the original inputs.
- **File readers/writers** – Identification formats such as mzIdentML, pepXML, or MGF can choose whichever format best matches their syntax and rely on the same conversion lifecycle.
- **Visualization / Spectral libraries** – Tools needing mass-shift only representations can parse once and serialize into multiple downstream formats with guaranteed consistency.

## Extending the System

1. **Add a schema (if needed)** – Define bracket characters, separators, allowed alphabets, and any format-specific tokens in a `SequenceFormatSchema` derivative.
2. **Implement parser and/or serializer** – Follow `ISequenceParser`/`ISequenceSerializer`, using the schema for tokenization. Emit `CanonicalSequence` objects via the builder and consume them via immutable accessors.
3. **Register with the service** – Update the format registry (usually inside `SequenceConversionService.Default` creation) so the system knows how to reach your implementation.
4. **Optional: custom lookup** – If a format references a specific ontology, implement `IModificationLookup` that encapsulates its resolution logic.
5. **Decide handling policy** – Surface the preferred `SequenceConversionHandlingMode` to callers, or allow them to override.

## Modification Resolution Flow

```
CanonicalSequence.Modifications ──► Serializer needs identifier
         │                               │
         │                               ▼
         └───── optional lookup.Resolve(CanonicalModification)
                                 │
                    adds UNIMOD ID / formulas / mzLib IDs
                                 │
                          Serializer emits value
```

Lookups can filter their candidate sets. For example, the Koina TMT model provides a curated list of allowed UNIMOD IDs to `UnimodModificationLookup`, preventing accidental enrichment of unsupported modifications.

## Usage Patterns

### Parse → Serialize

```csharp
var service = SequenceConversionService.Default;
var canonical = service.Parse("[Acetyl]PEP[Oxidation on M]TIDE", "mzLib");
var warnings = new ConversionWarnings();
var chronologer = service.Serialize(
    canonical!.Value,
    "Chronologer",
    warnings,
    SequenceConversionHandlingMode.ThrowException,
    lookup: MzLibModificationLookup.Instance);
```

### Source-to-target Conversion with Relaxed Policy

```csharp
var warnings = new ConversionWarnings();
var converted = SequenceConversionService.Default.Convert(
    input: "PEP[UnknownMod]TIDE",
    sourceFormat: "mzLib",
    targetFormat: "MassShift",
    warnings,
    handling: SequenceConversionHandlingMode.RemoveIncompatibleElements);

if (warnings.HasIncompatibleItems)
{
    Console.WriteLine(string.Join("; ", warnings.IncompatibleItems));
}
```

### Auto-detect and Integrate With External Pipelines

```csharp
var canonical = SequenceConversionService.Default.ParseAutoDetect(sequenceString);
if (canonical is null)
{
    // consult warnings for diagnostic text
    return;
}

// Downstream pipeline can inspect canonical.Modifications regardless of source syntax
```

## Implementation Status & Roadmap

| Area | Status |
|------|--------|
| Core IR, builder, schemas | ✅ complete |
| mzLib parser/serializer | ✅ complete |
| Mass Shift parser | ✅ complete |
| Mass Shift serializer | ⏳ planned |
| Chronologer serializer | ✅ complete (one-way) |
| Lookup infrastructure | ✅ (mzLib + UNIMOD) |
| Service orchestration | ✅ |
| Unit tests | ⏳ expand coverage |
| Additional formats (ProForma, PSI-MOD, etc.) | 🚧 backlog |

## Project Map

```
SequenceConversion/
├── CanonicalSequence(.cs)               // IR definitions
├── Parsers/                             // Format-specific readers
│   ├── MzLibSequenceParser.cs
│   └── MassShiftSequenceParser.cs
├── Serializers/
│   ├── MzLibSequenceSerializer.cs
│   └── ChronologerSequenceSerializer.cs
├── Schema/
│   ├── MzLibSequenceFormatSchema.cs
│   ├── MassShiftSequenceFormatSchema.cs
│   └── ChronologerSequenceFormatSchema.cs
├── ModificationLookup/
│   ├── ModificationLookupBase.cs
│   ├── MzLibModificationLookup.cs
│   └── UnimodModificationLookup.cs
├── Util/
│   ├── ConversionWarnings.cs
│   ├── SequenceConversionHandlingMode.cs
│   ├── ConversionFailureReason.cs
│   └── SequenceConversionException.cs
└── SequenceConversionService.cs         // Registry + orchestration
```

## Contributors & Context

SequenceConversion is part of mzLib’s shared infrastructure for proteomics and transcriptomics. It powers everything from Koina model input sanitation to file reader interoperability. Contributions should favor composability: add formats or lookups rather than embedding special cases, and rely on `ConversionWarnings` + handling modes to keep higher-level components (search engines, prediction clients, visualization tools) informed without breaking their control flow.

### Version Notes

- **Current** – mzLib, MassShift (parse), Chronologer (serialize) available; Koina models rely on these pathways.
- **Planned** – MassShift serializer, ProForma support, deeper unit coverage, and tighter integration with domain objects once adoption stabilizes.
