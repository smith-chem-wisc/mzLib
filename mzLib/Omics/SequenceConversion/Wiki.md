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
                                    │
                    Domain bridges (`SequenceConversionExtensions`)
```

Every conversion flows through the immutable canonical intermediate representation (IR). Parsers do only enough to populate `CanonicalSequence`, serializers focus on expressing that IR for their format, and `IModificationLookup` implementations lazily enrich missing metadata. `SequenceConversionService` orchestrates the flow, while `SequenceConversionExtensions` convert `IBioPolymer`/`IBioPolymerWithSetMods` objects into (and back from) the canonical representation without callers needing to construct format strings.

## Conversion Lifecycle

0. **Domain bridge (optional)** – `SequenceConversionExtensions.ToCanonicalSequence` and `CanonicalSequenceBuilderExtensions` turn mzLib biopolymers into canonical IR without touching strings.
1. **Format detection / selection** – Callers pass a format key or defer to `SequenceConversionService.ParseAutoDetect`, which probes registered parsers via `CanParse`.
2. **Parsing** – The selected parser (usually derived from `SequenceParserBase`) tokenizes according to its `SequenceFormatSchema` and drives a `CanonicalSequenceBuilder`.
3. **Canonicalization** – The builder emits an immutable `CanonicalSequence` containing the base sequence, the source format, and a sorted list of `CanonicalModification` values.
4. **Lookup enrichment (optional)** – Serializers invoke `SequenceSerializerBase.EnrichModificationsIfNeeded`, which asks their configured `IModificationLookup` for only the modifications where `ShouldResolveMod` returns true. `ModificationLookupBase` supplies caching, tolerance handling, and fuzzy matching.
5. **Handling incompatible elements** – `SequenceConversionHandlingMode` and `ConversionWarnings` capture how failures propagate (throw, return null, strip mods, or fall back to the primary sequence). All public APIs write into `ConversionWarnings` before acting on the mode.
6. **Serialization** – The serializer rebuilds the format-specific string using its schema. Many serializers add format-specific policies (Chronologer length caps, UniProt suppression rules, Essential mod whitelists, etc.).
7. **Domain re-hydration (optional)** – `SequenceConversionExtensions.ConvertModifications` applies lookup-enriched modifications back onto `IBioPolymer`/`IBioPolymerWithSetMods` instances, honoring the same handling modes.

## Core Concepts & Their Interactions

### Canonical IR

| Type | Interaction Highlights |
|------|------------------------|
| `CanonicalSequence` | Immutable hub storing `BaseSequence`, `SourceFormat`, and an ordered array of `CanonicalModification`. Offers helpers like `GetModificationAt`, `NTerminalModification`, `HasModifications`, and `WithModifications` to encourage pure transformations. |
| `CanonicalModification` | Captures `ModificationPositionType`, residue index, target residue, original representation, optional mass / formula / IDs, and a resolved `Modification`. `WithResolvedModification` enriches UNIMOD IDs when lookups succeed. |
| `CanonicalSequenceBuilder` | Parser-owned mutable object that streams characters and mods before emitting the immutable IR. Includes validation helpers and overloads for residue / N-term / C-term modifications. |
| `CanonicalSequenceBuilderExtensions` | Bridge methods for `IBioPolymerWithSetMods` that honor the `AllModsOneIsNterminus` indexing convention, plus overloads that infer `originalRepresentation` from mzLib IDs or masses. |

### Format Abstractions

| Component | Role in the flow |
|-----------|------------------|
| `SequenceFormatSchema` | Describes bracket characters, separators, casing, and any additional knobs a format needs (e.g., Chronologer state tokens, MassShift decimal precision). Multiple parsers/serializers can share a schema. |
| `SequenceParserBase` | Supplies bracket-aware parsing, balanced bracket checks, terminal handling, and consistent error reporting via `SequenceConversionHelpers`. Format-specific subclasses only implement `ParseModificationString` and `CanParse`. |
| `SequenceSerializerBase` | Handles bracketed output, optional terminal separators, and modification enrichment. Subclasses override `GetModificationString`, `CanSerialize`, and `ShouldResolveMod` (plus `SerializeInternal` when formats do not emit bracket syntax, e.g., Chronologer). |
| `ISequenceConverter` / `SequenceConverter` | Compose a parser + serializer for the `{sourceFormat}-{targetFormat}` pipeline. Converter instances can be pre-registered or synthesized on demand. |
| `SequenceConversionService` | Registry + orchestrator located in `SequenceConversionService.cs`. Exposes `Parse`, `Serialize`, `Convert`, `ParseAutoDetect`, `ConvertAutoDetect`, format discovery APIs, and dynamic converter construction. |
| `SequenceConversionExtensions` | Extension hub (`SequenceConversionExtensions.cs`) containing helpers for serializing `IBioPolymer`/`IBioPolymerWithSetMods`, converting modifications in place, and choosing sensible defaults (e.g., mzLib schema when no source format is provided). |

## Integration Helpers

- **`SequenceConversionExtensions`** – Convert domain objects to canonical IR, serialize them with a target format key or `ISequenceSerializer`, and push lookup-enriched modifications back into peptides/proteins/nucleic acids. Handling modes map to concrete actions (clearing mods, removing incompatible entries, etc.).
- **`SequenceConversionHelpers`** – Shared error utilities for parsers/serializers so `ConversionWarnings` always capture the failure before modes are applied.
- **`ConversionWarnings`** – Aggregates warnings, incompatible elements, and the final `ConversionFailureReason`. Every public API in the service accepts an optional instance and will instantiate one when callers pass `null`.

## Modification Lookups

`IModificationLookup` sits between raw IR and serialized output. Parsers keep whatever descriptors they receive (names, partial IDs, masses). Serializers decide which modifications must be resolved and defer to a lookup. `ModificationLookupBase` centralizes caching, fuzzy matching, tolerance checks, motif filters, and identifier normalization.

| Lookup | Backing catalog | Primary usage |
|--------|-----------------|----------------|
| `MzLibModificationLookup` | mzLib / MetaMorpheus mod dictionaries | Default for mzLib parsing/serialization; preserves canonical mzLib IDs and motifs. |
| `GlobalModificationLookup` | `Mods.AllKnownMods` (proteomic + nucleic acid sources) | Blanket fallback for serializers that only require masses or general identifiers (MassShift, Chronologer, Essential, Unimod serializer enrichment). |
| `UnimodModificationLookup` | `Mods.UnimodModifications` | Resolves `UNIMOD:*` identifiers, mass-based lookups, or strings referencing the UNIMOD ontology. Prioritized by `UnimodSequenceSerializer`. |
| `UniProtModificationLookup` | `Mods.UniprotModifications` | Extensive name normalization (residue-specific variants, motif expansions) to match UniProt PTM vocabulary before serialization. |

Serializers specify which lookup they need via their constructor and may override `ShouldResolveMod` to keep enrichment cost proportional to the target format’s requirements.

## Supported Formats & Default Registrations

`SequenceConversionService.Default` registers parsers, serializers, and a handful of ready-made converters (see `SequenceConversionService.CreateDefault`). Additional converter combinations are synthesized dynamically whenever both the parser and serializer exist. Built-in formats:

| Format | Parser | Serializer | Schema | Default lookup | Notes |
|--------|--------|------------|--------|----------------|-------|
| **mzLib** | `MzLibSequenceParser` | `MzLibSequenceSerializer` | `MzLibSequenceFormatSchema` | `MzLibModificationLookup` | Native mzLib syntax for `PeptideWithSetModifications`/`OligoWithSetMods`. Serves as both primary source and target format. |
| **MassShift** | `MassShiftSequenceParser` | `MassShiftSequenceSerializer` | `MassShiftSequenceFormatSchema` | `GlobalModificationLookup` | Mass-only representation with configurable decimal precision. Used by spectral libraries and workflows that only know delta masses. |
| **Chronologer** | – (one-way) | `ChronologerSequenceSerializer` | `ChronologerSequenceFormatSchema` | `GlobalModificationLookup` | Produces fixed-length single-character encodings for ML models (Koina). Enforces amino acid whitelist and length caps. |
| **Unimod** | – | `UnimodSequenceSerializer` | `UnimodSequenceFormatSchema` | `GlobalModificationLookup` → `UnimodModificationLookup` | Emits `[UNIMOD:id]` tokens. Drops or errors on modifications that cannot be mapped to an accession. |
| **UniProt** | – | `UniProtSequenceSerializer` | `UniProtSequenceSchema` | `UniProtModificationLookup` | Applies UniProt naming rules, suppresses certain fixed mods, and replaces residues when required. Default handling mode removes incompatible elements so consumers can continue. |
| **Essential** | – | `EssentialSequenceSerializer` | `EssentialSequenceFormatSchema` | `GlobalModificationLookup` | Mirrors MetaMorpheus “Essential Sequence” exports with a pruned set of modification types (glycosylation, UniProt, biological, etc.). |

The service exposes `AvailableSourceFormats`, `AvailableTargetFormats`, and `AvailableConverters`. When a converter is not pre-registered, `GetConverter` automatically composes a new `SequenceConverter` so long as the source parser and target serializer exist.

## Service Layer Responsibilities

- Maintain parser/serializer/converter registries keyed by format name (case-insensitive) and expose the registries via read-only collections.
- Provide `RegisterParser`, `RegisterSerializer`, and `RegisterConverter` so hosts can extend the system at runtime.
- Parse (`Parse`) and serialize (`Serialize`) sequences for explicit format keys, handling null/empty inputs and unknown formats with `ConversionWarnings` + `SequenceConversionHandlingMode`.
- Convert in a single call (`Convert`), including convenience methods `ParseAutoDetect` and `ConvertAutoDetect` that iterate through registered parsers using their `CanParse` heuristics.
- Detect formats (`DetectFormat`) and expose `GetParser` / `GetSerializer` helpers for host tooling.
- Dynamically create converters on demand (cross product of available parsers + serializers) and cache them for future calls.

## Interaction With the Rest of mzLib

- **Proteomics / Transcriptomics domains** – `SequenceConversionExtensions` allow `PeptideWithSetModifications`, `Protein`, `OligoWithSetMods`, and related types to adopt canonical conversions when serializing/deserializing sequences.
- **Prediction clients (Koina)** – Each model supplies an `ISequenceConverter` (e.g., `mzLib-Chronologer` or `mzLib-Unimod`) to enforce allowed ontologies, re-label modifications, and surface warnings to users.
- **MetaMorpheus essential exports** – Use `EssentialSequenceSerializer` to emit curated modification sets consistent with search parameters.
- **File readers / writers** – Identification formats (mzIdentML, pepXML, MGF, etc.) can target whichever sequence format best matches their syntax while relying on the same lifecycle.
- **Visualization / spectral libraries** – Mass-only or projection-friendly formats (MassShift, Chronologer) provide deterministic conversions for scoring, clustering, or UI highlights.

## Extending the System

1. **Add or re-use a schema** – Derive from `SequenceFormatSchema` if your format uses new delimiters, separators, or specials. Schemas can be shared across parsers/serializers.
2. **Implement parser and/or serializer** – Subclass `SequenceParserBase` or `SequenceSerializerBase` (or implement the interfaces directly for exotic formats) and rely on `CanonicalSequenceBuilder` for IR construction.
3. **Register with the service** – Update `SequenceConversionService.Default` (or call the registration methods at runtime) so the service can discover your implementation.
4. **Optionally provide a custom lookup** – Implement `IModificationLookup` if your format depends on a specific ontology or naming convention.
5. **Decide handling policy** – Surface the preferred `SequenceConversionHandlingMode` or expose knobs so callers can specify the policy per invocation.
6. **Integrate with domain helpers** – If domain objects need extra metadata (e.g., new modification dictionaries), extend `SequenceConversionExtensions` accordingly.

## Modification Resolution Flow

```
CanonicalSequence.Modifications
        │
        ├─ serializer.ShouldResolveMod(mod)?
        │        │
        │        └─► IModificationLookup (MzLib / Global / Unimod / UniProt)
        │                  │
        │                  └─ returns enriched CanonicalModification (with masses / IDs)
        │
SequenceSerializerBase writes tokens ▸ SequenceConversionExtensions.ConvertModifications can apply the same lookup output back onto domain objects.
```

Lookups can filter their candidate sets (e.g., the Koina TMT model supplies a curated list of allowed UNIMOD IDs) to keep resolutions deterministic.

## Usage Patterns

### Parse → Serialize (Mass Shift example)

```csharp
var warnings = new ConversionWarnings();
var canonical = SequenceConversionService.Default.Parse(
    input: "[Acetyl]-PEP[Oxidation on M]TIDE",
    sourceFormat: "mzLib",
    warnings);

var massShift = SequenceConversionService.Default.Serialize(
    canonical!.Value,
    targetFormat: "MassShift",
    warnings,
    mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);
```

### Source-to-target conversion with relaxed policy

```csharp
var warnings = new ConversionWarnings();
var chronologer = SequenceConversionService.Default.Convert(
    input: "PEP[+15.995]TIDE",
    sourceFormat: "MassShift",
    targetFormat: "Chronologer",
    warnings,
    mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);
```

### Auto-detect and convert to a lookup-dependent format

```csharp
var warnings = new ConversionWarnings();
var canonical = SequenceConversionService.Default.ParseAutoDetect(sequenceString, warnings);
if (canonical is null)
{
    Console.WriteLine(string.Join("; ", warnings.Warnings));
}
else
{
    var unimod = SequenceConversionService.Default.Serialize(
        canonical.Value,
        targetFormat: "Unimod",
        warnings,
        mode: SequenceConversionHandlingMode.ReturnNull);
}
```

### Domain integration via extensions

```csharp
// Serialize directly from a peptide with set mods
var essential = peptideWithMods.Serialize(
    targetFormat: "Essential",
    warnings: new ConversionWarnings(),
    mode: SequenceConversionHandlingMode.RemoveIncompatibleElements);

// Force UniProt-compatible modifications on an IBioPolymerWithSetMods instance
peptideWithMods.ConvertModifications(UniProtSequenceSerializer.Instance);

// Convert a protein’s original modification dictionaries using a lookup
protein.ConvertModifications(GlobalModificationLookup.Instance);
```

## Project Map

```
SequenceConversion/
├── BaseClasses/                       // Shared parser/serializer/lookup bases
├── CanonicalSequence/                 // IR definitions + builder extensions
├── Implementations/
│   ├── Chronologer/                   // One-way serializer + schema
│   ├── Essential/                     // Essential sequence schema + serializer
│   ├── Global/                        // GlobalModificationLookup
│   ├── MassShift/                     // Parser, serializer, schema
│   ├── MzLib/                         // Parser, serializer, schema, lookup
│   ├── Unimod/                        // Serializer, schema, lookup
│   └── Uniprot/                       // Serializer, schema, lookup
├── Util/                              // ConversionWarnings, helpers, enums
├── Interfaces (IModificationLookup, ISequence*)
├── SequenceConversionExtensions.cs    // Domain bridges & helpers
├── SequenceConversionService.cs       // Registry + orchestration
├── SequenceConverter.cs               // Parser/serializer composition
└── Wiki.md                          // This wiki
```

## Contributors & Context

SequenceConversion underpins mzLib’s shared infrastructure for proteomics and transcriptomics. It powers Koina model input sanitation, MetaMorpheus export paths, format interop for file readers, and visualization tooling. Contributions should prefer composability: add schemas, parsers/serializers, or lookups instead of embedding special cases, and surface incompatibilities via `ConversionWarnings` so higher-level components stay informed without losing control of their error-handling strategy.