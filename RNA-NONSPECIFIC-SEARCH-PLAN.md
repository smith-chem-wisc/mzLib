# RNA Non-Specific Search Plan

## Motivation

Recent digestion changes made protein non-specific searching explicit through
`SearchModeType`, `FragmentationTerminus`, `SpecificProtease`, and the
`singleN`/`singleC` digestion agents. MetaMorpheus now relies on this behavior
when it clones digestion parameters and searches from each terminus.

RNA digestion exposes the same general concepts, but the implementation is not
currently symmetric:

- `RnaDigestionParams.SearchModeType` is not used to select the effective RNase.
- `RnaDigestionParams.Clone` does not preserve all search-mode semantics.
- RNA has no equivalent of `SpecificProtease` for retaining the user-selected
  RNase when the effective agent changes to `singleN` or `singleC`.
- `OligoWithSetMods` does not resolve `CleavageSpecificity.Unknown` in the same
  way as `PeptideWithSetModifications`.
- MetaMorpheus non-specific search still contains protein-specific assumptions.

The goal is to support RNA non-specific searching while unifying the common
protein/RNA workflow. The implementation should not create a second parallel
API or redesign the search engine. Existing concrete classes and downstream
solutions should continue to be the primary entry points.

## Cross-Solution Gotchas

### Public API and source compatibility

`mzLib` is consumed by MetaMorpheus and other solutions. Changes to
`IDigestionParams`, `DigestionAgent`, constructors, or serialization behavior
can cause failures outside the mzLib repository.

Before changing a public member:

- Search all solutions for implementations of `IDigestionParams`.
- Search all uses of `DigestionAgent`, `Protease`, and `Rnase`.
- Search TOML serializers and deserializers for digestion-parameter property
  names.
- Check whether any consumers construct `DigestionParams` or
  `RnaDigestionParams` positionally.
- Check binary serialization and cached index assumptions.

Prefer optional constructor parameters and additive properties. If an existing
property is generalized, retain strongly typed aliases such as
`SpecificProtease` where downstream code already uses them.

### `DigestionAgent` is shared infrastructure

`DigestionAgent` is the common parent of proteases and RNases. Adding an
abstract member to it could break custom agents in downstream solutions.

Prefer one of these approaches, in order:

1. Add a virtual member with a safe existing default where possible.
2. Add a shared method whose input is the existing digestion-product context.
3. Only make the member abstract if all in-repository and downstream agent
   implementations can be updated deliberately.

Do not introduce a second unrelated `SpecificProtease`/`SpecificRnase`
workflow. The common concept should be represented by the existing
`DigestionAgent` abstraction.

### Serialization and cloning

MetaMorpheus clones digestion parameters for separate terminus passes. A clone
must preserve:

- the named digestion agent;
- `SearchModeType`;
- missed-cleavage and length limits;
- modification limits;
- RNA or protein-specific settings; and
- the new terminus.

The effective agent may change during cloning, but the specific agent must not.
For example, cloning a non-specific trypsin search from N to C should change
`singleN` to `singleC` while retaining trypsin as the specific agent.

The same invariant must hold for RNase T1 and the RNA `singleN`/`singleC`
agents.

### Terminology is not interchangeable

Protein termini use N/C terminology. RNA uses five-prime/three-prime
terminology. They represent analogous search directions but not identical
chemical behavior.

Shared parameter logic may unify direction selection, but product creation,
terminus formulas, and fragment-ion types must remain domain-specific.

Do not replace RNA `FivePrime`/`ThreePrime` values with protein N/C values.

### `None` does not mean “return every non-specific peptide”

The current protein contract uses `SearchModeType.None` with N/C termini to
return anchored seeds. MetaMorpheus decides the unresolved end later from the
precursor mass.

RNA must follow the same seed contract. A change that makes
`NucleicAcid.Digest` emit every possible RNA subsequence would create a large
index and would not match the existing non-specific search algorithm.

### Cleavage specificity and FDR

MetaMorpheus uses `CleavageSpecificityForFdrCategory` to classify matches. A
trimmed oligo created with `Unknown` specificity can reach FDR classification
and fail because the existing RNA path does not resolve it.

The resolution workflow should be moved to the shared digestion-product
parent. The actual site calculation remains in the concrete digestion agent.
This keeps the lifecycle common without pretending that protein cleavage and
RNA cleavage use the same chemistry.

### Fragmentation types

The non-specific engine currently contains peptide-oriented assumptions when
selecting terminus-specific product types. RNA uses five-prime/three-prime
product collections and different chemical formulas.

The search engine must select the product-type collection from the digestion
parameters or analyte type. It must not assume that a protein product-type
dictionary can handle RNA terminus values.

### Terminal modifications

Protein terminal modifications and RNA five-prime/three-prime termini have
different indexing and mass behavior. Shared code should coordinate the
operation, but the concrete peptide and oligo classes must continue to own
their terminal chemistry and modification placement.

### Existing RNA refusal behavior

MetaMorpheus currently refuses RNA non-specific searches deliberately. That
guard should remain until mzLib and the engine are ready together. Removing it
early would replace a clear refusal with runtime failures involving casts,
product types, or FDR classification.

## Phase 1 Inventory: Completed

The initial contract inventory was performed before changing any public types.

### `IDigestionParams` implementations

Production mzLib contains exactly two implementations:

- `mzLib/Proteomics/ProteolyticDigestion/DigestionParams.cs`
- `mzLib/Transcriptomics/Digestion/RnaDigestionParams.cs`

MetaMorpheus contains two test-only implementations:

- `MetaMorpheus/Test/DigestionAgentNameTests.cs`:
  `AgentlessDigestionParams`
- `MetaMorpheus/Test/ParameterTest.cs`:
  `BadDigestionParams`

Adding `SpecificDigestionAgent` to `IDigestionParams` therefore has a small,
known implementation surface, but both downstream test doubles must be updated
at the same time. Their purpose is specifically to exercise unknown or invalid
parameter implementations, so they should not receive fake protein/RNA
behavior. They should return `null` for the new property where appropriate.

The interface is consumed broadly as a parameter type by digestion, decoy
generation, quantification, spectra matching, parsimony, file-specific
parameters, and post-search analysis. The new member must therefore be
read-only and semantically stable; it must not replace `DigestionAgent`.

### Production digestion agents

The production hierarchy currently has:

- `Omics/Digestion/DigestionAgent.cs`
- `Proteomics/ProteolyticDigestion/Protease.cs`
- `Transcriptomics/Digestion/Rnase.cs`

The base `DigestionAgent` already contains shared cleavage-site logic through
`GetDigestionSiteIndices` and `GetCleavageSpecificity`. This is an important
reuse point. RNA-specific behavior should extend or correct this existing
operation rather than create a second specificity classifier.

There is also one mzLib test subclass,
`Test/Omics/Modifications/TestDigestionMotif.cs:649`, which derives directly
from `DigestionAgent`. Any new abstract agent member would require this test
agent to implement it. A virtual/default implementation is safer unless the
new behavior cannot have a meaningful base default.

### Parameter construction and cloning

Protein construction and cloning are implemented in
`Proteomics/ProteolyticDigestion/DigestionParams.cs`:

- `SpecificProtease` is assigned before the effective `Protease` is changed.
- `SearchModeType.None` selects `singleN` or `singleC`.
- `Clone` reconstructs the object and preserves the named protease.

RNA construction and cloning are implemented in
`Transcriptomics/Digestion/RnaDigestionParams.cs`:

- `Rnase` is currently the only agent property.
- `SearchModeType` defaults to `Full` but is not a constructor argument.
- `Clone` reconstructs from `Rnase.Name` and currently omits
  `SearchModeType`.
- RNA currently does not retain a separate named agent when an effective
  single-terminus agent is selected.

The shared base implementation must preserve the protein behavior exactly and
make RNA follow the same named-agent/effective-agent invariant.

### Serialization and settings paths

The mzLib source does not contain the primary TOML parameter configuration.
MetaMorpheus owns that behavior in `TaskLayer/MetaMorpheusTask.cs`:

- `ConfigureType<IDigestionParams>` selects `DigestionParams` or
  `RnaDigestionParams` during TOML deserialization.
- `SetAllFileSpecificCommonParams` reconstructs either concrete type manually.
- digestion parameters are written to `DigestionParameters.toml` for indexes.

The manual RNA reconstruction currently passes the effective
`DigestionAgent.Name`, not a retained specific RNase. This must be reviewed
when `SpecificRnase` and shared agent naming are introduced.

The MetaMorpheus GUI also constructs `RnaDigestionParams` in
`GuiFunctions/ViewModels/FragmentationParamsViewModel.cs:259-260`, and the
search-task GUI constructs digestion parameters in
`GUI/TaskWindows/SearchTaskWindow.xaml.cs`. Both paths must be checked when the
new search-mode constructor argument is added.

### Clone consumers

The important downstream clone paths are:

- `EngineLayer/CommonParameters.cs:283`, where file-specific parameters clone
  digestion parameters for a terminus.
- `TaskLayer/SearchTask/SearchTask.cs`, where semi-specific searches may run
  multiple terminus passes.
- `Test/DigestionAgentNameTests.cs`, which explicitly verifies that a cloned
  non-specific protein search reports its original protease.

The new RNA behavior must add equivalent coverage without changing the
existing protein assertions.

### Equality and hash-code dependencies

Both concrete parameter classes implement value equality and hash codes. They
are used in `HashSet<IDigestionParams>` collections throughout MetaMorpheus,
including parsimony and post-search analysis.

The specific agent must participate consistently in equality and hashing. The
effective agent alone is insufficient because two non-specific searches using
different named RNases or proteases must not compare equal merely because both
use `singleN`.

Changing equality fields can affect persisted or in-memory parameter grouping,
so existing protein equality tests and RNA equality tests must be expanded
before changing the implementation.

### Indexing dependencies

`MetaMorpheus/EngineLayer/Indexing/IndexingEngine.cs` currently detects
single-terminus indexing behavior through:

```csharp
CommonParameters.DigestionParams.DigestionAgent.Name.Contains("single")
```

This is intentionally based on the effective agent because it controls how
seed products are indexed. It should not be changed to
`SpecificDigestionAgent` without separately reviewing precursor and terminal
modification indexing.

The indexing results are typed as `IBioPolymerWithSetMods`, but several helper
methods and product-type lookups remain peptide-oriented. These are integration
issues for a later phase, not reasons to alter the digestion-parameter
contract during the inventory phase.

### Inventory conclusion

The safest implementation boundary is now clear:

1. Extend the existing `IDigestionParams` contract additively.
2. Add a shared parameter implementation base for the two production
   parameter classes.
3. Keep `DigestionAgent` as the common parent and prefer virtual/default
   additions over new abstract members.
4. Update the two MetaMorpheus test doubles and the one mzLib test agent in the
   same change if required by the selected contract shape.
5. Treat TOML reconstruction, equality/hash-code behavior, and indexing as
   explicit compatibility gates before proceeding to RNA digestion changes.

## Target Design

### Extend the existing digestion-parameter contract

Add a common property to `IDigestionParams`:

```csharp
DigestionAgent SpecificDigestionAgent { get; }
```

This is an extension of the existing contract, not a new parallel interface.
It represents the agent selected by the user, while `DigestionAgent` remains
the effective agent used to generate digestion products.

Concrete implementations retain typed properties where useful:

```csharp
DigestionParams.SpecificProtease
RnaDigestionParams.SpecificRnase
```

Both typed properties should refer to the same underlying specific-agent
concept.

### Add a small shared parameter base class

Introduce an implementation base class for the two existing parameter types,
for example `DigestionParamsBase`. It should contain only logic common to both
domains:

- common parameter storage;
- `SearchModeType` and `FragmentationTerminus` handling;
- specific-agent storage;
- effective-agent selection flow;
- clone invariants; and
- common equality/hash-code participation where safe.

The existing concrete classes remain the public domain-specific types:

```text
DigestionParams      : DigestionParamsBase
RnaDigestionParams   : DigestionParamsBase
```

The base class should use narrow protected extension points rather than
protein/RNA conditionals:

```csharp
protected abstract DigestionAgent GetSingleTerminusAgent(
    FragmentationTerminus terminus);

protected abstract IDigestionParams CreateClone(
    FragmentationTerminus terminus);
```

The base class owns the algorithm for preserving the named agent and selecting
the effective single-terminus agent. Each concrete parameter class resolves
its own dictionary entries.

### Generalize specificity resolution through existing parents

Move the common `Unknown` resolution workflow into `DigestionProduct`:

1. Detect `CleavageSpecificity.Unknown`.
2. Obtain `SpecificDigestionAgent` from the digestion parameters.
3. Ask the agent to classify the product.
4. Store the resulting FDR category and description.

`Protease` and `Rnase` provide the domain-specific classification logic.

This avoids duplicating the current `PeptideWithSetModifications` logic in
`OligoWithSetMods`. It also ensures future digestion products follow the same
classification lifecycle.

Preserve existing typed `Protease.GetCleavageSpecificity` behavior through a
compatibility overload or adapter if necessary. The protein implementation
must continue honoring initiator-methionine behavior.

### RNA parameter behavior

`RnaDigestionParams` should:

- accept `searchModeType` as an optional constructor argument;
- retain the selected RNase in `SpecificRnase`;
- select effective `singleN` or `singleC` for `SearchModeType.None`;
- preserve `SearchModeType` in `Clone`; and
- preserve all existing RNA-specific defaults and serialization behavior.

`NucleicAcid.Digest` should then use the effective RNase without changing the
existing full-digestion path.

The initial implementation should focus on `None`, because that is the
behavior required by MetaMorpheus non-specific search. Semi-specific RNA
behavior should be added only if the current search workflows require it; the
shared parameter design should not force an unimplemented RNA semi-specific
algorithm.

## Implementation Phases

### Phase 1: Contract inventory

- Identify every `IDigestionParams` implementation across the repository and
  downstream solutions.
- Identify all concrete `DigestionAgent` implementations.
- Identify serialization, TOML, equality, and index-cache dependencies.
- Record current constructor and clone behavior for protein and RNA.

### Phase 2: Shared parameter semantics (Completed)

- Added `SpecificDigestionAgent` to `IDigestionParams`.
- Added the small `DigestionParamsBase` implementation layer.
- Moved common effective-agent selection state and selection flow into the base.
- Adapted `DigestionParams` without changing protein search behavior.
- Adapted `RnaDigestionParams` with `SpecificRnase`, search-mode construction,
  and clone preservation.
- Updated the known MetaMorpheus test-only `IDigestionParams` implementations.
- Added RNA regression coverage for effective-agent selection and cloning.

The base intentionally does not replace the existing concrete scalar parameter
properties or public concrete types. This keeps serialization and downstream
construction stable while establishing the shared named-agent/effective-agent
contract needed by later phases.

### Phase 3: Shared product classification

- Move `Unknown` specificity resolution into `DigestionProduct`.
- Add the shared agent-facing specificity operation.
- Implement RNase cleavage-specificity classification.
- Preserve protein initiator-methionine behavior.
- Verify peptide construction and deserialization behavior.

### Phase 4: RNA digestion support

- Make RNA digestion use the effective RNase selected by its parameters.
- Verify singleN and singleC seed generation.
- Verify five-prime and three-prime termini and mass formulas.
- Add RNA clone and specificity regression tests.

### Phase 5: MetaMorpheus integration

- Pass `SearchModeType` when constructing `RnaDigestionParams`.
- Replace peptide-only product-type lookup in non-specific search.
- Ensure trimmed products use the shared specificity workflow.
- Update digestion-agent naming to use `SpecificDigestionAgent`.
- Remove RNA refusal guards only after the complete path works.
- Update RNA GUI parameter selection to use `RnaseDictionary`.

### Phase 6: Cross-solution verification

- Build mzLib and all affected downstream solutions.
- Run existing protein digestion and non-specific search tests.
- Run existing RNA digestion and modern-search tests.
- Add an RNA non-specific integration test in MetaMorpheus.
- Verify TOML round-tripping and cloned file-specific parameters.
- Verify no cached peptide index is accidentally reused for oligos.

## Required Tests

### mzLib

- `RnaDigestionParams` preserves `SearchModeType` through cloning.
- Non-specific RNA parameters retain `SpecificRnase` after cloning.
- Effective RNA agent changes from `singleN` to `singleC` when the terminus
  changes.
- Existing full RNase digestion output is unchanged.
- RNA singleN and singleC seeds respect minimum and maximum lengths.
- Oligos constructed with `Unknown` specificity are classified correctly.
- Existing protein `Unknown` specificity behavior is unchanged.
- Existing equality and hash-code tests continue to pass.

### MetaMorpheus

- RNA non-specific search is accepted after the guard is removed.
- The index contains `OligoWithSetMods` products.
- Five-prime and three-prime product types are selected correctly.
- Trimmed oligos receive valid cleavage-specificity categories.
- The named RNase is preserved in output and digestion-agent comparisons.
- Protein non-specific search results remain unchanged.
- RNA search parameters round-trip through TOML.

## Completion Criteria

The work is complete when RNA non-specific searches use the same parameter,
clone, effective-agent, and specificity-resolution workflow as protein searches,
while retaining RNA-specific sequence rules, termini, formulas, and product
types.

The change should be considered unsuccessful if it requires MetaMorpheus to
maintain separate protein and RNA non-specific search algorithms or if it
breaks existing consumers that only use the current digestion parameter and
agent APIs.
