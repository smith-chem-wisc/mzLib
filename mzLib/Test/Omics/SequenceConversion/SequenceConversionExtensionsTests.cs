using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Omics.SequenceConversion;
using System;
using System.Collections.Generic;
using System.Linq;
using Test.Omics;

namespace Test.Omics.SequenceConversion;

[TestFixture]
    public class SequenceConversionExtensionsTests
    {
    [Test]
    public void ToCanonicalSequence_WithSetMods_HonorsSpecialIndices()
    {
        var mods = new Dictionary<int, Modification>
        {
            { 1, CreateModification("Nterm", "P", "N-terminal.") },
            { 4, CreateModification("Residue", "P") },
            { "PEPTIDE".Length + 2, CreateModification("Cterm", "E", "C-terminal.") }
        };

        var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE", mods: mods);

        var canonical = peptide.ToCanonicalSequence();

        Assert.That(canonical.NTerminalModification, Is.Not.Null);
        Assert.That(canonical.CTerminalModification, Is.Not.Null);
        Assert.That(canonical.GetModificationAt(2), Is.Not.Null);
    }

    [Test]
    public void ToCanonicalSequence_BioPolymerFallsBackToLocalizedMods()
    {
        var protein = new MockBioPolymer("PEPTIDE", "P12345");
        protein.OneBasedPossibleLocalizedModifications[3] = new List<Modification>
        {
            CreateModification("Residue", "P")
        };

        var warnings = new ConversionWarnings();
        var canonical = protein.ToCanonicalSequence(warnings);

        Assert.That(canonical.HasModifications, Is.True);
        Assert.That(warnings.HasWarnings, Is.True);
        Assert.That(warnings.Warnings.Single(), Does.Contain("Falling back"));
    }

    [Test]
    public void ConvertModifications_WithUsePrimarySequenceClearsDictionary()
    {
        var mods = new Dictionary<int, Modification>
        {
            { 2, CreateModification("Residue", "P") }
        };
        var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE", mods: mods);
        var lookup = new RecordingLookup(true);
        var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.UsePrimarySequence);

        peptide.ConvertModifications(serializer);

        Assert.That(peptide.AllModsOneIsNterminus, Is.Empty);
    }

    [Test]
    public void ConvertModifications_RemovesUnresolvedEntriesInRemoveMode()
    {
        var mods = new Dictionary<int, Modification>
        {
            { 2, CreateModification("Residue", "P") }
        };
        var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE", mods: mods);
        var lookup = new RecordingLookup(resolveSuccessfully: false);
        var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        peptide.ConvertModifications(serializer);

        Assert.That(peptide.AllModsOneIsNterminus, Is.Empty);
    }

    [Test]
        public void ConvertModifications_ReplacesEntriesWhenLookupResolves()
        {
            var originalMod = CreateModification("Residue", "P");
            var mods = new Dictionary<int, Modification> { { 2, originalMod } };
            var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE", mods: mods);
            var lookup = new RecordingLookup(resolveSuccessfully: true);
            var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.RemoveIncompatibleElements);

            peptide.ConvertModifications(serializer);

            Assert.That(peptide.AllModsOneIsNterminus[2], Is.Not.SameAs(originalMod));
            Assert.That(peptide.AllModsOneIsNterminus[2].ModificationType, Is.EqualTo("Resolved"));
        }

        [Test]
        public void ToCanonicalSequence_BioPolymerPrefersOriginalModsAndWarnsOnDuplicates()
        {
            var protein = new MockBioPolymer("PEPT", "P12345");
            protein.OriginalNonVariantModifications[1] = new List<Modification>
            {
                CreateModification("NTerm", "P", "N-terminal.")
            };
            protein.OriginalNonVariantModifications[protein.BaseSequence.Length] = new List<Modification>
            {
                CreateModification("CTerm", "T", "C-terminal.")
            };
            protein.OriginalNonVariantModifications[2] = new List<Modification>
            {
                CreateModification("Primary", "E"),
                CreateModification("Secondary", "E")
            };

            var warnings = new ConversionWarnings();
            var canonical = protein.ToCanonicalSequence(warnings, "unit-test");

            Assert.That(canonical.NTerminalModification, Is.Not.Null);
            Assert.That(canonical.NTerminalModification!.Value.OriginalRepresentation, Does.Contain("NTerm"));
            Assert.That(canonical.CTerminalModification, Is.Not.Null);
            Assert.That(canonical.CTerminalModification!.Value.OriginalRepresentation, Does.Contain("CTerm"));
            Assert.That(canonical.GetModificationAt(1)?.OriginalRepresentation, Does.Contain("Primary"));
            Assert.That(warnings.HasWarnings, Is.True);
            Assert.That(warnings.Warnings.Any(w => w.Contains("Multiple modifications")), Is.True);
        }

        [Test]
        public void Serialize_BioPolymerUsesProvidedServiceAndForwardsMode()
        {
            var protein = new MockBioPolymer("PEPT", "P12345");
            var serializer = new RecordingSequenceSerializer();
            var service = new SequenceConversionService();
            service.RegisterSerializer(serializer);

            var warnings = new ConversionWarnings();
            var result = protein.Serialize(serializer.FormatName, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements, service);

            Assert.That(serializer.Calls, Has.Count.EqualTo(1));
            Assert.That(serializer.Calls.Single().BaseSequence, Is.EqualTo("PEPT"));
            Assert.That(serializer.RecordedWarnings.Single(), Is.SameAs(warnings));
            Assert.That(serializer.RecordedModes.Single(), Is.EqualTo(SequenceConversionHandlingMode.RemoveIncompatibleElements));
            Assert.That(result, Is.EqualTo("recording:PEPT:RemoveIncompatibleElements"));
        }

        [Test]
        public void Serialize_WithSetModsUsesSerializerDirectly()
        {
            var peptide = new MockBioPolymerWithSetMods("PEPT", "PEPT");
            var serializer = new RecordingSequenceSerializer();
            var warnings = new ConversionWarnings();

            var result = peptide.Serialize(serializer, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

            Assert.That(serializer.Calls, Has.Count.EqualTo(1));
            Assert.That(serializer.Calls.Single().BaseSequence, Is.EqualTo("PEPT"));
            Assert.That(serializer.RecordedWarnings.Single(), Is.SameAs(warnings));
            Assert.That(serializer.RecordedModes.Single(), Is.EqualTo(SequenceConversionHandlingMode.RemoveIncompatibleElements));
            Assert.That(result, Is.EqualTo("recording:PEPT:RemoveIncompatibleElements"));
        }

        [Test]
        public void ConvertModifications_BioPolymerUsePrimarySequenceClearsAllCollections()
        {
            var protein = new MockBioPolymer("PEPTIDE", "P12345");
            protein.OneBasedPossibleLocalizedModifications[2] = new List<Modification> { CreateModification("Possible", "E") };
            protein.OriginalNonVariantModifications[3] = new List<Modification> { CreateModification("Original", "P") };

            var seqVar = new SequenceVariation(1, 1, "P", "A", "seq", new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { CreateModification("Variant", "P") } }
            });
            var appliedVar = new SequenceVariation(1, 1, "P", "A", "applied", new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { CreateModification("Applied", "P") } }
            });

            protein.SequenceVariations.Add(seqVar);
            protein.AppliedSequenceVariations.Add(appliedVar);

            var lookup = new RecordingLookup(true);
            var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.UsePrimarySequence);

            protein.ConvertModifications(serializer);

            Assert.That(protein.OneBasedPossibleLocalizedModifications, Is.Empty);
            Assert.That(protein.OriginalNonVariantModifications, Is.Empty);
            Assert.That(seqVar.OneBasedModifications, Is.Empty);
            Assert.That(appliedVar.OneBasedModifications, Is.Empty);
        }

        [Test]
        public void ConvertModifications_BioPolymerRemoveModeDropsFailedLookups()
        {
            var protein = new MockBioPolymer("PEPTIDE", "P12345");
            protein.OriginalNonVariantModifications[2] = new List<Modification>
            {
                CreateModification("First", "E"),
                CreateModification("Second", "E")
            };

            var lookup = new ScriptedLookup(new[] { false, true });
            var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.RemoveIncompatibleElements);

            protein.ConvertModifications(serializer);

            Assert.That(protein.OriginalNonVariantModifications[2], Has.Count.EqualTo(1));
            Assert.That(protein.OriginalNonVariantModifications[2][0].ModificationType, Is.EqualTo("Resolved"));
        }

        [Test]
        public void ConvertModifications_BioPolymerThrowModeThrowsWhenLookupFails()
        {
            var protein = new MockBioPolymer("PEPTIDE", "P12345");
            protein.OriginalNonVariantModifications[2] = new List<Modification> { CreateModification("Fail", "E") };

            var lookup = new RecordingLookup(false);
            var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.ThrowException);

            Assert.That(() => protein.ConvertModifications(serializer), Throws.TypeOf<SequenceConversionException>());
        }

        [Test]
        public void ConvertModifications_BioPolymerResolvesVariantPositions()
        {
            var protein = new MockBioPolymer("PEPTIDE", "P12345");
            protein.SequenceVariations.Add(new SequenceVariation(1, 1, "P", "A", "seq", new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { CreateModification("VariantMod", "A") } }
            }));

            var lookup = new RecordingLookup(true);
            var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.RemoveIncompatibleElements);

            protein.ConvertModifications(serializer);

            var variantMods = protein.SequenceVariations.Single().OneBasedModifications[1];
            Assert.That(variantMods.Single().ModificationType, Is.EqualTo("Resolved"));
        }

    [Test]
    public void ConvertModifications_WithSetModsThrowModeThrowsWhenLookupFails()
    {
        var mods = new Dictionary<int, Modification> { { 2, CreateModification("Fail", "P") } };
            var peptide = new MockBioPolymerWithSetMods("PEPTIDE", "PEPTIDE", mods: mods);
            var lookup = new RecordingLookup(false);
            var serializer = new LookupBackedSerializer(lookup, SequenceConversionHandlingMode.ThrowException);

        Assert.That(() => peptide.ConvertModifications(serializer), Throws.TypeOf<SequenceConversionException>());
    }

    private static Modification CreateModification(string id, string residue, string location = "Anywhere.")
    {
        if (!ModificationMotif.TryGetMotif(residue, out var motif))
        {
            ModificationMotif.TryGetMotif("A", out motif);
        }
        return new Modification(id, null, "Test", null, motif, location, null, 10);
    }

    private sealed class RecordingLookup : IModificationLookup
    {
        private readonly bool _resolveSuccessfully;

        public RecordingLookup(bool resolveSuccessfully)
        {
            _resolveSuccessfully = resolveSuccessfully;
        }

        public string Name => "recording";

        public CanonicalModification? TryResolve(CanonicalModification mod)
        {
            if (!_resolveSuccessfully)
            {
                return null;
            }

            if (!ModificationMotif.TryGetMotif(mod.TargetResidue?.ToString() ?? "X", out var motif))
            {
                ModificationMotif.TryGetMotif("A", out motif);
            }
            var resolved = new Modification("Resolved", null, "Resolved", null, motif, "Anywhere.", null, 15.99);
            return mod.WithResolvedModification(resolved, mod.ResidueIndex, mod.PositionType);
        }

        public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, Chemistry.ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
            => throw new System.NotImplementedException();
    }

    private sealed class LookupBackedSerializer : ISequenceSerializer
    {
        public LookupBackedSerializer(IModificationLookup lookup, SequenceConversionHandlingMode mode)
        {
            ModificationLookup = lookup;
            HandlingMode = mode;
        }

        public string FormatName => "lookup";
        public SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;
        public IModificationLookup? ModificationLookup { get; }
        public SequenceConversionHandlingMode HandlingMode { get; }
        public bool CanSerialize(CanonicalSequence sequence) => true;
        public bool ShouldResolveMod(CanonicalModification mod) => false;
        public string? Serialize(CanonicalSequence sequence, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => sequence.BaseSequence;
        public Dictionary<int, Modification> ToOneIsNterminusModificationDictionary(
            CanonicalSequence sequence,
            Dictionary<string, Modification>? knownMods = null,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            if (mode == SequenceConversionHandlingMode.UsePrimarySequence)
            {
                return new Dictionary<int, Modification>();
            }

            var projected = new Dictionary<int, Modification>();
            foreach (var modification in sequence.Modifications)
            {
                var resolved = ModificationLookup?.TryResolve(modification);
                if (!resolved.HasValue || resolved.Value.MzLibModification == null)
                {
                    if (mode == SequenceConversionHandlingMode.ThrowException)
                    {
                        throw new SequenceConversionException(
                            $"Unable to resolve modification {modification}",
                            ConversionFailureReason.IncompatibleModifications,
                            new[] { modification.ToString() });
                    }

                    continue;
                }

                var index = modification.PositionType switch
                {
                    ModificationPositionType.NTerminus => 1,
                    ModificationPositionType.CTerminus => sequence.BaseSequence.Length + 2,
                    _ => modification.ResidueIndex.GetValueOrDefault() + 2
                };

                projected[index] = resolved.Value.MzLibModification;
            }

            return projected;
        }
    }

    private sealed class RecordingSequenceSerializer : ISequenceSerializer
    {
        public RecordingSequenceSerializer(string formatName = "recording", SequenceConversionHandlingMode handlingMode = SequenceConversionHandlingMode.ThrowException)
        {
            FormatName = formatName;
            HandlingMode = handlingMode;
        }

        public string FormatName { get; }
        public SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;
        public IModificationLookup? ModificationLookup => null;
        public SequenceConversionHandlingMode HandlingMode { get; }
        public List<CanonicalSequence> Calls { get; } = new();
        public List<ConversionWarnings?> RecordedWarnings { get; } = new();
        public List<SequenceConversionHandlingMode> RecordedModes { get; } = new();
        public bool CanSerialize(CanonicalSequence sequence) => true;
        public bool ShouldResolveMod(CanonicalModification mod) => false;

        public string? Serialize(CanonicalSequence sequence, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            Calls.Add(sequence);
            RecordedWarnings.Add(warnings);
            RecordedModes.Add(mode);
            return $"{FormatName}:{sequence.BaseSequence}:{mode}";
        }

        public Dictionary<int, Modification> ToOneIsNterminusModificationDictionary(
            CanonicalSequence sequence,
            Dictionary<string, Modification>? knownMods = null,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => new();
    }

    private sealed class StubProjectionParser : ISequenceParser
    {
        private readonly Func<string, bool> _canParse;
        private readonly CanonicalSequence _result;

        public StubProjectionParser(string formatName, Func<string, bool> canParse, CanonicalSequence result)
        {
            FormatName = formatName;
            Schema = MzLibSequenceFormatSchema.Instance;
            _canParse = canParse;
            _result = result;
        }

        public string FormatName { get; }

        public SequenceFormatSchema Schema { get; }

        public bool CanParse(string input) => _canParse(input);

        public CanonicalSequence? Parse(
            string input,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => _result;
    }

    private sealed class StubProjectionSerializer : ISequenceSerializer
    {
        public StubProjectionSerializer(string formatName)
        {
            FormatName = formatName;
        }

        public string FormatName { get; }

        public SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;

        public IModificationLookup? ModificationLookup => null;

        public SequenceConversionHandlingMode HandlingMode => SequenceConversionHandlingMode.ThrowException;

        public bool CanSerialize(CanonicalSequence sequence) => true;

        public bool ShouldResolveMod(CanonicalModification mod) => false;

        public string? Serialize(
            CanonicalSequence sequence,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException) => sequence.BaseSequence;

        public Dictionary<int, Modification> ToOneIsNterminusModificationDictionary(
            CanonicalSequence sequence,
            Dictionary<string, Modification>? knownMods = null,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            if (knownMods == null)
            {
                return new Dictionary<int, Modification>();
            }

            return new Dictionary<int, Modification>
            {
                { 42, knownMods["Carbamidomethyl on C"] }
            };
        }
    }

    private sealed class ScriptedLookup : IModificationLookup
    {
        private readonly Queue<bool> _outcomes;

        public ScriptedLookup(IEnumerable<bool> outcomes)
        {
            _outcomes = new Queue<bool>(outcomes);
        }

        public string Name => "scripted";

        public CanonicalModification? TryResolve(CanonicalModification mod)
        {
            var shouldResolve = _outcomes.Count == 0 || _outcomes.Dequeue();
            if (!shouldResolve)
            {
                return null;
            }

            if (!ModificationMotif.TryGetMotif(mod.TargetResidue?.ToString() ?? "X", out var motif))
            {
                ModificationMotif.TryGetMotif("A", out motif);
            }

            var resolved = new Modification("Resolved", null, "Resolved", null, motif, "Anywhere.", null, 15.99);
            return mod.WithResolvedModification(resolved, mod.ResidueIndex, mod.PositionType);
        }

        public CanonicalModification? TryResolve(string originalRepresentation, char? targetResidue = null, Chemistry.ChemicalFormula? chemicalFormula = null, ModificationPositionType? positionType = null)
            => throw new NotImplementedException();
    }
}
