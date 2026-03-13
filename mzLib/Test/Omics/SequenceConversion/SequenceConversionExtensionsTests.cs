using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
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
    }
}
