using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;

namespace Quantification;

/// <summary>
/// The peptide row for spectral matches that are ambiguous only in where their modifications sit: one base sequence,
/// the same modifications, two or more sites (oxidation on M3 or on M5). It stands for all of those sites at once and
/// claims none of them, so a match the spectrum could not localize is still quantified, and its row is never mistaken
/// for one of the localized forms.
/// </summary>
/// <remarks>
/// <see cref="FullSequence"/> is the sites' full sequences, each once, in ordinal order, joined by '|' (the separator
/// MetaMorpheus's Full Sequence column uses for an ambiguous match). Every match with the same base sequence and
/// modifications is one row, named by every site those matches leave open, so two ambiguous matches of one peptide are
/// not two peptides. A reader holding this row among other peptides never gets one site's value as if it were the
/// row's:
/// <list type="bullet">
/// <item>Members the constructor makes the same for every site (base sequence, length, mass, formula, the number of
/// modifications) come from the sites.</item>
/// <item>Members that can differ between sites (parent, residues, missed cleavages, digestion, description, the fixed
/// and variable counts) are given only when every site agrees, and otherwise throw
/// <see cref="NotSupportedException"/> naming the member; read them from each of <see cref="Sites"/> instead. They do
/// differ on real inputs: with <c>RespectCleavageBlockingModifications</c> a peptide's missed cleavages depend on where
/// its modification sits, and one site may be known only in one protein and the next only in another.</item>
/// <item><see cref="AllModsOneIsNterminus"/> holds only the modifications whose site every candidate agrees on.</item>
/// <item>A localization-ambiguous peptide has no single set of fragments, so <see cref="Fragment"/>,
/// <see cref="FragmentInternally"/> and <see cref="Localize"/> throw; fragment its <see cref="Sites"/> instead.</item>
/// </list>
/// Each site is one copy per full sequence (<c>QuantificationEngine.OnePerFullSequence</c>), so agreement is among those
/// copies: as for any peptide row, a sequence also found in another protein is not listed again.
/// </remarks>
public sealed class LocalizationAmbiguousPeptide : IBioPolymerWithSetMods
{
    /// <param name="sites">One copy per full sequence, all of one base sequence and one set of modifications; at
    /// least two.</param>
    public LocalizationAmbiguousPeptide(IEnumerable<IBioPolymerWithSetMods> sites)
    {
        Sites = sites.OrderBy(s => s.FullSequence, StringComparer.Ordinal).ToList();
        if (Sites.Count < 2 || Sites.Select(s => s.FullSequence).Distinct(StringComparer.Ordinal).Count() != Sites.Count)
        {
            throw new ArgumentException("A localization-ambiguous peptide needs two or more sites, one copy of each.", nameof(sites));
        }
        if (Sites.Select(SiteFreeIdentity).Distinct().Count() != 1 || SiteFreeIdentity(Sites[0]) == null)
        {
            throw new ArgumentException("Every site must be one base sequence carrying the same modifications.", nameof(sites));
        }

        FullSequence = string.Join("|", Sites.Select(s => s.FullSequence));
        var first = Sites[0];
        AllModsOneIsNterminus = first.AllModsOneIsNterminus
            .Where(kvp => Sites.All(s => s.AllModsOneIsNterminus.TryGetValue(kvp.Key, out var mod) && mod.Equals(kvp.Value)))
            .ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
    }

    /// <summary>The localized forms this peptide may be, one copy per full sequence, in ordinal order.</summary>
    public IReadOnlyList<IBioPolymerWithSetMods> Sites { get; }

    public string FullSequence { get; }
    public string BaseSequence => Sites[0].BaseSequence;
    public double MonoisotopicMass => Sites[0].MonoisotopicMass;
    public ChemicalFormula ThisChemicalFormula => Sites[0].ThisChemicalFormula;
    public double MostAbundantMonoisotopicMass => Sites[0].MostAbundantMonoisotopicMass;
    public string SequenceWithChemicalFormulas => string.Join("|", Sites.Select(s => s.SequenceWithChemicalFormulas));
    public int OneBasedStartResidue => Agreed(s => s.OneBasedStartResidue, nameof(OneBasedStartResidue));
    public int OneBasedEndResidue => Agreed(s => s.OneBasedEndResidue, nameof(OneBasedEndResidue));
    public int MissedCleavages => Agreed(s => s.MissedCleavages, nameof(MissedCleavages));
    public string Description => Agreed(s => s.Description, nameof(Description));
    public CleavageSpecificity CleavageSpecificityForFdrCategory
    {
        get => Agreed(s => s.CleavageSpecificityForFdrCategory, nameof(CleavageSpecificityForFdrCategory));
        set => throw new NotSupportedException("Set the cleavage specificity on the sites of a localization-ambiguous peptide.");
    }
    public char PreviousResidue => Agreed(s => s.PreviousResidue, nameof(PreviousResidue));
    public char NextResidue => Agreed(s => s.NextResidue, nameof(NextResidue));
    public IDigestionParams DigestionParams => Agreed(s => s.DigestionParams, nameof(DigestionParams));
    public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
    public int NumMods => Sites[0].NumMods;
    public int NumFixedMods => Agreed(s => s.NumFixedMods, nameof(NumFixedMods));
    public int NumVariableMods => Agreed(s => s.NumVariableMods, nameof(NumVariableMods));
    public int Length => Sites[0].Length;
    public IBioPolymer Parent => Agreed(s => s.Parent, nameof(Parent));

    /// <summary>
    /// A member's value when every site gives the same one. A member that can differ between sites is read through this,
    /// never from <see cref="Sites"/>[0], so nobody takes one site's value for the row's without knowing.
    /// </summary>
    /// <exception cref="NotSupportedException">The sites disagree, so the row has no single value.</exception>
    private T Agreed<T>(Func<IBioPolymerWithSetMods, T> member, string memberName)
    {
        T value = member(Sites[0]);
        if (Sites.Skip(1).All(site => EqualityComparer<T>.Default.Equals(member(site), value)))
        {
            return value;
        }
        throw new NotSupportedException(
            $"The sites of localization-ambiguous peptide {FullSequence} disagree on {memberName}, so it has no single value; read it from each of Sites.");
    }

    public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
        List<Product> products, IFragmentationParams? fragmentationParams = null) =>
        throw new NotSupportedException("A localization-ambiguous peptide has no single set of fragments; fragment its Sites.");

    public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
        List<Product> products, IFragmentationParams? fragmentationParams = null) =>
        throw new NotSupportedException("A localization-ambiguous peptide has no single set of fragments; fragment its Sites.");

    public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize) =>
        throw new NotSupportedException("A localization-ambiguous peptide has no single site to localize a mass on.");

    public bool Equals(IBioPolymerWithSetMods? other) =>
        other is LocalizationAmbiguousPeptide ambiguous && string.Equals(FullSequence, ambiguous.FullSequence, StringComparison.Ordinal);

    public override bool Equals(object? obj) => Equals(obj as IBioPolymerWithSetMods);

    public override int GetHashCode() => StringComparer.Ordinal.GetHashCode(FullSequence);

    public override string ToString() => FullSequence;

    /// <summary>
    /// What a peptide is with its modification sites left out: the base sequence and each modification by what it is,
    /// sorted. A modification is what it is by its type, its name without the residue
    /// (<see cref="Modification.OriginalId"/>, so Phosphorylation on S and on T are one modification) and its chemical
    /// formula (by its mass when it has none). Null when a modification has neither a formula nor a mass. Two positional
    /// isomers give the same string. A formula alone does not say what a modification is: Hydroxylation on D and
    /// Oxidation on M are both one oxygen and give different strings, as two oxidations and one dioxidation do.
    /// MetaMorpheus's <c>SpectralMatch.GetSiteFreeIdentity</c> judges a match the same way.
    /// </summary>
    public static string? SiteFreeIdentity(IBioPolymerWithSetMods peptide)
    {
        var keys = new List<string>();
        foreach (var modification in peptide.AllModsOneIsNterminus.Values)
        {
            string? composition = modification.ChemicalFormula?.Formula
                ?? modification.MonoisotopicMass?.ToString("F5", CultureInfo.InvariantCulture);
            if (composition == null)
            {
                return null;
            }
            keys.Add(modification.ModificationType + ":" + (modification.OriginalId ?? modification.IdWithMotif) + " " + composition);
        }
        keys.Sort(StringComparer.Ordinal);
        return peptide.BaseSequence + " " + string.Join(";", keys);
    }

    /// <summary>
    /// The site-free identity a spectral match is ambiguous over, when it is ambiguous only in where its modifications
    /// sit: its candidates have two or more full sequences, and one <see cref="SiteFreeIdentity"/>. Null otherwise --
    /// one full sequence, or candidates that differ in base sequence or in their modifications.
    /// </summary>
    public static string? SiteFreeIdentityOrNull(ISpectralMatch spectralMatch)
    {
        var candidates = spectralMatch.GetIdentifiedBioPolymersWithSetMods()?.Where(bp => bp != null).ToList();
        if (candidates == null || candidates.Select(c => c.FullSequence).Distinct(StringComparer.Ordinal).Take(2).Count() < 2)
        {
            return null;
        }

        string? identity = SiteFreeIdentity(candidates[0]);
        return identity != null && candidates.All(c => SiteFreeIdentity(c) == identity) ? identity : null;
    }
}
