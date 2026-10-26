using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Omics;
using Proteomics;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// SETTEST1's member resolution (QuantProject <c>DEF-DIFF-SET-MEMBER</c>, up to <c>shared_group</c>): which protein group
/// stands for a member, by UniProt primary gene name or frozen accession, inside the stratum's row set; the all-members
/// exception; and the outcome-blind choice rule. The H1 and H3 members are resolved against the primary names their entries
/// carry in PXD057407's own search database (<c>ReferenceData/settest_genes_PXD057407.tsv</c>).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SetTestMembersTests
{
    private static readonly string[] Samples = { "adult_1", "adult_2", "old_1", "old_2" };

    private static SetGeneIndex RealGenes()
    {
        var lines = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
            "ReferenceData", "settest_genes_PXD057407.tsv")).Skip(1).Where(l => l.Length > 0);
        return SetGeneIndex.FromPrimaryNames(lines.Select(l => l.Split('\t'))
            .Select(c => (c[0], (IEnumerable<string>)c[1].Split(';', StringSplitOptions.RemoveEmptyEntries))));
    }

    private static SetGeneIndex Genes(params (string Accession, string Genes)[] entries) =>
        SetGeneIndex.FromPrimaryNames(entries.Select(e =>
            (e.Accession, (IEnumerable<string>)e.Genes.Split(';', StringSplitOptions.RemoveEmptyEntries))));

    /// <summary>A protein group with its peptides: each peptide's values in the four samples (NaN = none).</summary>
    private sealed record Group(string Name, double QValue = 0.001, bool Contaminant = false, bool Decoy = false)
    {
        public List<(string Sequence, string Groups, double[] Values)> Peptides { get; } = new();

        public Group With(string sequence, params double[] values) => WithShared(sequence, Name, values);

        public Group WithShared(string sequence, string groups, params double[] values)
        {
            Peptides.Add((sequence, groups, values));
            return this;
        }
    }

    private static ObservationTable Table(params Group[] groups)
    {
        var runs = Samples.Select(s => new ObservationRun(s + ".raw", s)).ToList();
        var peptides = groups.SelectMany(g => g.Peptides)
            .GroupBy(p => p.Sequence).Select(p => p.First())
            .Select(p => new PeptideRunValues(p.Sequence, p.Sequence, p.Groups.Split(';'),
                Samples.Select((s, i) => (s, p.Values[i])).Where(x => !double.IsNaN(x.Item2))
                    .ToDictionary(x => x.s + ".raw", x => new RunValue(x.Item2, ObservationState.Quantified))))
            .ToList();
        var infos = groups.Select(g => new ProteinGroupInfo(g.Name, g.QValue, g.Decoy, g.Contaminant)).ToList();
        return ObservationTable.Build(runs, runs.Select(r => r.FileName), peptides, infos, QuantBasis.MsmsOnly);
    }

    private static readonly double[] All = { 1e6, 2e6, 3e6, 4e6 };
    private static double[] Only(params int[] samples) =>
        Enumerable.Range(0, 4).Select(i => samples.Contains(i) ? 1e6 : double.NaN).ToArray();
    private static readonly double[] None = { double.NaN, double.NaN, double.NaN, double.NaN };

    private static SetDatasetContrast Contrast(params string[] members)
    {
        string header = "hypothesis_id\tkind\tprediction\tmember\tspecies\tdataset\tfactor\ttest_level\treference_level\tcontexts\tmrna_log2\tindependence\tpost_hoc";
        var rows = members.Select(m => $"H9\tset\tup\t{m}\tmouse\tPXD1\tage group\told\tadult\t\t\tindependent\tno");
        return SetTestFreeze.Read(Encoding.UTF8.GetBytes(string.Join('\n', new[] { header }.Concat(rows)) + "\n"), "t")
            .Parts.Single().Contrasts.Single();
    }

    private static IReadOnlyList<SetMemberAssignment> Resolve(SetDatasetContrast contrast, ObservationTable table,
        SetGeneIndex genes, IReadOnlyList<string>? partMembers = null) =>
        SetMemberResolver.Resolve(contrast, partMembers ?? contrast.Members, table, Samples, genes);

    private static SetMemberAssignment One(IReadOnlyList<SetMemberAssignment> all, string member) =>
        all.Single(a => a.Member == member);

    [Test]
    public void TheRealH3Members_TakeEachStateUpToSharedGroup()
    {
        // H3 up on PXD057407, from the real freeze; the groups are synthetic, named by the members' real accessions.
        var h3 = SetTestFreeze.Read(SetTestFreezeTests.RealExcerpt(), "hypotheses-v4").Parts[1];
        var table = Table(
            new Group("P04441").With("CD74PEP", All),                                    // Cd74: its own group
            new Group("P01921").With("AB1PEP", All),                                     // H2-Ab1: its own group
            new Group("P01921|P04230").With("ABEBPEP", All),                             // H2-Ab1 + H2-Eb1, both members
            new Group("P01910", QValue: 0.05).With("AAPEP", All),                        // H2-Aa: above 1%
            new Group("P28063", Contaminant: true).With("PSMB8PEP", All));               // Psmb8: a contaminant
        // Psmb9 (P28076) is in the database but in no group.

        var resolved = Resolve(h3.Contrasts.Single(), table, RealGenes(), h3.Members);

        Assert.That(resolved.Select(a => a.Member), Is.EqualTo(h3.Members));
        Assert.That(One(resolved, "Cd74").State, Is.Null);
        Assert.That(One(resolved, "Cd74").ProteinGroup, Is.EqualTo("P04441"));
        Assert.That(One(resolved, "Cd74").GroupGenes, Is.EqualTo(new[] { "Cd74" }));
        Assert.That(One(resolved, "Cd74").MemberAccessionsUsed, Is.EqualTo(new[] { "P04441" }));
        Assert.That(One(resolved, "H2-Ab1").ProteinGroup, Is.EqualTo("P01921"));
        // H2-Ab1 has a group of its own, so the H2-Ab1/H2-Eb1 group does not enter (ruled 2026-10-25): H2-Eb1 is shared.
        Assert.That(One(resolved, "H2-Eb1").State, Is.EqualTo(SetMemberState.SharedGroup));
        Assert.That(One(resolved, "H2-Aa").State, Is.EqualTo("not_in_rowset:above_1pct_q"));
        Assert.That(One(resolved, "Psmb8").State, Is.EqualTo("not_in_rowset:contaminant"));
        Assert.That(One(resolved, "Psmb9").State, Is.EqualTo(SetMemberState.NoGroup));
    }

    [Test]
    public void AMemberNotInTheDatabases_IsUnmapped_AndAMemberMatchesInAnyCase()
    {
        var table = Table(new Group("P04441").With("CD74PEP", All));
        var resolved = Resolve(Contrast("Hsd3b5", "CD74"), table, RealGenes());
        Assert.That(One(resolved, "Hsd3b5").State, Is.EqualTo(SetMemberState.Unmapped));
        Assert.That(One(resolved, "CD74").ProteinGroup, Is.EqualTo("P04441"));
    }

    [Test]
    public void AGeneLessAccession_CountsAsAnotherGene_UnlessTheMemberIsFrozenByAccession()
    {
        // Igha's reviewed entry P01878 has no gene name in PXD057407's database (AGING-D8).
        var table = Table(
            new Group("A0A0A6YXW6|P01878").With("IGHAPEP", All),
            new Group("P01878").With("IGHACONST", All));

        var byName = Resolve(Contrast("Igha"), table, RealGenes());
        Assert.That(One(byName, "Igha").State, Is.EqualTo(SetMemberState.SharedGroup));

        var excerpt = SetTestFreeze.Read(SetTestFreezeTests.RealExcerpt(), "hypotheses-v4");
        var pxd010091 = excerpt.Parts[0].Contrasts.Single(c => c.Dataset == "PXD010091");
        var byAccession = SetMemberResolver.Resolve(pxd010091, excerpt.Parts[0].Members, table, Samples, RealGenes());
        Assert.That(One(byAccession, "Igha").ProteinGroup, Is.EqualTo("P01878"));
        Assert.That(One(byAccession, "Igha").MemberAccessionsUsed, Is.EqualTo(new[] { "P01878" }));
        Assert.That(One(byAccession, "Igha").OtherGroups, Is.Empty);
    }

    [Test]
    public void AGroupPassingTheFilterWithNoValueInTheStratum_IsOutsideTheRowSet()
    {
        var table = Table(new Group("Q9D8I1").With("MZB1PEP", None));
        Assert.That(One(Resolve(Contrast("Mzb1"), table, RealGenes()), "Mzb1").State,
            Is.EqualTo("not_in_rowset:no_value_in_stratum"));
    }

    [Test]
    public void TheReasonGiven_IsTheGroupThatCameClosestToTheRowSet()
    {
        var genes = Genes(("A1", "GeneA"), ("A2", "GeneA"), ("A3", "GeneA"));
        var contaminantAndAbove = Table(new Group("A1", Contaminant: true).With("P1", All),
            new Group("A2", QValue: 0.2).With("P2", All));
        Assert.That(One(Resolve(Contrast("GeneA"), contaminantAndAbove, genes), "GeneA").State,
            Is.EqualTo("not_in_rowset:above_1pct_q"));

        var aboveAndEmpty = Table(new Group("A2", QValue: 0.2).With("P2", All), new Group("A3").With("P3", None));
        Assert.That(One(Resolve(Contrast("GeneA"), aboveAndEmpty, genes), "GeneA").State,
            Is.EqualTo("not_in_rowset:no_value_in_stratum"));
    }

    [Test]
    public void ADecoyGroup_NeverStandsForAMember()
    {
        var genes = Genes(("A1", "GeneA"));
        var table = Table(new Group("A1", Decoy: true).With("P1", All));
        Assert.That(One(Resolve(Contrast("GeneA"), table, genes), "GeneA").State, Is.EqualTo(SetMemberState.NoGroup));
    }

    [Test]
    public void AMultiGeneEntry_StandsForNoSingleMember_ButEntersOnceWhenAllItsGenesAreMembers()
    {
        // One histone H4 entry, encoded by two genes.
        var genes = Genes(("H4X", "H4c1;H4c2"), ("B1", "GeneB"));
        var table = Table(new Group("H4X").With("H4PEP", All), new Group("B1").With("BPEP", All));

        var alone = Resolve(Contrast("H4c1", "GeneB"), table, genes);
        Assert.That(One(alone, "H4c1").State, Is.EqualTo(SetMemberState.SharedGroup));

        var both = Resolve(Contrast("H4c1", "H4c2", "GeneB"), table, genes);
        Assert.That(One(both, "H4c1").ProteinGroup, Is.EqualTo("H4X"));
        Assert.That(One(both, "H4c2").ProteinGroup, Is.EqualTo("H4X"));
        Assert.That(One(both, "H4c1").GroupGenes, Is.EqualTo(new[] { "H4c1", "H4c2" }));
        Assert.That(One(both, "H4c1").ByAllMembersException, Is.True);
        Assert.That(One(both, "GeneB").ByAllMembersException, Is.False);
    }

    [Test]
    public void TheExceptionReadsThePartsMembers_NotOnlyTheContrasts()
    {
        var genes = Genes(("H4X", "H4c1;H4c2"));
        var table = Table(new Group("H4X").With("H4PEP", All));
        // The contrast lists only H4c1, but H4c2 is a member of the part in another dataset.
        var resolved = Resolve(Contrast("H4c1"), table, genes, partMembers: new[] { "H4c1", "H4c2" });
        Assert.That(One(resolved, "H4c1").ProteinGroup, Is.EqualTo("H4X"));
    }

    [Test]
    public void TheChoice_TakesTheMostSamplesWithAUniquePeptideValue()
    {
        var genes = Genes(("A1", "GeneA"), ("A2", "GeneA"));
        var table = Table(
            new Group("A1").With("A1P1", Only(0, 1)).With("A1P2", Only(0)).With("A1P3", Only(1)),
            new Group("A2").With("A2P1", Only(0, 1, 2)));
        var a = One(Resolve(Contrast("GeneA"), table, genes), "GeneA");
        Assert.That(a.ProteinGroup, Is.EqualTo("A2"));
        Assert.That(a.OtherGroups, Is.EqualTo(new[] { "A1" }));
    }

    [Test]
    public void TheChoice_ThenTakesTheMostUniquePeptides_CountingOnlyUniqueOnes()
    {
        var genes = Genes(("A1", "GeneA"), ("A2", "GeneA"), ("C1", "GeneC"));
        var table = Table(
            new Group("A1").With("A1P1", All).With("A1P2", All),
            new Group("A2").With("A2P1", All).WithShared("SHARED1", "A2;C1", All).WithShared("SHARED2", "A2;C1", All),
            new Group("C1").With("C1P1", All));
        Assert.That(One(Resolve(Contrast("GeneA"), table, genes), "GeneA").ProteinGroup, Is.EqualTo("A1"));
    }

    [Test]
    public void TheChoice_CountsOnlyValuesInTheContrastsSamples()
    {
        var genes = Genes(("A1", "GeneA"), ("A2", "GeneA"));
        var table = Table(
            new Group("A1").With("A1P1", Only(0, 1, 2, 3)),
            new Group("A2").With("A2P1", Only(0, 1)).With("A2P2", Only(0)));
        // Restricted to adult_1 and adult_2, both have two samples, and A2 has more unique peptides there.
        var a = SetMemberResolver.Resolve(Contrast("GeneA"), new[] { "GeneA" }, table, new[] { "adult_1", "adult_2" }, genes)
            .Single();
        Assert.That(a.ProteinGroup, Is.EqualTo("A2"));
    }

    [Test]
    public void TheChoice_LastBreaksATieByOrdinalName()
    {
        var genes = Genes(("A1", "GeneA"), ("a1", "GeneA"));
        var table = Table(new Group("a1").With("LOWER", All), new Group("A1").With("UPPER", All));
        Assert.That(One(Resolve(Contrast("GeneA"), table, genes), "GeneA").ProteinGroup, Is.EqualTo("A1"));
    }

    [Test]
    public void AGroupWithNoUniquePeptide_CanStillBeChosen_AndAGroupOutsideTheRowSetIsNeverOffered()
    {
        var genes = Genes(("A1", "GeneA"), ("A2", "GeneA"), ("C1", "GeneC"));
        var table = Table(
            new Group("A1").WithShared("SH", "A1;C1", All),
            new Group("A2", QValue: 0.5).With("A2P", All),
            new Group("C1").With("C1P", All));
        var a = One(Resolve(Contrast("GeneA"), table, genes), "GeneA");
        Assert.That(a.ProteinGroup, Is.EqualTo("A1"));
        Assert.That(a.OtherGroups, Is.Empty);
    }

    [Test]
    public void FromBioPolymers_ReadsOnlyPrimaryNames_OfEveryGeneElement()
    {
        var histone = new Protein("PEPTIDE", "H4X", geneNames: new List<Tuple<string, string>>
        {
            new("primary", "H4c1"), new("synonym", "Hist1h4a"), new("primary", "H4c2"),
        });
        var none = new Protein("PEPTIDEK", "P01878");
        var index = SetGeneIndex.FromBioPolymers(new IBioPolymer[] { histone, none });

        Assert.That(index.GenesOf("H4X"), Is.EqualTo(new[] { "H4c1", "H4c2" }));
        Assert.That(index.AccessionsOf("hist1h4a"), Is.Empty);
        Assert.That(index.AccessionsOf("h4C2"), Is.EqualTo(new[] { "H4X" }));
        Assert.That(index.Contains("P01878"), Is.True);
        Assert.That(index.GenesOf("P01878"), Is.Empty);
    }

    [Test]
    public void ASampleTheStratumDoesNotHave_IsRefused() =>
        Assert.That(() => SetMemberResolver.Resolve(Contrast("GeneA"), new[] { "GeneA" },
                Table(new Group("A1").With("P", All)), new[] { "elsewhere" }, Genes(("A1", "GeneA"))),
            Throws.ArgumentException.With.Message.Contains("'elsewhere'"));
}
