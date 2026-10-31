using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using NUnit.Framework;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// The frozen <c>tests.tsv</c> reader of SETTEST1 (QuantProject <c>DEF-DIFF-SET-INPUT</c>): columns by header name, the
/// value vocabularies, and the flags read per member row and folded per dataset contrast (<c>DEF-DIFF-SET-POOL</c>,
/// <c>-POSTHOC</c>, <c>-MATURATION</c>). The real rows are an excerpt of aging's hypotheses-v4 freeze
/// (<c>ReferenceData/PROVENANCE_settest.txt</c>).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SetTestFreezeTests
{
    private const string Header =
        "hypothesis_id\tkind\tprediction\tmember\tspecies\tdataset\tfactor\ttest_level\treference_level\tcontexts\tmrna_log2\tindependence\tpost_hoc";

    internal static byte[] RealExcerpt() =>
        File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential",
            "ReferenceData", "settest_hypotheses_v4_excerpt.tsv"));

    private static string Row(string hypothesis = "H9", string kind = "set", string prediction = "up", string member = "GeneA",
        string dataset = "PXD1", string contexts = "", string mrna = "", string independence = "independent",
        string postHoc = "no") =>
        string.Join('\t', hypothesis, kind, prediction, member, "mouse", dataset, "age group", "old", "young", contexts, mrna,
            independence, postHoc);

    private static SetTestFreeze Read(string header, params string[] rows) =>
        SetTestFreeze.Read(Encoding.UTF8.GetBytes(string.Join('\n', new[] { header }.Concat(rows)) + "\n"), "test-freeze");

    [Test]
    public void TheRealExcerpt_ReadsEveryRowIntoItsParts()
    {
        var freeze = SetTestFreeze.Read(RealExcerpt(), "hypotheses-v4");

        Assert.That(freeze.Rows, Has.Count.EqualTo(16));
        Assert.That(freeze.Parts.Select(p => $"{p.HypothesisId} {SetTestFreeze.PredictionName(p.Prediction)}"),
            Is.EqualTo(new[] { "H1 up", "H3 up", "H5 unchanged", "H21 down" }));
        var h1 = freeze.Parts[0];
        Assert.That(h1.Members, Is.EqualTo(new[] { "Jchain", "Ighm", "Igkc", "Igha", "Mzb1" }));
        Assert.That(h1.IsSet, Is.True);
        Assert.That(h1.Contrasts.Select(c => c.Dataset), Is.EqualTo(new[] { "PXD057407", "PXD031782", "PXD010091" }));
        Assert.That(freeze.Parts[1].Members, Is.EqualTo(new[] { "Cd74", "H2-Ab1", "H2-Eb1", "H2-Aa", "Psmb8", "Psmb9" }));
        // A label kind with one member is tested as a single (AGING-D9).
        Assert.That(freeze.Parts[2].Kind, Is.EqualTo("concordance"));
        Assert.That(freeze.Parts[2].IsSet, Is.False);
        Assert.That(freeze.Rows[0].Line, Is.EqualTo(2));
    }

    [Test]
    public void TheRealExcerpt_FoldsEachContrastsFlagsOverItsMemberRows()
    {
        var freeze = SetTestFreeze.Read(RealExcerpt(), "hypotheses-v4");
        var h1 = freeze.Parts[0].Contrasts;

        // PXD057407: Ighm and Igkc reproduce [114], so the contrast is a reproduction for H1 although three rows are independent.
        Assert.That(h1[0].Members, Is.EqualTo(new[] { "Jchain", "Ighm", "Igkc", "Igha", "Mzb1" }));
        Assert.That(h1[0].Independence, Is.EqualTo(SetIndependence.Reproduction));
        Assert.That(h1[1].Independence, Is.EqualTo(SetIndependence.Independent));
        Assert.That(freeze.Parts[1].Contrasts.Single().Independence, Is.EqualTo(SetIndependence.Independent));
        Assert.That(freeze.Parts[2].Contrasts.Single().PostHoc, Is.True);
        Assert.That(freeze.Parts[1].Contrasts.Single().PostHoc, Is.False);
        Assert.That(freeze.Parts[3].Contrasts.Single().Independence, Is.EqualTo(SetIndependence.Reproduction));
        Assert.That(freeze.Rows.All(r => r.Primary == true), Is.True);
        Assert.That(freeze.Rows.Any(r => r.Maturation), Is.False);
    }

    [Test]
    public void TheRealExcerpt_ParsesContextsAndMemberAccessions()
    {
        var h1 = SetTestFreeze.Read(RealExcerpt(), "hypotheses-v4").Parts[0];

        var treatment = h1.Contrasts[1].Contexts.Single();
        Assert.That(treatment.Factor, Is.EqualTo("treatment"));
        Assert.That(treatment.Levels, Is.EqualTo(new[] { "UBX0101", "vehicle" }));
        Assert.That(h1.Contrasts[0].Contexts, Is.Empty);
        Assert.That(h1.Contrasts[2].Rows.Single().MemberAccessions, Is.EqualTo(new[] { "P01878" }));
        Assert.That(h1.Contrasts[0].Rows.Single(r => r.Member == "Igha").MemberAccessions, Is.Empty);
        Assert.That(h1.Contrasts[0].Rows[0].Measure, Is.EqualTo(SetTestFreeze.ProteinMeasure));
    }

    [Test]
    public void Sha256_IsTakenOverTheBytesRead()
    {
        byte[] lf = Encoding.UTF8.GetBytes(Header + "\n" + Row() + "\n");
        byte[] crlf = Encoding.UTF8.GetBytes(Header + "\r\n" + Row() + "\r\n");

        var a = SetTestFreeze.Read(lf, "t");
        var b = SetTestFreeze.Read(crlf, "t");

        Assert.That(a.Sha256, Is.EqualTo(Convert.ToHexString(SHA256.HashData(lf)).ToLowerInvariant()));
        Assert.That(b.Sha256, Is.Not.EqualTo(a.Sha256));
        Assert.That(b.Rows.Single().PostHoc, Is.EqualTo(a.Rows.Single().PostHoc));
        Assert.That(b.Rows.Single().IndependenceText, Is.EqualTo("independent"));
        Assert.That(a.Freeze, Is.EqualTo("t"));
    }

    [Test]
    public void AByteOrderMark_IsSkipped()
    {
        byte[] withBom = new byte[] { 0xEF, 0xBB, 0xBF }.Concat(Encoding.UTF8.GetBytes(Header + "\n" + Row() + "\n")).ToArray();
        Assert.That(SetTestFreeze.Read(withBom, "t").Rows.Single().HypothesisId, Is.EqualTo("H9"));
    }

    [TestCase("independent", SetIndependence.Independent)]
    [TestCase("independent (sole liver dataset)", SetIndependence.Independent)]
    [TestCase("reproduces [102]", SetIndependence.Reproduction)]
    [TestCase("probably reproduces [102]", SetIndependence.Reproduction)]
    [TestCase("weak (not a test)", SetIndependence.NotPooled)]
    [TestCase("not suited (secretome)", SetIndependence.NotPooled)]
    public void Independence_IsReadAsRuled(string value, SetIndependence expected) =>
        Assert.That(SetTestFreeze.ReadIndependence(value), Is.EqualTo(expected));

    [Test]
    public void AContrastIsIndependentOnlyIfEveryMemberRowIs_AndPostHocOrMaturationIfAnyIs()
    {
        string header = Header + "\tnote";
        var freeze = Read(header,
            Row(member: "GeneA") + "\t",
            Row(member: "GeneB", independence: "weak (not a test)", postHoc: "yes") + "\tseen",
            Row(member: "GeneC") + "\tMATURATION CONTRAST, reported apart",
            Row(member: "GeneA", dataset: "PXD2") + "\tnote",
            Row(member: "GeneB", dataset: "PXD2") + "\tmaturation contrast, not the prefix");

        var first = freeze.Parts.Single().Contrasts[0];
        var second = freeze.Parts.Single().Contrasts[1];
        Assert.That(first.Independence, Is.EqualTo(SetIndependence.NotPooled));
        Assert.That(first.PostHoc, Is.True);
        Assert.That(first.Maturation, Is.True);
        Assert.That(first.Rows.Select(r => r.Maturation), Is.EqualTo(new[] { false, false, true }));
        Assert.That(second.Independence, Is.EqualTo(SetIndependence.Independent));
        Assert.That(second.PostHoc, Is.False);
        // The prefix is matched exactly, in its case.
        Assert.That(second.Maturation, Is.False);
    }

    [Test]
    public void AReproductionAmongNonIndependentRows_MakesTheContrastAReproduction()
    {
        var freeze = Read(Header,
            Row(member: "GeneA", independence: "weak (not a test)"),
            Row(member: "GeneB", independence: "reproduces [9]"));
        Assert.That(freeze.Parts.Single().Contrasts.Single().Independence, Is.EqualTo(SetIndependence.Reproduction));
    }

    [Test]
    public void OptionalColumns_AreReadWhenPresent()
    {
        string header = Header + "\tmember_accessions\tmargin_log2\tattenuation_fraction\tmargin_reason\tprimary\tmeasure\tisoform_accession\tmarker_region";
        var row = SetTestFreeze.Read(Encoding.UTF8.GetBytes(header + "\n" +
                Row(prediction: "attenuated_down", mrna: "-2.21") + "\tP1; P2\t0.25\t0.4\tsmall effect\tno\texon_peptide\tP12023-2\tP12023-2:289-302\n"),
            "t").Rows.Single();

        Assert.That(row.Prediction, Is.EqualTo(SetPrediction.AttenuatedDown));
        Assert.That(row.MrnaLog2, Is.EqualTo(-2.21));
        Assert.That(row.MemberAccessions, Is.EqualTo(new[] { "P1", "P2" }));
        Assert.That(row.MarginLog2, Is.EqualTo(0.25));
        Assert.That(row.AttenuationFraction, Is.EqualTo(0.4));
        Assert.That(row.MarginReason, Is.EqualTo("small effect"));
        Assert.That(row.Primary, Is.False);
        Assert.That(row.Measure, Is.EqualTo("exon_peptide"));
        Assert.That(row.IsoformAccession, Is.EqualTo("P12023-2"));
        Assert.That(row.MarkerRegion, Is.EqualTo("P12023-2:289-302"));
    }

    [Test]
    public void AnAbsentOptionalColumn_GivesItsDefault()
    {
        var row = Read(Header, Row()).Rows.Single();
        Assert.That(row.Primary, Is.Null);
        Assert.That(row.MarginLog2, Is.Null);
        Assert.That(row.MrnaLog2, Is.Null);
        Assert.That(row.MemberAccessions, Is.Empty);
        Assert.That(row.Measure, Is.EqualTo("protein"));
        Assert.That(row.Maturation, Is.False);
    }

    [TestCase("up", SetPrediction.Up)]
    [TestCase("down", SetPrediction.Down)]
    [TestCase("unchanged", SetPrediction.Unchanged)]
    [TestCase("not_up", SetPrediction.NotUp)]
    [TestCase("not_down", SetPrediction.NotDown)]
    [TestCase("attenuated_down", SetPrediction.AttenuatedDown)]
    [TestCase("attenuated_up", SetPrediction.AttenuatedUp)]
    public void EveryPrediction_IsReadAndNamedBack(string text, SetPrediction expected)
    {
        Assert.That(Read(Header, Row(prediction: text)).Rows.Single().Prediction, Is.EqualTo(expected));
        Assert.That(SetTestFreeze.PredictionName(expected), Is.EqualTo(text));
    }

    [Test]
    public void ASingleKindWithSeveralMembers_IsStillASingle_AndASetKindWithOneIsASet()
    {
        var freeze = Read(Header,
            Row(hypothesis: "H1", kind: "single", member: "GeneA"),
            Row(hypothesis: "H1", kind: "single", member: "GeneB"),
            Row(hypothesis: "H2", kind: "set", member: "GeneA"));
        Assert.That(freeze.Parts[0].IsSet, Is.False);
        Assert.That(freeze.Parts[1].IsSet, Is.True);
    }

    [TestCase("", new string[0])]
    [TestCase("treatment=normal", new[] { "treatment:normal" })]
    [TestCase("treatment=UBX0101|treatment=vehicle", new[] { "treatment:UBX0101,vehicle" })]
    [TestCase("treatment=UBX0101|vehicle", new[] { "treatment:UBX0101,vehicle" })]
    [TestCase("ischemia-reperfusion=pre-IRI; preparation=proteome", new[] { "ischemia-reperfusion:pre-IRI", "preparation:proteome" })]
    public void Contexts_AreParsed(string text, string[] expected) =>
        Assert.That(SetTestFreeze.ParseContexts(text).Select(c => $"{c.Factor}:{string.Join(",", c.Levels)}"),
            Is.EqualTo(expected));

    [TestCase("treatment=a|sex=male", "joins levels of one factor")]
    [TestCase("treatment", "names no factor")]
    [TestCase("treatment=a; treatment=b", "names 'treatment' twice")]
    [TestCase("treatment=a;;sex=f", "empty group")]
    [TestCase("treatment=", "empty factor or level")]
    public void AMalformedContextsCell_IsRefused(string text, string reason) =>
        Assert.That(() => SetTestFreeze.ParseContexts(text), Throws.TypeOf<FormatException>().With.Message.Contains(reason));

    [Test]
    public void AMissingRequiredColumn_IsRefusedByName() =>
        Assert.That(() => Read(Header.Replace("\tpost_hoc", ""), Row().Replace("\tno", "")),
            Throws.TypeOf<FormatException>().With.Message.Contains("'post_hoc'"));

    [Test]
    public void AColumnNamedTwice_IsRefused() =>
        Assert.That(() => Read(Header + "\tmember", Row() + "\tGeneB"),
            Throws.TypeOf<FormatException>().With.Message.Contains("named twice"));

    [Test]
    public void ARowWithTheWrongNumberOfCells_IsRefusedWithItsLine() =>
        Assert.That(() => Read(Header, Row(), Row() + "\textra"),
            Throws.TypeOf<FormatException>().With.Message.StartsWith("Line 3:"));

    [TestCase("prediction", "rises", "not one of")]
    [TestCase("post_hoc", "maybe", "write yes or no")]
    [TestCase("mrna_log2", "-1,5", "not a number")]
    [TestCase("member", "", "'member' is empty")]
    public void ABadCell_IsRefusedWithItsLine(string column, string value, string reason)
    {
        var cells = Row().Split('\t');
        cells[Header.Split('\t').ToList().IndexOf(column)] = value;
        Assert.That(() => Read(Header, string.Join('\t', cells)),
            Throws.TypeOf<FormatException>().With.Message.StartsWith("Line 2:").And.Message.Contains(reason));
    }

    [Test]
    public void OnePartWithTwoKinds_IsRefused() =>
        Assert.That(() => Read(Header, Row(kind: "set"), Row(kind: "concordance", member: "GeneB")),
            Throws.TypeOf<FormatException>().With.Message.Contains("one part has one kind"));

    [Test]
    public void OneContrastWithTwoPrimaryMarks_IsRefused() =>
        Assert.That(() => Read(Header + "\tprimary", Row() + "\tyes", Row(member: "GeneB") + "\tno"),
            Throws.TypeOf<FormatException>().With.Message.Contains("different 'primary' marks"));

    [Test]
    public void AnEmptyFileOrAnUnnamedFreeze_IsRefused()
    {
        Assert.That(() => SetTestFreeze.Read(Array.Empty<byte>(), "t"), Throws.TypeOf<FormatException>());
        Assert.That(() => SetTestFreeze.Read(Encoding.UTF8.GetBytes(Header + "\n"), " "), Throws.ArgumentException);
    }
}
