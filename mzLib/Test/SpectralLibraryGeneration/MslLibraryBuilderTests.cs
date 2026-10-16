using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using MassSpectrometry;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.SpectralLibraryGeneration;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.SpectralLibrary;
using UsefulProteomicsDatabases;

namespace Test.SpectralLibraryGeneration;

/// <summary>
/// Builds a predicted <c>.msl</c> library from proteins: digest, generate decoys, predict fragments and iRT with Koina,
/// and label every entry. Offline tests answer Koina from a canned transport. Every peptide gets the same three
/// fragments, and its iRT is 3 × its residue count, so each entry can be checked against its own peptide.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MslLibraryBuilderTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "MslLibraryBuilderTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown() => Directory.Delete(_directory, true);

    private static readonly DigestionParams Trypsin = new(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 4, maxPeptideLength: 50);

    private static MslLibraryBuildParameters Parameters(List<Modification>? fixedMods = null, int chunk = 50_000) =>
        new(Trypsin, fixedMods ?? [], [], PrecursorCharges: [2, 3], Nce: 27, DissociationType: DissociationType.HCD,
            PredictionChunkSize: chunk);

    private static MslLibraryBuilder Builder(CannedKoina? koina = null)
    {
        koina ??= new CannedKoina();
        return new MslLibraryBuilder(koina.Intensity, koina.Irt);
    }

    private static string Strip(string sequence) => Regex.Replace(sequence, @"\[[^\]]*\]", "");

    /// <summary>Each target peptide, at each requested charge, becomes one labelled entry.</summary>
    [Test]
    public void EveryTargetPeptideAtEveryChargeBecomesALabelledEntry()
    {
        var proteins = new List<Protein> { new("MPEPTIDEKAAGGLLRWWSSK", "P1"), new("MNNQQKTTTTVVR", "P2") };

        var entries = Builder().Build(proteins, Parameters(), out var result);

        var expected = proteins.SelectMany(p => p.Digest(Trypsin, [], []))
            .Where(p => p.Length <= 30)
            .SelectMany(p => new[] { 2, 3 }.Select(z => (p.FullSequence, z))).ToHashSet();
        var targets = entries.Where(e => !e.IsDecoy).ToList();
        Assert.That(targets.Select(e => (e.FullSequence, e.ChargeState)).ToHashSet(), Is.EquivalentTo(expected));
        Assert.That(result.TargetPrecursors, Is.EqualTo(expected.Count));
        foreach (var entry in targets)
        {
            Assert.That(entry.Nce, Is.EqualTo(27));
            Assert.That(entry.DissociationType, Is.EqualTo(DissociationType.HCD));
            Assert.That(entry.RetentionTime, Is.EqualTo(3.0 * Strip(entry.FullSequence).Length).Within(1e-4), "RT holds the predicted iRT");
            Assert.That(entry.MatchedFragmentIons, Is.Not.Empty);
            Assert.That(entry.ProteinAccession, Is.EqualTo(proteins.Single(p => p.BaseSequence.Contains(entry.BaseSequence)).Accession));
        }
    }

    /// <summary>Decoys are the reversed proteins' peptides, with DECOY_ accessions, predicted like the targets.</summary>
    [Test]
    public void DecoysComeFromReversedProteinsWithDecoyAccessions()
    {
        var proteins = new List<Protein> { new("MPEPTIDEKAAGGLLRWWSSK", "P1") };

        var entries = Builder().Build(proteins, Parameters(), out var result);

        var decoyProteins = DecoyProteinGenerator.GenerateDecoys(proteins, DecoyType.Reverse);
        var expected = decoyProteins.SelectMany(p => p.Digest(Trypsin, [], []))
            .SelectMany(p => new[] { 2, 3 }.Select(z => (p.FullSequence, z))).ToHashSet();
        var decoys = entries.Where(e => e.IsDecoy).ToList();
        Assert.That(decoys.Select(e => (e.FullSequence, e.ChargeState)).ToHashSet(), Is.EquivalentTo(expected));
        Assert.That(decoys.All(e => e.ProteinAccession == "DECOY_P1"));
        Assert.That(decoys.All(e => e.MatchedFragmentIons.Count > 0 && e.Nce == 27));
        Assert.That(result.DecoyPrecursors, Is.EqualTo(decoys.Count));
    }

    /// <summary>
    /// A decoy peptide identical to a target peptide would be a target in disguise. The target keeps the sequence, and the
    /// decoy is dropped and counted. The reverse of MSSSSKQQQQK keeps its initiator M and yields QQQQK again.
    /// </summary>
    [Test]
    public void ADecoyWithATargetSequenceIsDroppedAndCounted()
    {
        var proteins = new List<Protein> { new("MSSSSKQQQQK", "P1") };

        var entries = Builder().Build(proteins, Parameters(), out var result);

        Assert.That(entries.Count(e => e.BaseSequence == "QQQQK"), Is.EqualTo(2), "one target entry per charge, no decoy");
        Assert.That(entries.Where(e => e.BaseSequence == "QQQQK").All(e => !e.IsDecoy));
        Assert.That(result.DecoysDroppedAsTargetSequences, Is.EqualTo(2));
    }

    /// <summary>
    /// I and L have the same mass, so a decoy equal to a target once I = L has the target's spectrum: it too is a target in
    /// disguise, and when that target is in the sample it scores like a real identification. Reversed, P1 yields QIQQK
    /// (P2's QLQQK) and P2 yields QQLQK (P1's QQIQK); both decoys are dropped and counted.
    /// </summary>
    [Test]
    public void ADecoyEqualToATargetWithIEqualToLIsDroppedAndCounted()
    {
        var proteins = new List<Protein> { new("MSSSSKQQIQK", "P1"), new("MGGGGKQLQQK", "P2") };

        var entries = Builder().Build(proteins, Parameters(), out var result);

        var targets = entries.Where(e => !e.IsDecoy).Select(e => e.BaseSequence.Replace('I', 'L')).ToHashSet();
        Assert.That(entries.Where(e => e.IsDecoy).Select(e => e.BaseSequence.Replace('I', 'L')), Has.None.AnyOf(targets.ToArray()));
        Assert.That(entries.Where(e => !e.IsDecoy).Select(e => e.BaseSequence), Is.SupersetOf(new[] { "QQIQK", "QLQQK" }), "targets keep their sequences");
        Assert.That(result.DecoysDroppedAsTargetSequences, Is.EqualTo(4), "two decoy peptides at two charges");
    }

    /// <summary>A peptide in several proteins lists every accession, sorted and '|'-joined, and that survives a save.</summary>
    [Test]
    public void ASharedPeptideListsEveryAccessionAndSurvivesASave()
    {
        var proteins = new List<Protein> { new("MAAAAKSHAMEDPEPR", "P2"), new("MGGGGKSHAMEDPEPR", "P1") };

        var entries = Builder().Build(proteins, Parameters(), out _);

        var shared = entries.Where(e => !e.IsDecoy && e.BaseSequence == "SHAMEDPEPR").ToList();
        Assert.That(shared, Is.Not.Empty);
        Assert.That(shared.All(e => e.ProteinAccession == "P1|P2"));
        string path = Path.Combine(_directory, "shared.msl");
        MslLibrary.Save(path, entries);
        using var library = MslLibrary.Load(path);
        Assert.That(library.GetAllEntries().Where(e => !e.IsDecoy && e.BaseSequence == "SHAMEDPEPR").Select(e => e.ProteinAccession),
            Is.All.EqualTo("P1|P2"));
    }

    /// <summary>Fixed modifications reach the library with their labels and masses (Koina's modified-peptide path).</summary>
    [Test]
    public void ModifiedPeptidesKeepTheirModificationsAndMasses()
    {
        var carbamidomethyl = Mods.AllKnownProteinModsDictionary["Carbamidomethyl on C"];
        var proteins = new List<Protein> { new("MPEPCTIDEKAAGGLLR", "P1") };

        var entries = Builder().Build(proteins, Parameters([carbamidomethyl]), out _);

        var digested = proteins[0].Digest(Trypsin, [carbamidomethyl], []).Single(p => p.BaseSequence == "PEPCTIDEK");
        var entry = entries.Single(e => !e.IsDecoy && e.BaseSequence == "PEPCTIDEK" && e.ChargeState == 2);
        Assert.That(entry.FullSequence, Is.EqualTo(digested.FullSequence));
        Assert.That(entry.PrecursorMz, Is.EqualTo(Chemistry.ClassExtensions.ToMz(digested.MonoisotopicMass, 2)).Within(1e-4));
    }

    /// <summary>A peptide the model cannot take (longer than 30) is counted, never thrown.</summary>
    [Test]
    public void PeptidesTheModelCannotTakeAreCountedNotThrown()
    {
        string tooLong = new string('A', 34) + "K";
        var proteins = new List<Protein> { new("MPEPTIDEK" + tooLong + "GGLLR", "P1") };

        var entries = Builder().Build(proteins, Parameters(), out var result);

        Assert.That(entries.Any(e => e.BaseSequence == tooLong), Is.False);
        Assert.That(result.NotPredicted, Is.GreaterThanOrEqualTo(2), "the long peptide at charges 2 and 3");
    }

    /// <summary>Entrapment proteins from the hook are digested, decoyed and predicted like targets, keeping their accessions.</summary>
    [Test]
    public void EntrapmentProteinsFromTheHookAreBuiltLikeTargets()
    {
        var proteins = new List<Protein> { new("MPEPTIDEKAAGGLLR", "P1") };
        var builder = new MslLibraryBuilder(new CannedKoina().Intensity, new CannedKoina().Irt)
        {
            EntrapmentProteins = targets => targets.Select(p => new Protein("MWWYYKFFHHR", "Random_" + p.Accession)).ToList(),
        };

        var entries = builder.Build(proteins, Parameters(), out _);

        Assert.That(entries.Any(e => !e.IsDecoy && e.ProteinAccession == "Random_P1" && e.BaseSequence == "FFHHR"));
        Assert.That(entries.Any(e => e.IsDecoy && e.ProteinAccession == "DECOY_Random_P1"));
    }

    /// <summary>Koina is called in chunks, so a proteome-scale build never sends one enormous session.</summary>
    [Test]
    public void PredictionsAreRequestedInChunks()
    {
        var koina = new CannedKoina();
        var proteins = new List<Protein> { new("MPEPTIDEKAAGGLLRWWSSKNNQQKTTTTVVR", "P1") };

        var chunked = Builder(koina).Build(proteins, Parameters(chunk: 3), out _);
        int chunkedCalls = koina.IntensityCalls;
        var whole = Builder().Build(proteins, Parameters(), out _);

        Assert.That(chunkedCalls, Is.GreaterThan(3));
        Assert.That(chunked.Select(e => (e.FullSequence, e.ChargeState, e.IsDecoy)), Is.EquivalentTo(whole.Select(e => (e.FullSequence, e.ChargeState, e.IsDecoy))));
    }

    [Test]
    public void ACancelledBuildStopsBeforeCallingKoina()
    {
        var koina = new CannedKoina();
        using var cancelled = new CancellationTokenSource();
        cancelled.Cancel();

        Assert.Throws<OperationCanceledException>(() =>
            Builder(koina).Build([new Protein("MPEPTIDEKAAGGLLR", "P1")], Parameters(), out _, cancelled.Token));
        Assert.That(koina.IntensityCalls, Is.Zero);
    }

    [Test]
    public void ArgumentsAreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => new MslLibraryBuilder(null!, new CannedKoina().Irt));
        Assert.Throws<ArgumentNullException>(() => new MslLibraryBuilder(new CannedKoina().Intensity, null!));
        Assert.Throws<ArgumentException>(() => Builder().Build([new Protein("MPEPTIDEK", "P1")], Parameters() with { PrecursorCharges = [] }, out _));
        Assert.Throws<ArgumentOutOfRangeException>(() => Builder().Build([new Protein("MPEPTIDEK", "P1")], Parameters(chunk: 0), out _));
    }

    /// <summary>Live: one small protein through the real Prosit 2025 40PTM models.</summary>
    [Test]
    [Category("ExternalService")]
    [Category("Koina")]
    public async Task LiveBuildAgainstKoina()
    {
        await ExternalServiceTestHelper.RunAsync("Koina", () =>
        {
            var carbamidomethyl = Mods.AllKnownProteinModsDictionary["Carbamidomethyl on C"];
            var builder = new MslLibraryBuilder(new Prosit2025Intensity40PTM(), new Prosit2025iRT40PTM());
            var entries = builder.Build([new Protein("MSLPDKAGTCVLVEGDRPYVLNSR", "P1")], Parameters([carbamidomethyl]), out var result);

            Assert.That(result.TargetPrecursors, Is.GreaterThan(0));
            Assert.That(entries.Where(e => !e.IsDecoy), Has.All.Matches<MslLibraryEntry>(e => e.MatchedFragmentIons.Count > 0 && double.IsFinite(e.RetentionTime)));
            return Task.CompletedTask;
        });
    }

    /// <summary>Koina answered from memory: three fragments per peptide, and iRT = 3 × residue count.</summary>
    private sealed class CannedKoina
    {
        public int IntensityCalls;
        public FakeIntensity Intensity { get; }
        public FakeIrt Irt { get; }

        public CannedKoina()
        {
            Intensity = new FakeIntensity(this);
            Irt = new FakeIrt();
        }

        public static string[] Peptides(Dictionary<string, object> request) =>
            JObject.Parse(JsonConvert.SerializeObject(request))["inputs"]![0]!["data"]!.Select(t => t.ToString()).ToArray();

        public sealed class FakeIntensity(CannedKoina owner) : Prosit2025Intensity40PTM
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
            {
                Interlocked.Increment(ref owner.IntensityCalls);
                int n = Peptides(request).Length;
                string Repeat(string items) => string.Join(",", Enumerable.Repeat(items, n));
                return Task.FromResult("{\"outputs\":[" +
                    $"{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{n},3],\"data\":[{Repeat("\"y1+1\",\"y2+1\",\"b2+1\"")}]}}," +
                    $"{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{n},3],\"data\":[{Repeat("100.0,200.0,300.0")}]}}," +
                    $"{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{n},3],\"data\":[{Repeat("0.2,1.0,0.5")}]}}]}}");
            }
        }

        public sealed class FakeIrt : Prosit2025iRT40PTM
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
            {
                var irts = Peptides(request).Select(p => (3.0 * Strip(p).Replace("-", "").Length).ToString(System.Globalization.CultureInfo.InvariantCulture));
                var data = string.Join(",", irts);
                return Task.FromResult($"{{\"outputs\":[{{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[{irts.Count()}],\"data\":[{data}]}}]}}");
            }
        }
    }
}
