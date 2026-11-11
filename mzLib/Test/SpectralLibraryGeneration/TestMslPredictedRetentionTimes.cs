using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using MassSpectrometry;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SequenceConversion;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.SpectralLibraryGeneration;

namespace Test.SpectralLibraryGeneration;

/// <summary>
/// Tests for <see cref="MslPredictedRetentionTimes"/>: predicted iRT for a library's entries (both modes),
/// the Koina adapter's ProForma route, and Chronologer mapped onto the iRT scale.
/// </summary>
[TestFixture]
public sealed class TestMslPredictedRetentionTimes
{
	private const string Phospho = "[Common Biological:Phosphorylation on T]";
	private const MslFormat.RetentionTimeStatus None = MslFormat.RetentionTimeStatus.None;

	private static MslLibraryEntry Entry(string fullSequence, string baseSequence, int charge = 2, double rt = double.NaN,
		MslFormat.RetentionTimeStatus status = None, bool struck = false, bool rtIsCalibrated = false) => new()
	{
		FullSequence = fullSequence,
		BaseSequence = baseSequence,
		ChargeState = charge,
		PrecursorMz = 500,
		RetentionTime = rt,
		RtStatus = status,
		RtStruckByPrediction = struck,
		RtIsCalibrated = rtIsCalibrated,
		Source = MslFormat.SourceType.Empirical,
		MatchedFragmentIons = new List<MslFragmentIon>
		{
			new() { Mz = 300f, Intensity = 1f, ProductType = ProductType.y, FragmentNumber = 3, Charge = 1 },
		}
	};

	/// <summary>A predictor that records what it is asked and answers from a function.</summary>
	private sealed class Recorder(Func<string, double> answer)
	{
		public readonly List<IReadOnlyList<string>> Calls = new();
		public IReadOnlyList<double> Predict(IReadOnlyList<string> sequences)
		{
			Calls.Add(sequences.ToList());
			return sequences.Select(answer).ToList();
		}
	}

	private static readonly string ShortPeptide = "PEPTIDESK";                                   // 9
	private static readonly string ShortModified = "PEPT" + Phospho + "IDESK";
	private static readonly string LongPeptide = new string('A', 20) + new string('L', 14) + "K"; // 35
	private static readonly string LongModified = "PEPT" + Phospho + new string('A', 30) + "K";   // 35
	private static readonly string VeryLong = new string('G', 54) + "K";                           // 55

	// ── Assign ─────────────────────────────────────────────────────────────────

	[Test]
	public void ReplaceAll_GivesEveryEntryItsPrediction_ByLengthAndModification()
	{
		var entries = new List<MslLibraryEntry>
		{
			Entry(ShortPeptide, ShortPeptide, charge: 2, rt: 37.0, status: MslFormat.RetentionTimeStatus.Confirmed, rtIsCalibrated: true),
			Entry(ShortPeptide, ShortPeptide, charge: 3, rt: 37.0),
			Entry(ShortModified, ShortPeptide, rt: 12.0, status: MslFormat.RetentionTimeStatus.Challenged, struck: true),
			Entry(LongPeptide, LongPeptide),
			Entry(LongModified, "PEPT" + new string('A', 30) + "K"),
			Entry(VeryLong, VeryLong),
		};
		var irt = new Recorder(s => s.Length);            // a stand-in value per sequence
		var longIrt = new Recorder(s => 1000);

		var output = MslPredictedRetentionTimes.Assign(entries, irt.Predict, MslPredictedIrtMode.ReplaceAll, out var result, longIrt.Predict);

		Assert.That(irt.Calls.Single(), Is.EquivalentTo(new[] { ShortPeptide, ShortModified }), "one request, one value per full sequence");
		Assert.That(longIrt.Calls.Single(), Is.EqualTo(new[] { LongPeptide }), "only unmodified 31-50 residues");
		Assert.That(output.Select(e => e.RetentionTime),
			Is.EqualTo(new[] { ShortPeptide.Length, ShortPeptide.Length, ShortModified.Length, 1000, double.NaN, double.NaN }));
		Assert.That(output.Select(e => e.RtStatus), Is.All.EqualTo(None));
		Assert.That(output.Select(e => e.RtStruckByPrediction), Is.All.False);
		Assert.That(output.Select(e => e.RtIsCalibrated), Is.All.False);
		Assert.That(output.Select(e => e.Source), Is.All.EqualTo(MslFormat.SourceType.Empirical), "the source is the caller's to set");
		Assert.That(result, Is.EqualTo(new MslPredictedIrtResult(6, 3, 1, 2, 3)));

		// The input entries are not modified
		Assert.That(entries[0].RetentionTime, Is.EqualTo(37.0));
		Assert.That(entries[0].RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.Confirmed));
	}

	[Test]
	public void FillMissing_TouchesOnlyEntriesWithNoTimeAndNoStatus_AndCanBeRepeated()
	{
		var entries = new List<MslLibraryEntry>
		{
			Entry(ShortPeptide, ShortPeptide),                                                        // no time: filled
			Entry(ShortModified, ShortPeptide),                                                       // no time, no prediction yet
			Entry("MEASUREDK", "MEASUREDK", rt: 50, status: MslFormat.RetentionTimeStatus.ObservedOnce), // measured: untouched
			Entry("PREDICTEDK", "PREDICTEDK", rt: 20),                                                // a predicted time: untouched
		};
		var outage = new Recorder(s => s == ShortPeptide ? 11.0 : double.NaN);

		var first = MslPredictedRetentionTimes.Assign(entries, outage.Predict, MslPredictedIrtMode.FillMissing, out var r1);
		Assert.That(outage.Calls.Single(), Is.EquivalentTo(new[] { ShortPeptide, ShortModified }));
		Assert.That(first.Select(e => e.RetentionTime), Is.EqualTo(new[] { 11.0, double.NaN, 50, 20 }));
		Assert.That(first[2].RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.ObservedOnce));
		Assert.That(ReferenceEquals(first[2], entries[2]) && ReferenceEquals(first[3], entries[3]), Is.True);
		Assert.That(r1, Is.EqualTo(new MslPredictedIrtResult(2, 1, 0, 1, 2)));

		// Repeated after the outage: only the entry still without a time is asked for
		var recovered = new Recorder(_ => 33.0);
		var second = MslPredictedRetentionTimes.Assign(first, recovered.Predict, MslPredictedIrtMode.FillMissing, out var r2);
		Assert.That(recovered.Calls.Single(), Is.EqualTo(new[] { ShortModified }));
		Assert.That(second.Select(e => e.RetentionTime), Is.EqualTo(new[] { 11.0, 33.0, 50, 20 }));
		Assert.That(r2.Considered, Is.EqualTo(1));
	}

	[Test]
	public void WithoutALongPeptidePredictor_LongPeptidesHaveNoTime()
	{
		var irt = new Recorder(_ => 5.0);
		var output = MslPredictedRetentionTimes.Assign(new List<MslLibraryEntry> { Entry(LongPeptide, LongPeptide) },
			irt.Predict, MslPredictedIrtMode.ReplaceAll, out var result);

		Assert.That(irt.Calls, Is.Empty, "nothing to ask the iRT predictor");
		Assert.That(output.Single().RetentionTime, Is.NaN);
		Assert.That(result.NoPrediction, Is.EqualTo(1));
	}

	[Test]
	public void Assign_ChecksItsArguments_AndThePredictorsContract()
	{
		var entries = new List<MslLibraryEntry> { Entry(ShortPeptide, ShortPeptide) };
		Func<IReadOnlyList<string>, IReadOnlyList<double>> ok = s => s.Select(_ => 1.0).ToList();

		Assert.That(() => MslPredictedRetentionTimes.Assign(null!, ok, MslPredictedIrtMode.ReplaceAll, out _), Throws.ArgumentNullException);
		Assert.That(() => MslPredictedRetentionTimes.Assign(entries, null!, MslPredictedIrtMode.ReplaceAll, out _), Throws.ArgumentNullException);
		Assert.That(() => MslPredictedRetentionTimes.Assign(entries, ok, (MslPredictedIrtMode)9, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
		Assert.That(() => MslPredictedRetentionTimes.Assign(new List<MslLibraryEntry> { null! }, ok, MslPredictedIrtMode.ReplaceAll, out _),
			Throws.ArgumentException);
		Assert.That(() => MslPredictedRetentionTimes.Assign(entries, s => new double[s.Count + 1], MslPredictedIrtMode.ReplaceAll, out _),
			Throws.InvalidOperationException);
		Assert.That(() => MslPredictedRetentionTimes.Assign(entries, _ => null!, MslPredictedIrtMode.ReplaceAll, out _),
			Throws.InvalidOperationException);
	}

	// ── FromKoina ──────────────────────────────────────────────────────────────

	/// <summary>Prosit 2025 iRT answered from memory: 3 x the residues of each sequence it is sent; records what it was sent.</summary>
	private sealed class FakeProsit : Prosit2025iRT40PTM
	{
		public readonly List<string> Sent = new();

		protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
		{
			string[] peptides = JObject.Parse(JsonConvert.SerializeObject(request))["inputs"]![0]!["data"]!.Select(t => t.ToString()).ToArray();
			lock (Sent) Sent.AddRange(peptides);
			var irts = peptides.Select(p => (3.0 * System.Text.RegularExpressions.Regex.Replace(p, @"\[[^\]]*\]", "").Replace("-", "").Length)
				.ToString(System.Globalization.CultureInfo.InvariantCulture));
			return Task.FromResult($"{{\"outputs\":[{{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[{peptides.Length}],\"data\":[{string.Join(",", irts)}]}}]}}");
		}
	}

	[Test]
	public void FromKoina_SendsProForma_AndReturnsOneValuePerSequenceInOrder()
	{
		var model = new FakeProsit();
		var predict = MslPredictedRetentionTimes.FromKoina(model);
		string unknownModification = "PEPT[Nonexistent:Nothing on T]IDEK";

		IReadOnlyList<double> values = predict(new[] { ShortModified, unknownModification, ShortPeptide, LongPeptide });

		Assert.That(values[0], Is.EqualTo(3.0 * ShortPeptide.Length), "the phosphopeptide");
		Assert.That(values[1], Is.NaN, "a modification ProForma cannot name is never sent");
		Assert.That(values[2], Is.EqualTo(3.0 * ShortPeptide.Length));
		Assert.That(values[3], Is.NaN, "over Prosit's 30 residues");
		Assert.That(model.Sent, Has.Some.Contains("[UNIMOD:21]"), "the phosphorylation reaches the model by its Unimod id");
		Assert.That(model.Sent, Has.None.Contains("Nonexistent"));
	}

	[Test]
	public void FromKoina_ChecksItsArguments()
	{
		Assert.That(() => MslPredictedRetentionTimes.FromKoina(null!), Throws.ArgumentNullException);
		Assert.That(() => MslPredictedRetentionTimes.FromKoina(new ChronologerRT()), Throws.ArgumentException, "a model of minutes, not iRT");
		Assert.That(() => MslPredictedRetentionTimes.FromKoina(new FakeProsit(), chunkSize: 0), Throws.TypeOf<ArgumentOutOfRangeException>());
	}

	// ── ChronologerOnIrtScale ──────────────────────────────────────────────────

	private sealed record Bare(string FullSequence, string BaseSequence) : IRetentionPredictable
	{
		public double MonoisotopicMass => 0;
		public string FullSequenceWithMassShifts => FullSequence;
	}

	private static double ChronologerAcn(string sequence)
	{
		using var chronologer = new ChronologerRetentionTimePredictor(SequenceConversionHandlingMode.ReturnNull);
		return chronologer.PredictRetentionTimeEquivalents(new IRetentionPredictable[] { new Bare(sequence, sequence) }).Single().PredictedValue
			?? double.NaN;
	}

	[Test]
	public void ChronologerOnIrtScale_MapsChronologerOntoTheKnownScale()
	{
		// 60 unmodified peptides whose "known" iRT is a line through Chronologer's own predictions
		const string residues = "ADEFGHIKLNPQRSTVWY";
		var peptides = Enumerable.Range(0, 60)
			.Select(i => new string(Enumerable.Range(0, 8 + i % 15).Select(k => residues[(i * 7 + k * 5 + k * k) % 18]).ToArray()) + "K")
			.Distinct()
			.ToList();
		using var chronologer = new ChronologerRetentionTimePredictor(SequenceConversionHandlingMode.ReturnNull);
		var acn = chronologer.PredictRetentionTimeEquivalents(peptides.Select(p => (IRetentionPredictable)new Bare(p, p)).ToList())
			.Where(r => r.PredictedValue is double v && double.IsFinite(v))
			.ToDictionary(r => r.Peptide.FullSequence, r => r.PredictedValue!.Value);
		var known = acn.ToDictionary(kv => kv.Key, kv => 2 * kv.Value + 5);
		known["PEPT" + Phospho + "IDEK"] = 1000; // modified: not used for the fit

		var predict = MslPredictedRetentionTimes.ChronologerOnIrtScale(known);
		string middle = acn.OrderBy(kv => kv.Value).ElementAt(acn.Count / 2).Key;
		IReadOnlyList<double> mapped = predict(new[] { middle, "PEPT" + Phospho + "IDEK", VeryLong });

		Assert.That(mapped[0], Is.EqualTo(2 * acn[middle] + 5).Within(1.0));
		Assert.That(mapped[1], Is.NaN, "modified");
		Assert.That(mapped[2], Is.NaN, "over 50 residues");

		// A 35-residue peptide, past Prosit's limit; the known scale is a line, so even past the fitted range
		// (the fit extrapolates along its end slope) it maps onto the same line
		double longAcn = ChronologerAcn(LongPeptide);
		Assert.That(longAcn, Is.Not.NaN);
		Assert.That(predict(new[] { LongPeptide })[0], Is.EqualTo(2 * longAcn + 5).Within(1.0));
	}

	[Test]
	public void ChronologerOnIrtScale_NeedsEnoughKnownPeptides()
	{
		Assert.That(() => MslPredictedRetentionTimes.ChronologerOnIrtScale(null!), Throws.ArgumentNullException);
		Assert.That(() => MslPredictedRetentionTimes.ChronologerOnIrtScale(new Dictionary<string, double>(), maximumFitPeptides: 5),
			Throws.TypeOf<ArgumentOutOfRangeException>());
		Assert.That(() => MslPredictedRetentionTimes.ChronologerOnIrtScale(new Dictionary<string, double> { ["PEPTIDEK"] = 10 }),
			Throws.ArgumentException);
	}
}
