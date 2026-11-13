using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Tests for <see cref="MslLibrary.UpdateRetentionTimesAndSave"/>. Each dataset is simulated: 400
/// unmodified peptides of 7-30 residues with a known true iRT, runs whose minutes follow a curved,
/// run-shifted map of iRT, and a predictor that returns the true iRT with a fixed error of up to 3 iRT
/// (so calibration behaves as on real data). Single peptides are then moved to exercise each rule.
/// </summary>
[TestFixture]
public sealed class TestMslRetentionTimeUpdate
{
	private const int PeptideCount = 400;
	private const string Phospho = "[Common Biological:Phosphorylation on T]";
	private const MslFormat.RetentionTimeStatus None = MslFormat.RetentionTimeStatus.None;
	private const MslFormat.RetentionTimeStatus ObservedOnce = MslFormat.RetentionTimeStatus.ObservedOnce;
	private const MslFormat.RetentionTimeStatus Confirmed = MslFormat.RetentionTimeStatus.Confirmed;
	private const MslFormat.RetentionTimeStatus Challenged = MslFormat.RetentionTimeStatus.Challenged;

	private string _dir = null!;
	private int _fileCounter;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_dir = Path.Combine(Path.GetTempPath(), "MslRetentionTimeUpdate_" + Guid.NewGuid().ToString("N"));
		Directory.CreateDirectory(_dir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		try { Directory.Delete(_dir, recursive: true); } catch (IOException) { }
	}

	private string NewPath(string name) => Path.Combine(_dir, $"{name}_{_fileCounter++}.msl");

	// ── The simulated peptides ─────────────────────────────────────────────────

	private const string Residues = "ADEFGHIKLNPQRSTVWY"; // no C or M: no fixed or variable modification

	/// <summary>A unique base sequence of 7 + (i mod 24) residues (7-30).</summary>
	private static string Peptide(int i)
	{
		int length = 7 + i % 24;
		var chars = new char[length];
		// The first three residues spell i in base 18, so every peptide is distinct
		chars[0] = Residues[i / 324 % 18];
		chars[1] = Residues[i / 18 % 18];
		chars[2] = Residues[i % 18];
		for (int k = 3; k < length - 1; k++)
			chars[k] = Residues[(i * 7 + k * 13 + k * k) % 18];
		chars[^1] = 'K';
		return new string(chars);
	}

	/// <summary>The true iRT of peptide i, spread over -20 to 120.</summary>
	private static double TrueIrt(int i) => -20 + 140 * Frac(i * 0.6180339887);

	private static double Frac(double x) => x - Math.Floor(x);

	/// <summary>A fixed per-sequence prediction error, uniform in [-3, 3] iRT (FNV-1a hash of the sequence).</summary>
	private static double PredictionError(string sequence)
	{
		uint h = 2166136261;
		foreach (char c in sequence)
		{
			h ^= c;
			h *= 16777619;
		}
		return 6 * (h / (double)uint.MaxValue - 0.5);
	}

	/// <summary>A run's minutes for an iRT: curved, and shifted per run.</summary>
	private static double Minutes(double irt, int run) =>
		8 + 1.5 * run + 0.45 * (irt + 25) + 0.0008 * (irt + 25) * (irt + 25);

	/// <summary>The 400 peptides and their true iRT.</summary>
	private static readonly IReadOnlyDictionary<string, double> Base =
		Enumerable.Range(0, PeptideCount).ToDictionary(Peptide, TrueIrt);

	private static MslLibraryEntry Entry(string fullSequence, string? baseSequence = null, double rt = double.NaN,
		MslFormat.RetentionTimeStatus status = None, bool struck = false, bool isDecoy = false) => new()
	{
		FullSequence = fullSequence,
		BaseSequence = baseSequence ?? fullSequence,
		PrecursorMz = 400 + (fullSequence.GetHashCode() & 0xFFFF) / 100.0,
		ChargeState = 2,
		RetentionTime = rt,
		RtStatus = status,
		RtStruckByPrediction = struck,
		IsDecoy = isDecoy,
		Source = MslFormat.SourceType.Predicted,
		DissociationType = DissociationType.HCD,
		MatchedFragmentIons = new List<MslFragmentIon>
		{
			new() { Mz = 300f, Intensity = 1f, ProductType = ProductType.y, FragmentNumber = 3, Charge = 1 },
			new() { Mz = 400f, Intensity = 0.5f, ProductType = ProductType.b, FragmentNumber = 4, Charge = 1 },
		}
	};

	/// <summary>
	/// 100 measured times that agree with the truth (peptides 200-299), so a test's one disagreeing time is
	/// a small share of the dataset's overlaps.
	/// </summary>
	private static MslLibraryEntry[] AgreeingBackground() =>
		Enumerable.Range(200, 100).Select(i => Entry(Peptide(i), rt: TrueIrt(i), status: ObservedOnce)).ToArray();

	/// <summary>A library of the 400 peptides with no times, with <paramref name="extra"/> replacing or adding entries.</summary>
	private string Library(params MslLibraryEntry[] extra)
	{
		var replaced = extra.Select(e => e.FullSequence).ToHashSet();
		var entries = Base.Keys.Where(p => !replaced.Contains(p)).Select(p => Entry(p)).ToList();
		entries.AddRange(extra);
		string path = NewPath("library");
		MslWriter.Write(path, entries);
		return path;
	}

	/// <summary>
	/// <paramref name="runCount"/> runs that see the 400 peptides at their true iRT (or where
	/// <paramref name="moved"/> puts them), plus <paramref name="extra"/> in every run and
	/// <paramref name="onlyInRun0"/> in the first run only.
	/// </summary>
	private static List<MslRunObservations> Dataset(int runCount = 3, Dictionary<string, double>? moved = null,
		Dictionary<string, (string Base, double Irt)>? extra = null, Dictionary<string, double>? onlyInRun0 = null)
	{
		var runs = new List<MslRunObservations>();
		for (int r = 0; r < runCount; r++)
		{
			var psms = new List<MslRunPsm>();
			foreach (var (seq, irt) in Base)
				psms.Add(new MslRunPsm(seq, seq, 2, Minutes(moved != null && moved.TryGetValue(seq, out double m) ? m : irt, r)));
			if (extra != null)
				foreach (var (full, (bas, irt)) in extra)
					psms.Add(new MslRunPsm(full, bas, 2, Minutes(irt, r)));
			if (r == 0 && onlyInRun0 != null)
				foreach (var (seq, irt) in onlyInRun0)
					psms.Add(new MslRunPsm(seq, seq, 2, Minutes(irt, r)));
			runs.Add(new MslRunObservations($"run{r}", psms, RunLengthMinutes: 100));
		}
		return runs;
	}

	/// <summary>Predicts the true iRT plus a fixed error; NaN for a sequence it does not know.</summary>
	private sealed class Predictor
	{
		public readonly Dictionary<string, double> Known = new(Base);
		public int Calls;
		public double ExtraNoise;

		public IReadOnlyList<double> Predict(IReadOnlyList<string> sequences)
		{
			Calls++;
			return sequences.Select(s => Known.TryGetValue(s, out double irt)
				? irt + PredictionError(s) + ExtraNoise * (2 * Frac(s.Length * 0.37 + s[0] * 0.11 + s[^2] * 0.07) - 1)
				: double.NaN).ToList();
		}
	}

	private (MslRetentionTimeUpdateResult Result, Dictionary<string, MslLibraryEntry> Entries, string Output) Update(
		string libraryPath, List<MslRunObservations> runs, Predictor? predictor = null,
		Func<IReadOnlyList<string>, IReadOnlyList<double>>? longPredictor = null,
		MslRetentionTimeUpdateOptions? options = null, bool indexOnly = false)
	{
		predictor ??= new Predictor();
		string output = NewPath("updated");
		MslRetentionTimeUpdateResult result;
		using (MslLibrary lib = indexOnly ? MslLibrary.LoadIndexOnly(libraryPath) : MslLibrary.Load(libraryPath))
			result = lib.UpdateRetentionTimesAndSave(runs, predictor.Predict, output, longPredictor, options);

		using MslLibrary updated = MslLibrary.Load(output);
		return (result, updated.GetAllEntries().ToDictionary(e => e.FullSequence), output);
	}

	// ── Calibration and the first dataset ──────────────────────────────────────

	[TestCase(false)]
	[TestCase(true)]
	public void FirstDataset_FillsEveryTimeNearTheTruth_AsObservedOnce(bool indexOnly)
	{
		var predictor = new Predictor();
		var (result, entries, _) = Update(Library(), Dataset(), predictor, indexOnly: indexOnly);

		Assert.That(result.Runs.Select(r => r.Calibrated), Is.All.True);
		Assert.That(result.Runs.Select(r => r.AnchorsSelected), Is.All.EqualTo(PeptideCount));
		Assert.That(result.Runs.Select(r => r.ResidualSd), Is.All.LessThan(3));
		Assert.That(result.Filled, Is.EqualTo(PeptideCount));
		Assert.That(result.EntriesChanged, Is.EqualTo(PeptideCount));
		Assert.That(predictor.Calls, Is.EqualTo(1), "every run's anchors in one request, and nothing else needed");

		foreach (var (seq, irt) in Base)
		{
			MslLibraryEntry e = entries[seq];
			Assert.That(e.RetentionTime, Is.EqualTo(irt).Within(1.0), seq);
			Assert.That(e.RtStatus, Is.EqualTo(ObservedOnce), seq);
			Assert.That(e.RtStruckByPrediction, Is.False, seq);
			Assert.That(e.RtIsCalibrated, Is.False, seq);
			Assert.That(e.Source, Is.EqualTo(MslFormat.SourceType.EmpiricalRefined), seq);
		}
	}

	[Test]
	public void TheLoadedLibrary_IsNotModified()
	{
		using MslLibrary lib = MslLibrary.Load(Library());
		lib.UpdateRetentionTimesAndSave(Dataset(), new Predictor().Predict, NewPath("out"));

		Assert.That(lib.GetAllEntries().Select(e => e.RtStatus), Is.All.EqualTo(None));
		Assert.That(lib.GetAllEntries().Select(e => e.RetentionTime), Is.All.NaN);
	}

	[Test]
	public void InputOrder_DoesNotChangeTheResult()
	{
		string path = Library();
		List<MslRunObservations> given = Dataset();
		var shuffled = given
			.Select(r => r with { Psms = r.Psms.OrderBy(p => p.FullSequence.GetHashCode() ^ 0x5bd1e995).ToList() })
			.Reverse()
			.ToList();

		// 150 anchors of 400 candidates, so which peptides are chosen depends on the selection rule alone
		var options = new MslRetentionTimeUpdateOptions(AnchorsPerRun: 150);
		var (_, _, a) = Update(path, given, options: options);
		var (_, _, b) = Update(path, shuffled, options: options);

		Assert.That(File.ReadAllBytes(b), Is.EqualTo(File.ReadAllBytes(a)));
	}

	[Test]
	public void AnchorsPerRun_CapsTheAnchors()
	{
		var (result, _, _) = Update(Library(), Dataset(), options: new MslRetentionTimeUpdateOptions(AnchorsPerRun: 150));
		Assert.That(result.Runs.Select(r => r.AnchorsSelected), Is.All.EqualTo(150));
		// Fewer anchors span a little less of the run: only peptides beyond the outermost anchors are refused
		Assert.That(result.Runs.Select(r => r.TimesAccepted + r.RefusedExtrapolated), Is.All.EqualTo(PeptideCount));
		Assert.That(result.Filled, Is.GreaterThan(PeptideCount * 0.97));
	}

	// ── Run gates and guards ───────────────────────────────────────────────────

	[Test]
	public void ARunWithTooFewAnchors_ContributesNothing()
	{
		var small = new MslRunObservations("small", Dataset(runCount: 1)[0].Psms.Take(50).ToList(), 100);
		var (result, _, _) = Update(Library(), new List<MslRunObservations> { small });

		Assert.That(result.Runs.Single().Gate, Is.EqualTo("too few anchors"));
		Assert.That(result.EntriesChanged, Is.Zero);
	}

	[Test]
	public void APoorCalibration_ContributesNothing()
	{
		var (result, entries, _) = Update(Library(), Dataset(), new Predictor { ExtraNoise = 40 });

		Assert.That(result.Runs.Select(r => r.Gate), Is.All.EqualTo("residual SD"));
		Assert.That(result.EntriesChanged, Is.Zero);
		Assert.That(entries.Values.Select(e => e.RtStatus), Is.All.EqualTo(None));
	}

	[Test]
	public void NoPredictions_ChangeNothing_AndTheLibraryIsStillWritten()
	{
		string output = NewPath("out");
		MslRetentionTimeUpdateResult result;
		using (MslLibrary lib = MslLibrary.Load(Library()))
			result = lib.UpdateRetentionTimesAndSave(Dataset(), seqs => seqs.Select(_ => double.NaN).ToList(), output);

		Assert.That(result.Runs.Select(r => r.Gate), Is.All.EqualTo("too few anchors"));
		Assert.That(result.PredictionsMissing, Is.EqualTo(result.PredictionsRequested));
		Assert.That(result.EntriesChanged, Is.Zero);
		Assert.That(File.Exists(output), Is.True);
	}

	[Test]
	public void Guards_RefuseDoublePeaksStreaksTheVoidAndExtrapolation()
	{
		const string doublePeak = "DOUBLEPEAKAK", streak = "STREAKINGAK", early = "VOIDEARLYAK", late = "BEYONDENDAK";
		string path = Library(Entry(doublePeak), Entry(streak), Entry(early), Entry(late));
		var runs = Dataset().Select((run, r) =>
		{
			var psms = run.Psms.ToList();
			psms.Add(new MslRunPsm(doublePeak, doublePeak, 2, Minutes(30, r)));
			psms.Add(new MslRunPsm(doublePeak, doublePeak, 2, Minutes(40, r)));      // a second peak 10 iRT later
			foreach (double irt in new[] { 50.0, 52, 54, 56 })
				psms.Add(new MslRunPsm(streak, streak, 2, Minutes(irt, r)));         // 6 iRT of tailing, no gap over 3
			psms.Add(new MslRunPsm(early, early, 2, 2.0));                            // in the first 4% of 100 min
			psms.Add(new MslRunPsm(late, late, 2, 97.0));                             // past the last anchor
			return run with { Psms = psms };
		}).ToList();

		var (result, entries, _) = Update(path, runs);

		Assert.That(result.Runs.Select(r => r.RefusedDoublePeak), Is.All.EqualTo(1));
		Assert.That(result.Runs.Select(r => r.RefusedStreak), Is.All.EqualTo(1));
		Assert.That(result.Runs.Select(r => r.RefusedVoid), Is.All.EqualTo(1));
		Assert.That(result.Runs.Select(r => r.RefusedExtrapolated), Is.All.EqualTo(1));
		foreach (string seq in new[] { doublePeak, streak, early, late })
			Assert.That(entries[seq].RtStatus, Is.EqualTo(None), seq);
		Assert.That(result.Filled, Is.EqualTo(PeptideCount));
	}

	// ── Across datasets ────────────────────────────────────────────────────────

	[Test]
	public void ASecondDatasetThatAgrees_Confirms_WithTheMeanOfTheTwo_AndAThirdKeepsIt()
	{
		string seq = Peptide(17);
		var (_, _, first) = Update(Library(), Dataset());

		var (second, secondEntries, secondPath) = Update(first, Dataset(moved: new() { [seq] = Base[seq] + 4 }));
		Assert.That(second.Confirmed, Is.EqualTo(PeptideCount));
		Assert.That(secondEntries[seq].RtStatus, Is.EqualTo(Confirmed));
		Assert.That(secondEntries[seq].RetentionTime, Is.EqualTo(Base[seq] + 2).Within(1.0), "the mean of the two datasets");

		var (third, thirdEntries, _) = Update(secondPath, Dataset(moved: new() { [seq] = Base[seq] + 6 }));
		Assert.That(third.AlreadyConfirmed, Is.EqualTo(PeptideCount));
		Assert.That(thirdEntries[seq].RetentionTime, Is.EqualTo(secondEntries[seq].RetentionTime).Within(1e-4));
	}

	[Test]
	public void AStrongDisagreement_Challenges_AndASecondReplaces()
	{
		string seq = Peptide(23);
		var (_, _, first) = Update(Library(), Dataset());
		var moved = new Dictionary<string, double> { [seq] = Base[seq] + 20 };

		var (second, secondEntries, secondPath) = Update(first, Dataset(moved: moved));
		Assert.That(second.Challenged, Is.EqualTo(1));
		Assert.That(secondEntries[seq].RtStatus, Is.EqualTo(Challenged));
		Assert.That(secondEntries[seq].RetentionTime, Is.EqualTo(Base[seq]).Within(1.0), "the first strike keeps the time");
		Assert.That(secondEntries[seq].RtStruckByPrediction, Is.False);

		var (third, thirdEntries, _) = Update(secondPath, Dataset(moved: moved));
		Assert.That(third.ReplacedSecondStrike, Is.EqualTo(1));
		Assert.That(thirdEntries[seq].RtStatus, Is.EqualTo(ObservedOnce));
		Assert.That(thirdEntries[seq].RetentionTime, Is.EqualTo(Base[seq] + 20).Within(1.0));
	}

	[Test]
	public void BetweenTheTolerances_NothingChanges()
	{
		string seq = Peptide(29);
		var (_, firstEntries, first) = Update(Library(), Dataset());
		var (result, entries, _) = Update(first, Dataset(moved: new() { [seq] = Base[seq] + 9.5 }));

		Assert.That(result.Inconclusive, Is.EqualTo(1));
		Assert.That(entries[seq].RtStatus, Is.EqualTo(ObservedOnce));
		Assert.That(entries[seq].RetentionTime, Is.EqualTo(firstEntries[seq].RetentionTime).Within(1e-4));
	}

	[TestCase(Challenged)]
	[TestCase(Confirmed)]
	public void AStoredTimeFarFromThePrediction_GivesWayToANewTimeCloseToIt(MslFormat.RetentionTimeStatus status)
	{
		string seq = Peptide(41);
		string path = Library(AgreeingBackground().Append(Entry(seq, rt: Base[seq] + 30, status: status)).ToArray());

		var (result, entries, _) = Update(path, Dataset());

		Assert.That(result.OnlyConfirmations, Is.False);
		Assert.That(result.ReplacedTiebreak, Is.EqualTo(1));
		Assert.That(entries[seq].RtStatus, Is.EqualTo(ObservedOnce));
		Assert.That(entries[seq].RetentionTime, Is.EqualTo(Base[seq]).Within(1.0));
	}

	[Test]
	public void ADatasetThatDisagreesTooOften_OnlyConfirms()
	{
		// Of 200 measured times, every other one is 20 iRT off; the dataset sees the truth
		var stored = Enumerable.Range(0, 200)
			.Select(i => Entry(Peptide(i), rt: TrueIrt(i) + (i % 2 == 0 ? 20 : 0), status: ObservedOnce)).ToArray();
		var (result, entries, _) = Update(Library(stored), Dataset());

		Assert.That(result.OnlyConfirmations, Is.True);
		Assert.That((result.Overlaps, result.OverlapsDisagreeing), Is.EqualTo((200, 100)));
		Assert.That(result.Confirmed, Is.EqualTo(100));
		Assert.That(result.Challenged + result.Filled + result.ReplacedTiebreak, Is.Zero);
		Assert.That(result.HeldBack, Is.EqualTo(PeptideCount - 100));
		Assert.That(entries[Peptide(0)].RtStatus, Is.EqualTo(ObservedOnce), "a disagreeing time is not challenged");
		Assert.That(entries[Peptide(300)].RtStatus, Is.EqualTo(None), "a missing time is not filled");
	}

	// ── Times from a single run ────────────────────────────────────────────────

	[Test]
	public void ASingleRunTime_ThePredictionAgreesWith_IsFilled()
	{
		const string single = "SINGLERUNAGREEK";
		var predictor = new Predictor { Known = { [single] = 44.0 } };
		var (result, entries, _) = Update(Library(Entry(single)), Dataset(onlyInRun0: new() { [single] = 44.0 }), predictor);

		Assert.That(entries[single].RtStatus, Is.EqualTo(ObservedOnce));
		Assert.That(entries[single].RetentionTime, Is.EqualTo(44.0).Within(1.0));
		Assert.That(result.FilledChallenged, Is.Zero);
	}

	[Test]
	public void ASingleRunTime_ThePredictionRefuses_IsKeptButStruck_AndALaterDatasetConfirmsIt()
	{
		const string single = "SINGLERUNFARK";
		var predictor = new Predictor { Known = { [single] = 44.0 } }; // the run sees it 30 iRT later

		var (result, entries, path) = Update(Library(Entry(single)), Dataset(onlyInRun0: new() { [single] = 74.0 }), predictor);
		Assert.That(result.FilledChallenged, Is.EqualTo(1));
		Assert.That(entries[single].RtStatus, Is.EqualTo(Challenged));
		Assert.That(entries[single].RtStruckByPrediction, Is.True);
		Assert.That(entries[single].RetentionTime, Is.EqualTo(74.0).Within(1.0), "the measured time is kept");

		// Another dataset sees it in every run where the first did: confirmed, and the strike cleared
		var (_, secondEntries, _) = Update(path, Dataset(extra: new() { [single] = (single, 74.0) }), predictor);
		Assert.That(secondEntries[single].RtStatus, Is.EqualTo(Confirmed));
		Assert.That(secondEntries[single].RtStruckByPrediction, Is.False);
	}

	[Test]
	public void AStruckTimeReplacedByOneThePredictionAlsoRefuses_IsStruckAgain()
	{
		string seq = Peptide(57);
		string path = Library(AgreeingBackground().Append(Entry(seq, rt: Base[seq] + 25, status: Challenged, struck: true)).ToArray());

		var (result, entries, _) = Update(path, Dataset(moved: new() { [seq] = Base[seq] + 50 }));

		Assert.That(result.ReplacedSecondStrike, Is.EqualTo(1));
		Assert.That(result.ReplacedAndStruckAgain, Is.EqualTo(1));
		Assert.That(entries[seq].RtStatus, Is.EqualTo(Challenged));
		Assert.That(entries[seq].RtStruckByPrediction, Is.True);
		Assert.That(entries[seq].RetentionTime, Is.EqualTo(Base[seq] + 50).Within(1.0));
	}

	[Test]
	public void ASingleRunTimeDisagreeingWithTheStoredTime_MustAgreeWithThePrediction()
	{
		const string seq = "STOREDSINGLEK";
		var predictor = new Predictor { Known = { [seq] = 60.0 } };
		// Stored and confirmed at 60; one run sees it at 100, far from the prediction too
		string path = Library(AgreeingBackground().Append(Entry(seq, rt: 60.0, status: Confirmed)).ToArray());

		var (result, entries, _) = Update(path, Dataset(onlyInRun0: new() { [seq] = 100.0 }), predictor);

		Assert.That(result.RefusedSingleRunTimes, Is.EqualTo(1));
		Assert.That(entries[seq].RtStatus, Is.EqualTo(Confirmed));
		Assert.That(entries[seq].RetentionTime, Is.EqualTo(60.0).Within(1e-4));
	}

	[Test]
	public void AModifiedSingleRunTime_IsNeverStruckByAPrediction()
	{
		const string bas = "PEPTIDEWAYK";
		string full = "PEPT" + Phospho + "IDEWAYK";
		var runs = Dataset();
		runs[0] = runs[0] with { Psms = runs[0].Psms.Append(new MslRunPsm(full, bas, 2, Minutes(80, 0))).ToList() };

		var (result, entries, _) = Update(Library(Entry(full, bas)), runs);

		Assert.That(entries[full].RtStatus, Is.EqualTo(ObservedOnce));
		Assert.That(result.FilledChallenged, Is.Zero);
	}

	[Test]
	public void ALongPeptide_IsCheckedOnlyAgainstTheLongPeptidePredictor()
	{
		string longPeptide = new string('A', 20) + new string('L', 14) + "K"; // 35 residues
		var runs = Dataset(onlyInRun0: new() { [longPeptide] = 90.0 });

		var (_, withoutPredictor, _) = Update(Library(Entry(longPeptide)), runs);
		Assert.That(withoutPredictor[longPeptide].RtStatus, Is.EqualTo(ObservedOnce), "no long-peptide predictor: not checked");

		var asked = new List<string>();
		var (_, withPredictor, _) = Update(Library(Entry(longPeptide)), runs,
			longPredictor: seqs => { asked.AddRange(seqs); return seqs.Select(_ => 50.0).ToList(); });
		Assert.That(asked, Is.EqualTo(new[] { longPeptide }));
		Assert.That(withPredictor[longPeptide].RtStatus, Is.EqualTo(Challenged));
		Assert.That(withPredictor[longPeptide].RtStruckByPrediction, Is.True);
	}

	// ── Scope and contract ─────────────────────────────────────────────────────

	[Test]
	public void PrecursorsNotInTheLibrary_AndDecoys_AreLeftAlone()
	{
		const string decoy = "DECOYSEQUENCEK", unknown = "NOTINLIBRARYK";
		string path = Library(Entry(decoy, isDecoy: true));
		var runs = Dataset(extra: new() { [decoy] = (decoy, 30.0), [unknown] = (unknown, 35.0) });

		var (result, entries, _) = Update(path, runs);

		Assert.That(result.NotInLibrary, Is.EqualTo(2));
		Assert.That(entries[decoy].RtStatus, Is.EqualTo(None));
		Assert.That(entries[decoy].RetentionTime, Is.NaN);
	}

	[Test]
	public void APredictorThatBreaksItsContract_IsReported_AndNothingIsWritten()
	{
		using MslLibrary lib = MslLibrary.Load(Library());

		string wrongCount = NewPath("wrongcount");
		Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), seqs => new double[seqs.Count + 1], wrongCount),
			Throws.InvalidOperationException);
		Assert.That(File.Exists(wrongCount), Is.False);

		string nullResult = NewPath("nullresult");
		Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), _ => null!, nullResult), Throws.InvalidOperationException);
		Assert.That(File.Exists(nullResult), Is.False);

		string throwing = NewPath("throwing");
		Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), _ => throw new TimeoutException(), throwing),
			Throws.TypeOf<TimeoutException>());
		Assert.That(File.Exists(throwing), Is.False);
	}

	[Test]
	public void Arguments_AreChecked()
	{
		MslLibrary lib = MslLibrary.Load(Library());
		var predictor = new Predictor();
		string output = NewPath("args");

		Assert.That(() => lib.UpdateRetentionTimesAndSave(null!, predictor.Predict, output), Throws.ArgumentNullException);
		Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), null!, output), Throws.ArgumentNullException);
		Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), predictor.Predict, null!), Throws.ArgumentNullException);
		Assert.That(() => lib.UpdateRetentionTimesAndSave(new List<MslRunObservations> { null! }, predictor.Predict, output),
			Throws.ArgumentException);
		Assert.That(() => lib.UpdateRetentionTimesAndSave(
				new List<MslRunObservations> { new("r", new List<MslRunPsm> { new(null!, "A", 2, 1) }) }, predictor.Predict, output),
			Throws.ArgumentException);

		foreach (var bad in new[]
		{
			new MslRetentionTimeUpdateOptions(AnchorsPerRun: 0),
			new MslRetentionTimeUpdateOptions(MinimumAnchors: 2),
			new MslRetentionTimeUpdateOptions(MaxResidualSd: double.NaN),
			new MslRetentionTimeUpdateOptions(AgreementTolerance: 0),
			new MslRetentionTimeUpdateOptions(AgreementTolerance: 8, StrongDisagreement: 7),
			new MslRetentionTimeUpdateOptions(PredictionTolerance: -1),
			new MslRetentionTimeUpdateOptions(LongPeptidePredictionTolerance: double.PositiveInfinity),
			new MslRetentionTimeUpdateOptions(DatasetDisagreementLimit: 0),
			new MslRetentionTimeUpdateOptions(DatasetDisagreementLimit: 1.5),
			new MslRetentionTimeUpdateOptions(VoidFraction: 1),
			new MslRetentionTimeUpdateOptions(DoublePeakGap: 0),
			new MslRetentionTimeUpdateOptions(StreakSpan: -3),
		})
			Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), predictor.Predict, output, options: bad),
				Throws.TypeOf<ArgumentOutOfRangeException>(), bad.ToString());
		Assert.That(File.Exists(output), Is.False);

		lib.Dispose();
		Assert.That(() => lib.UpdateRetentionTimesAndSave(Dataset(), predictor.Predict, output), Throws.TypeOf<ObjectDisposedException>());
	}
}
