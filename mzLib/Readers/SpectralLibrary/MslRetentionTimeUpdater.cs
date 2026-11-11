using Chromatography.RetentionTimeCalibration;
using MzLibUtil;
using Omics.SpectralMatch.MslSpectralLibrary;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;

namespace Readers.SpectralLibrary;

/// <summary>
/// The retention-time update behind <see cref="MslLibrary.UpdateRetentionTimesAndSave"/>: each run is
/// calibrated onto the iRT scale against predicted anchors, each run gives at most one guarded time per
/// precursor, the dataset's time is the median over its runs, and the status rules decide what happens to
/// the library's stored time. The rules are those tested on PXReprise's rat, mouse and human library
/// chains (2026-11); see the method's documentation.
/// </summary>
internal static class MslRetentionTimeUpdater
{
	private const string FixedCarbamidomethyl = "Common Fixed:Carbamidomethyl on C";
	private const string OxidationOnM = "Common Variable:Oxidation on M";

	/// <summary>A bracketed modification label, allowing one level of nested brackets.</summary>
	private static readonly Regex ModificationLabel = new(@"\[((?:[^\[\]]|\[[^\]]*\])*)\]", RegexOptions.Compiled);

	/// <summary>A third of the anchored span with fewer anchors than this is unreliable.</summary>
	private const int MinimumAnchorsPerThird = 10;

	/// <summary>A third of the anchored span whose median residual exceeds this many residual SDs is unreliable.</summary>
	private const double MaxThirdBiasInResidualSds = 1.0;

	/// <summary>The tiebreak replaces a stored time farther than this from the prediction ...</summary>
	private const double TiebreakFarFromPrediction = 15.0;

	/// <summary>... with a new time within this of it.</summary>
	private const double TiebreakNearPrediction = 5.0;

	/// <summary>
	/// Added to a predicted iRT to give the expected measured iRT, by base length (7-10, 11-15, 16-20,
	/// 21-25, 26-30 residues): the median of measured minus Prosit 2025 iRT on held-out mouse runs.
	/// </summary>
	private static readonly double[] PredictionOffsetByLength = { -3.54, 0.35, 2.39, 4.90, 4.91 };

	internal static (List<MslLibraryEntry> Entries, MslRetentionTimeUpdateResult Result) Update(
		IReadOnlyList<MslLibraryEntry> entries,
		IReadOnlyList<MslRunObservations> runs,
		Func<IReadOnlyList<string>, IReadOnlyList<double>> predictIrt,
		Func<IReadOnlyList<string>, IReadOnlyList<double>>? predictLongPeptideIrt,
		MslRetentionTimeUpdateOptions options,
		string outputPath)
	{
		var prosit = new Predictions(predictIrt);
		Predictions? longPeptide = predictLongPeptideIrt is null ? null : new Predictions(predictLongPeptideIrt);

		// ── Each run: its precursors and its anchors; one prediction request for every run's anchors ──
		var runPrecursors = runs.Select(Precursors).ToList();
		var runAnchors = runPrecursors.Select(p => SelectAnchors(p, options.AnchorsPerRun)).ToList();
		prosit.Ensure(runAnchors.SelectMany(a => a.Select(x => x.FullSequence)));

		// ── Each run: calibrate, then at most one guarded time per precursor ──
		var datasetTimes = new Dictionary<(string FullSequence, int Charge), (string BaseSequence, List<double> Times)>();
		var reports = new List<MslRunCalibrationReport>(runs.Count);
		for (int r = 0; r < runs.Count; r++)
			reports.Add(CalibrateAndObserve(runs[r], runPrecursors[r], runAnchors[r], prosit, options, datasetTimes));

		var counts = new Counts();
		List<MslLibraryEntry> updated = Merge(entries, datasetTimes, prosit, longPeptide, options, counts);

		var result = new MslRetentionTimeUpdateResult
		{
			OutputPath = outputPath,
			Runs = reports,
			Filled = counts.Filled,
			FilledChallenged = counts.FilledChallenged,
			Confirmed = counts.Confirmed,
			AlreadyConfirmed = counts.AlreadyConfirmed,
			Challenged = counts.Challenged,
			ReplacedSecondStrike = counts.ReplacedSecondStrike,
			ReplacedAndStruckAgain = counts.ReplacedAndStruckAgain,
			ReplacedTiebreak = counts.ReplacedTiebreak,
			Inconclusive = counts.Inconclusive,
			RefusedSingleRunTimes = counts.RefusedSingleRunTimes,
			HeldBack = counts.HeldBack,
			NotInLibrary = counts.NotInLibrary,
			Overlaps = counts.Overlaps,
			OverlapsDisagreeing = counts.OverlapsDisagreeing,
			OnlyConfirmations = counts.OnlyConfirmations,
			PredictionsRequested = prosit.Requested + (longPeptide?.Requested ?? 0),
			PredictionsMissing = prosit.Missing + (longPeptide?.Missing ?? 0),
			EntriesChanged = counts.EntriesChanged,
		};
		return (updated, result);
	}

	// ────────────────────────────────────────────────────────────────────────
	// One run
	// ────────────────────────────────────────────────────────────────────────

	private sealed class RunPrecursor
	{
		public required string FullSequence { get; init; }
		public required string BaseSequence { get; init; }
		public required int Charge { get; init; }
		public List<double> Times { get; } = new();
		public double MedianTime { get; set; }
	}

	/// <summary>The run's precursors (full sequence and charge), each with its sorted identification times.</summary>
	private static List<RunPrecursor> Precursors(MslRunObservations run)
	{
		var byKey = new Dictionary<(string, int), RunPrecursor>();
		foreach (MslRunPsm psm in run.Psms)
		{
			if (!double.IsFinite(psm.RetentionTimeMinutes) || psm.RetentionTimeMinutes < 0)
				continue;
			if (!byKey.TryGetValue((psm.FullSequence, psm.Charge), out RunPrecursor? p))
				byKey[(psm.FullSequence, psm.Charge)] = p = new RunPrecursor
				{
					FullSequence = psm.FullSequence,
					BaseSequence = psm.BaseSequence,
					Charge = psm.Charge
				};
			p.Times.Add(psm.RetentionTimeMinutes);
		}

		var list = byKey.Values.ToList();
		foreach (RunPrecursor p in list)
		{
			p.Times.Sort();
			p.MedianTime = Median(p.Times);
		}
		return list;
	}

	/// <summary>
	/// The run's calibration anchors (rule msl:DEF-RT-CAL v1). Candidates are unmodified peptides of 7-30
	/// residues (fixed carbamidomethyl C allowed), one per base sequence, at the median time of all its
	/// identifications in the run. They are split into 15 strata (5 length bins x 3 equal-width thirds of the
	/// candidates' time range); min(<paramref name="target"/>, candidates) anchors are shared among the strata
	/// in proportion, by largest remainder (ties to the earlier stratum); within a stratum the candidates are
	/// taken in the order of the lowercase hex SHA-256 of their base sequence, so runs that share peptides
	/// pick the same ones, whatever the input order.
	/// </summary>
	private static List<(string FullSequence, double Time)> SelectAnchors(List<RunPrecursor> precursors, int target)
	{
		var byBase = new Dictionary<string, (List<double> Times, string FullSequence)>(StringComparer.Ordinal);
		foreach (RunPrecursor p in precursors)
		{
			int length = p.BaseSequence.Length;
			if (length < 7 || length > 30 || !Describe(p.FullSequence, p.BaseSequence).Unmodified)
				continue;
			if (byBase.TryGetValue(p.BaseSequence, out var e))
			{
				e.Times.AddRange(p.Times);
				if (string.CompareOrdinal(p.FullSequence, e.FullSequence) < 0)
					byBase[p.BaseSequence] = (e.Times, p.FullSequence);
			}
			else
				byBase[p.BaseSequence] = (new List<double>(p.Times), p.FullSequence);
		}

		var candidates = byBase
			.Select(kv => (Base: kv.Key, kv.Value.FullSequence, Time: Median(kv.Value.Times)))
			.ToList();
		if (candidates.Count == 0)
			return new List<(string, double)>();

		double lo = candidates.Min(c => c.Time), hi = candidates.Max(c => c.Time);
		int Third(double t) => hi > lo ? Math.Min(2, (int)Math.Floor(3 * (t - lo) / (hi - lo))) : 0;

		var strata = new List<(string Base, string FullSequence, double Time)>[15];
		for (int i = 0; i < strata.Length; i++)
			strata[i] = new();
		foreach (var c in candidates)
			strata[LengthBin(c.Base.Length) * 3 + Third(c.Time)].Add(c);

		int n = Math.Min(target, candidates.Count);
		double[] exact = strata.Select(s => (double)s.Count * n / candidates.Count).ToArray();
		int[] quota = exact.Select(e => (int)Math.Floor(e)).ToArray();
		int left = n - quota.Sum();
		foreach (int i in Enumerable.Range(0, strata.Length).OrderByDescending(i => exact[i] - quota[i]).ThenBy(i => i).Take(left))
			quota[i]++;

		var picked = new List<(string, double)>(n);
		for (int i = 0; i < strata.Length; i++)
			picked.AddRange(strata[i]
				.OrderBy(c => Sha256Hex(c.Base), StringComparer.Ordinal)
				.Take(quota[i])
				.Select(c => (c.FullSequence, c.Time)));
		return picked;
	}

	private static MslRunCalibrationReport CalibrateAndObserve(
		MslRunObservations run,
		List<RunPrecursor> precursors,
		List<(string FullSequence, double Time)> selected,
		Predictions prosit,
		MslRetentionTimeUpdateOptions o,
		Dictionary<(string, int), (string BaseSequence, List<double> Times)> datasetTimes)
	{
		var anchors = new List<(double Time, double Irt)>();
		foreach (var (full, time) in selected)
		{
			double irt = prosit.Get(full);
			if (double.IsFinite(irt))
				anchors.Add((time, irt));
		}

		MslRunCalibrationReport Gated(string gate, double residualSd) =>
			new(run.RunName, selected.Count, anchors.Count, gate, residualSd, 0, 0, 0, 0, 0, 0);

		if (anchors.Count < o.MinimumAnchors)
			return Gated("too few anchors", double.NaN);

		IrtCalibrationModel model;
		try
		{
			model = IrtCalibration.Fit(anchors.Select(a => (new RtMinutes(a.Time), new Irt(a.Irt))).ToList());
		}
		catch (ArgumentException)
		{
			return Gated("fit failed", double.NaN);
		}
		if (model.AnchorCount < o.MinimumAnchors)
			return Gated("too few anchors", model.ResidualSd);
		if (!(model.ResidualSd <= o.MaxResidualSd))
			return Gated("residual SD", model.ResidualSd);

		double ToIrt(double minutes) => model.ToIrt(new RtMinutes(minutes)).Value;

		// A third of the anchored span with too few anchors, or a biased fit, refuses its times
		double a0 = model.FirstAnchor.Value, a1 = model.LastAnchor.Value;
		int AnchoredThird(double t) => Math.Clamp((int)((t - a0) / Math.Max(1e-9, a1 - a0) * 3), 0, 2);
		var residuals = new[] { new List<double>(), new List<double>(), new List<double>() };
		foreach (var (time, irt) in anchors)
			if (time >= a0 && time <= a1)
				residuals[AnchoredThird(time)].Add(irt - ToIrt(time));
		bool[] unreliable = residuals
			.Select(r => r.Count < MinimumAnchorsPerThird || Math.Abs(Median(r)) > MaxThirdBiasInResidualSds * model.ResidualSd)
			.ToArray();

		double runLength = double.IsFinite(run.RunLengthMinutes) && run.RunLengthMinutes > 0
			? run.RunLengthMinutes
			: precursors.Count > 0 ? precursors.Max(p => p.Times[^1]) : 0.0;
		double voidEnd = o.VoidFraction * runLength;

		int accepted = 0, doublePeak = 0, streak = 0, inVoid = 0, extrapolated = 0, unreliableStretch = 0;
		foreach (RunPrecursor p in precursors)
		{
			if (p.Times.Count >= 2)
			{
				double gap = 0;
				for (int i = 1; i < p.Times.Count; i++)
					gap = Math.Max(gap, ToIrt(p.Times[i]) - ToIrt(p.Times[i - 1]));
				if (gap > o.DoublePeakGap) { doublePeak++; continue; }
				if (ToIrt(p.Times[^1]) - ToIrt(p.Times[0]) > o.StreakSpan) { streak++; continue; }
			}

			double t = p.MedianTime;
			if (t < voidEnd) { inVoid++; continue; }
			if (model.IsExtrapolated(new RtMinutes(t))) { extrapolated++; continue; }
			if (unreliable[AnchoredThird(t)]) { unreliableStretch++; continue; }

			var key = (p.FullSequence, p.Charge);
			if (!datasetTimes.TryGetValue(key, out var d))
				datasetTimes[key] = d = (p.BaseSequence, new List<double>());
			d.Times.Add(ToIrt(t));
			accepted++;
		}

		return new MslRunCalibrationReport(run.RunName, selected.Count, anchors.Count, null, model.ResidualSd,
			accepted, doublePeak, streak, inVoid, extrapolated, unreliableStretch);
	}

	// ────────────────────────────────────────────────────────────────────────
	// The dataset against the library
	// ────────────────────────────────────────────────────────────────────────

	private sealed class Counts
	{
		public int Filled, FilledChallenged, Confirmed, AlreadyConfirmed, Challenged, ReplacedSecondStrike,
			ReplacedAndStruckAgain, ReplacedTiebreak, Inconclusive, RefusedSingleRunTimes, HeldBack, NotInLibrary,
			Overlaps, OverlapsDisagreeing, EntriesChanged;
		public bool OnlyConfirmations;
	}

	private static List<MslLibraryEntry> Merge(
		IReadOnlyList<MslLibraryEntry> entries,
		Dictionary<(string FullSequence, int Charge), (string BaseSequence, List<double> Times)> datasetTimes,
		Predictions prosit,
		Predictions? longPeptide,
		MslRetentionTimeUpdateOptions o,
		Counts counts)
	{
		var index = new Dictionary<(string, int), int>();
		for (int i = 0; i < entries.Count; i++)
			if (!entries[i].IsDecoy)
				index.TryAdd((entries[i].FullSequence, entries[i].ChargeState), i);

		// A stored time counts only when its status says it was measured
		static bool Measured(MslLibraryEntry e) =>
			e.RtStatus != MslFormat.RetentionTimeStatus.None && double.IsFinite(e.RetentionTime);

		// The dataset's time per precursor in the library, in a fixed order
		var work = new List<(int Index, string FullSequence, SequenceInfo Info, double Time, int Runs)>();
		foreach (var ((full, charge), (baseSequence, times)) in datasetTimes
			.OrderBy(kv => kv.Key.FullSequence, StringComparer.Ordinal).ThenBy(kv => kv.Key.Charge))
		{
			if (!index.TryGetValue((full, charge), out int i)) { counts.NotInLibrary++; continue; }
			work.Add((i, full, Describe(full, baseSequence), Median(times), times.Count));
		}

		// One prediction request for the checks below: single-run times to check, stored times that may
		// need the tiebreak; and one long-peptide request
		var need = new List<string>();
		var needLong = new List<string>();
		foreach (var (i, full, info, x, nRuns) in work)
		{
			MslLibraryEntry e = entries[i];
			bool measured = Measured(e);
			if (info.PrositClass)
			{
				bool singleToCheck = nRuns == 1 && !(measured && Math.Abs(x - e.RetentionTime) <= o.AgreementTolerance);
				bool storedToTiebreak = measured &&
					(e.RtStatus == MslFormat.RetentionTimeStatus.Challenged || Math.Abs(x - e.RetentionTime) > o.StrongDisagreement);
				if (singleToCheck || storedToTiebreak)
					need.Add(full);
			}
			if (longPeptide is not null && info.Unmodified && info.Length > 30 && info.Length <= 50 &&
				((nRuns == 1 && !measured) || (measured && e.RtStruckByPrediction)))
				needLong.Add(full);
		}
		prosit.Ensure(need);
		longPeptide?.Ensure(needLong);

		// A dataset that disagrees with too many measured library times only confirms
		foreach (var (i, _, _, x, _) in work)
		{
			MslLibraryEntry e = entries[i];
			if (!Measured(e)) continue;
			counts.Overlaps++;
			if (Math.Abs(x - e.RetentionTime) > o.AgreementTolerance) counts.OverlapsDisagreeing++;
		}
		counts.OnlyConfirmations = counts.Overlaps > 0 && counts.OverlapsDisagreeing > o.DatasetDisagreementLimit * counts.Overlaps;

		// The prediction's verdict on a measured time: Prosit (length-corrected) for unmodified 7-30 aa or
		// Oxidation-on-M-only peptides; the long-peptide predictor for unmodified ones over 30 aa; any other
		// peptide, or no prediction, is never struck
		bool PredictionAgrees(string full, SequenceInfo info, double x)
		{
			if (info.PrositClass)
			{
				double p = prosit.Get(full);
				return !double.IsFinite(p) || Math.Abs(x - (p + PredictionOffsetByLength[LengthBin(info.Length)])) <= o.PredictionTolerance;
			}
			if (info.Unmodified && info.Length > 30)
			{
				double c = longPeptide?.Get(full) ?? double.NaN;
				return !double.IsFinite(c) || Math.Abs(x - c) <= o.LongPeptidePredictionTolerance;
			}
			return true;
		}

		// The tiebreak: a stored time far from the (raw) prediction gives way to a new one close to it
		bool TiebreakReplaces(string full, SequenceInfo info, double stored, double x)
		{
			double p = info.Unmodified ? prosit.Get(full) : double.NaN;
			return double.IsFinite(p) && Math.Abs(stored - p) > TiebreakFarFromPrediction && Math.Abs(x - p) <= TiebreakNearPrediction;
		}

		var result = entries.ToList();
		void Set(int i, double time, MslFormat.RetentionTimeStatus status, bool struck)
		{
			MslLibraryEntry copy = entries[i].WithFragments(entries[i].MatchedFragmentIons);
			copy.RetentionTime = time;
			copy.RtStatus = status;
			copy.RtStruckByPrediction = struck;
			copy.RtIsCalibrated = false; // the update writes iRT
			if (copy.Source == MslFormat.SourceType.Predicted)
				copy.Source = MslFormat.SourceType.EmpiricalRefined;
			result[i] = copy;
			counts.EntriesChanged++;
		}

		foreach (var (i, full, info, x, nRuns) in work)
		{
			MslLibraryEntry e = entries[i];
			bool measured = Measured(e);
			MslFormat.RetentionTimeStatus status = measured ? e.RtStatus : MslFormat.RetentionTimeStatus.None;
			double stored = e.RetentionTime;
			bool agrees = measured && Math.Abs(x - stored) <= o.AgreementTolerance;

			if (counts.OnlyConfirmations && !agrees) { counts.HeldBack++; continue; }

			if (nRuns == 1 && !measured)
			{
				// A time from a single run: a disagreeing prediction is its first strike, but the measured time is kept
				if (!PredictionAgrees(full, info, x))
				{
					Set(i, x, MslFormat.RetentionTimeStatus.Challenged, struck: true);
					counts.FilledChallenged++;
					continue;
				}
			}
			else if (measured && e.RtStruckByPrediction)
			{
				// A later dataset's time for a prediction-struck entry is merged whatever the prediction says
			}
			else if (!agrees && nRuns == 1 && info.PrositClass)
			{
				// A single-run time that disagrees with the stored time must at least agree with the prediction
				double p = prosit.Get(full);
				if (double.IsFinite(p) && Math.Abs(x - (p + PredictionOffsetByLength[LengthBin(info.Length)])) > o.PredictionTolerance)
				{
					counts.RefusedSingleRunTimes++;
					continue;
				}
			}

			if (!measured)
			{
				Set(i, x, MslFormat.RetentionTimeStatus.ObservedOnce, struck: false);
				counts.Filled++;
				continue;
			}

			if (status == MslFormat.RetentionTimeStatus.Challenged && TiebreakReplaces(full, info, stored, x))
			{
				Set(i, x, MslFormat.RetentionTimeStatus.ObservedOnce, struck: false);
				counts.ReplacedTiebreak++;
				continue;
			}

			double d = Math.Abs(x - stored);
			if (d <= o.AgreementTolerance)
			{
				// The first confirmation takes the mean of the two times; a confirmed time is kept
				if (status != MslFormat.RetentionTimeStatus.Confirmed)
				{
					Set(i, 0.5 * (stored + x), MslFormat.RetentionTimeStatus.Confirmed, struck: false);
					counts.Confirmed++;
				}
				else
					counts.AlreadyConfirmed++;
			}
			else if (d > o.StrongDisagreement)
			{
				if (TiebreakReplaces(full, info, stored, x))
				{
					Set(i, x, MslFormat.RetentionTimeStatus.ObservedOnce, struck: false);
					counts.ReplacedTiebreak++;
				}
				else if (status == MslFormat.RetentionTimeStatus.Challenged)
				{
					// The second strike replaces the time; a prediction-struck entry's replacement is checked again
					bool struckAgain = e.RtStruckByPrediction && !PredictionAgrees(full, info, x);
					Set(i, x, struckAgain ? MslFormat.RetentionTimeStatus.Challenged : MslFormat.RetentionTimeStatus.ObservedOnce, struckAgain);
					counts.ReplacedSecondStrike++;
					if (struckAgain) counts.ReplacedAndStruckAgain++;
				}
				else
				{
					// The first strike: the time is kept, marked challenged
					Set(i, stored, MslFormat.RetentionTimeStatus.Challenged, struck: false);
					counts.Challenged++;
				}
			}
			else
				counts.Inconclusive++;
		}

		return result;
	}

	// ────────────────────────────────────────────────────────────────────────
	// Helpers
	// ────────────────────────────────────────────────────────────────────────

	private readonly record struct SequenceInfo(int Length, bool Unmodified, bool OxidationOnly)
	{
		/// <summary>Checked against Prosit: unmodified 7-30 residues, or with Oxidation on M only.</summary>
		public bool PrositClass => Length >= 7 && Length <= 30 && (Unmodified || OxidationOnly);
	}

	/// <summary>Length is the base sequence's; fixed carbamidomethyl C counts as unmodified.</summary>
	private static SequenceInfo Describe(string fullSequence, string baseSequence)
	{
		var labels = ModificationLabel.Matches(fullSequence)
			.Select(m => m.Groups[1].Value)
			.Where(l => l != FixedCarbamidomethyl)
			.Distinct(StringComparer.Ordinal)
			.ToList();
		return new SequenceInfo(baseSequence.Length, labels.Count == 0, labels.Count > 0 && labels.All(l => l == OxidationOnM));
	}

	private static int LengthBin(int length) => length <= 10 ? 0 : length <= 15 ? 1 : length <= 20 ? 2 : length <= 25 ? 3 : 4;

	private static string Sha256Hex(string s) =>
		Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(s))).ToLowerInvariant();

	private static double Median(List<double> values)
	{
		if (values.Count == 0) return double.NaN;
		var a = values.OrderBy(v => v).ToArray();
		return a.Length % 2 == 1 ? a[a.Length / 2] : 0.5 * (a[a.Length / 2 - 1] + a[a.Length / 2]);
	}

	/// <summary>Predictions by full sequence, each requested once.</summary>
	private sealed class Predictions
	{
		private readonly Func<IReadOnlyList<string>, IReadOnlyList<double>> _predict;
		private readonly Dictionary<string, double> _values = new(StringComparer.Ordinal);

		public Predictions(Func<IReadOnlyList<string>, IReadOnlyList<double>> predict) => _predict = predict;

		public int Requested { get; private set; }
		public int Missing { get; private set; }

		public void Ensure(IEnumerable<string> fullSequences)
		{
			var need = fullSequences.Where(s => !_values.ContainsKey(s)).Distinct(StringComparer.Ordinal).ToList();
			if (need.Count == 0)
				return;

			IReadOnlyList<double> values = _predict(need)
				?? throw new InvalidOperationException("The retention-time predictor returned null.");
			if (values.Count != need.Count)
				throw new InvalidOperationException(
					$"The retention-time predictor returned {values.Count} values for {need.Count} sequences.");

			for (int i = 0; i < need.Count; i++)
			{
				double v = double.IsFinite(values[i]) ? values[i] : double.NaN;
				_values[need[i]] = v;
				Requested++;
				if (double.IsNaN(v)) Missing++;
			}
		}

		public double Get(string fullSequence) =>
			_values.TryGetValue(fullSequence, out double v) ? v : double.NaN;
	}
}
