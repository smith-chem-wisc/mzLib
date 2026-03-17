using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary;

// ─────────────────────────────────────────────────────────────────────────────
// DeconvolutedPeak
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// A single deconvoluted fragment peak from a top-down MS2 spectrum.
/// The neutral mass is computed from the observed m/z and charge state by the caller's
/// deconvolution algorithm; this struct is charge-state-agnostic.
/// </summary>
public readonly struct DeconvolutedPeak
{
	/// <summary>Neutral monoisotopic mass (Da).</summary>
	public readonly double NeutralMass;

	/// <summary>Relative intensity (0–1, normalized by caller).</summary>
	public readonly float Intensity;

	/// <param name="neutralMass">Neutral monoisotopic mass in daltons.</param>
	/// <param name="intensity">Relative intensity, normalized 0–1 by the caller.</param>
	public DeconvolutedPeak(double neutralMass, float intensity)
	{
		NeutralMass = neutralMass;
		Intensity = intensity;
	}
}

// ─────────────────────────────────────────────────────────────────────────────
// MslProteoformScoringResult
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Scoring result for one proteoform library candidate against an experimental spectrum.
/// </summary>
public record MslProteoformScoringResult
{
	/// <summary>The library entry that was scored.</summary>
	public required MslLibraryEntry LibraryEntry { get; init; }

	/// <summary>
	/// Composite score = SpectralAngle × log2(1 + MatchedFragmentCount).
	/// Higher is better. Range: 0 to ~(1 × log2(1 + MaxFragments)).
	/// </summary>
	public required double CompositeScore { get; init; }

	/// <summary>Spectral angle component (cosine similarity of intensity vectors). Range 0–1.</summary>
	public required double SpectralAngle { get; init; }

	/// <summary>Number of library fragment ions matched to experimental peaks.</summary>
	public required int MatchedFragmentCount { get; init; }

	/// <summary>Total number of fragment ions in the library entry.</summary>
	public required int TotalLibraryFragments { get; init; }

	/// <summary>
	/// Fraction of library fragments matched (MatchedFragmentCount / TotalLibraryFragments).
	/// Useful as a secondary ranking dimension.
	/// </summary>
	public double FragmentCoverage =>
		TotalLibraryFragments > 0
			? (double)MatchedFragmentCount / TotalLibraryFragments
			: 0.0;

	/// <summary>
	/// The matched fragment ions, each carrying both the library and experimental m/z and
	/// intensity. Useful for output annotation and visualization.
	/// </summary>
	public required IReadOnlyList<MatchedFragmentIon> MatchedFragments { get; init; }
}

// ─────────────────────────────────────────────────────────────────────────────
// MslProteoformScorer
// ─────────────────────────────────────────────────────────────────────────────

/// <summary>
/// Scores top-down proteoform library candidates against a deconvoluted experimental
/// MS2 spectrum. Uses a composite score of spectral angle × log2(1 + matched fragments).
///
/// <para>
/// The composite score rewards both spectral similarity <em>and</em> breadth of fragment
/// coverage. A library entry covering 30% of the proteoform at high confidence will
/// outscore one covering 5% at the same spectral angle — unlike spectral angle alone,
/// which would treat them identically. The log2 scale means: 1 match → 1.0, 3 → 2.0,
/// 7 → 3.0, 15 → 4.0.
/// </para>
///
/// <para>
/// Fragment matching is performed in neutral-mass space (charge-state-independent) using
/// binary search on the sorted experimental peak list, giving O(F log N) complexity per
/// candidate where F = library fragment count and N = experimental peak count.
/// </para>
///
/// Typical usage from MetaMorpheus:
/// <code>
/// // 1. Deconvolute the MS2 scan once (outside the scoring loop)
/// var deconvPeaks = DeconvoluteScan(scan)
///     .Select(p => new DeconvolutedPeak(p.NeutralMass, (float)p.NormalizedIntensity))
///     .ToList();
///
/// // 2. Compute precursor neutral mass from scan header
/// double expNeutralMass =
///     (scan.SelectedIonMZ.Value * scan.SelectedIonChargeStateGuess.Value)
///     - (scan.SelectedIonChargeStateGuess.Value * 1.007276);
///
/// // 3. Query the proteoform index
/// var candidates = library.QueryProteoformMassWindow(
///     expNeutralMass - massTolerance,
///     expNeutralMass + massTolerance);
///
/// // 4. Score all candidates
/// var results = MslProteoformScorer.Score(
///     candidates, library, deconvPeaks,
///     fragmentMassTolerance: 20e-6,
///     minMatchedFragments: 3);
///
/// // 5. Use top result
/// var topHit = results.FirstOrDefault();
/// </code>
/// </summary>
public static class MslProteoformScorer
{
	/// <summary>Proton mass in daltons used for neutral mass computation.</summary>
	private const double ProtonMass = 1.007276;

	// ── Public API ────────────────────────────────────────────────────────────

	/// <summary>
	/// Scores a list of proteoform candidates against a deconvoluted spectrum.
	/// </summary>
	/// <param name="candidates">
	///   Index entries from <see cref="MslProteoformIndex.QueryMassWindow"/>.
	///   An empty span returns an empty list without error.
	/// </param>
	/// <param name="library">
	///   The library to retrieve full entries from. Must not be null.
	/// </param>
	/// <param name="deconvolutedPeaks">
	///   Deconvoluted peaks from the experimental MS2 spectrum. Must not be null.
	///   Peaks need not be pre-sorted; an internal sorted copy is always made.
	///   An empty list returns an empty result list.
	/// </param>
	/// <param name="fragmentMassTolerance">
	///   Fragment matching mass tolerance as a fraction (ppm × 1e-6).
	///   Default: 20e-6 (20 ppm). Wider than typical bottom-up search to accommodate
	///   deconvolution uncertainty on long top-down fragments.
	/// </param>
	/// <param name="minMatchedFragments">
	///   Minimum number of matched fragments required to include a candidate in results.
	///   Default: 3. Candidates below this threshold are excluded entirely.
	/// </param>
	/// <returns>
	///   List of <see cref="MslProteoformScoringResult"/> sorted by
	///   <see cref="MslProteoformScoringResult.CompositeScore"/> descending.
	///   Empty when no candidates pass the minimum-matched-fragment filter.
	/// </returns>
	/// <exception cref="ArgumentNullException">
	///   <paramref name="library"/> or <paramref name="deconvolutedPeaks"/> is null.
	/// </exception>
	public static List<MslProteoformScoringResult> Score(
		ReadOnlySpan<MslProteoformIndexEntry> candidates,
		MslLibrary library,
		IReadOnlyList<DeconvolutedPeak> deconvolutedPeaks,
		double fragmentMassTolerance = 20e-6,
		int minMatchedFragments = 3)
	{
		ArgumentNullException.ThrowIfNull(library);
		ArgumentNullException.ThrowIfNull(deconvolutedPeaks);

		if (candidates.IsEmpty || deconvolutedPeaks.Count == 0)
			return new List<MslProteoformScoringResult>();

		// Sort experimental peaks by neutral mass once; all candidate scoring reuses this.
		DeconvolutedPeak[] sortedPeaks = SortPeaks(deconvolutedPeaks);

		var results = new List<MslProteoformScoringResult>(candidates.Length);

		for (int i = 0; i < candidates.Length; i++)
		{
			MslLibraryEntry? entry = library.GetEntry(candidates[i].OrdinalIndex);
			if (entry is null || entry.Fragments is null || entry.Fragments.Count == 0)
				continue;

			MslProteoformScoringResult? result = ScoreCandidate(
				entry, sortedPeaks, fragmentMassTolerance, minMatchedFragments);

			if (result is not null)
				results.Add(result);
		}

		// Sort descending by composite score.
		results.Sort(static (a, b) => b.CompositeScore.CompareTo(a.CompositeScore));
		return results;
	}

	/// <summary>
	/// Convenience overload that deconvolutes the scan internally before scoring.
	/// Uses mzLib's <see cref="ClassicDeconvolutionParameters"/> with default settings.
	/// Prefer the <see cref="Score(ReadOnlySpan{MslProteoformIndexEntry}, MslLibrary, IReadOnlyList{DeconvolutedPeak}, double, int)"/>
	/// overload if the caller has already deconvoluted, to avoid redundant work.
	/// </summary>
	/// <param name="candidates">Index entries from <see cref="MslProteoformIndex.QueryMassWindow"/>.</param>
	/// <param name="library">Library to retrieve full entries from. Must not be null.</param>
	/// <param name="ms2Scan">The raw MS2 scan to deconvolute internally. Must not be null.</param>
	/// <param name="fragmentMassTolerance">Fragment matching tolerance as a fraction (ppm × 1e-6). Default: 20e-6.</param>
	/// <param name="minMatchedFragments">Minimum matched fragments to include a candidate. Default: 3.</param>
	/// <returns>
	///   List of <see cref="MslProteoformScoringResult"/> sorted by composite score descending.
	/// </returns>
	/// <exception cref="ArgumentNullException">
	///   <paramref name="library"/> or <paramref name="ms2Scan"/> is null.
	/// </exception>
	public static List<MslProteoformScoringResult> Score(
		ReadOnlySpan<MslProteoformIndexEntry> candidates,
		MslLibrary library,
		MsDataScan ms2Scan,
		double fragmentMassTolerance = 20e-6,
		int minMatchedFragments = 3)
	{
		ArgumentNullException.ThrowIfNull(library);
		ArgumentNullException.ThrowIfNull(ms2Scan);

		if (candidates.IsEmpty)
			return new List<MslProteoformScoringResult>();

		List<DeconvolutedPeak> deconvPeaks = DeconvoluteScanInternal(ms2Scan);
		return Score(candidates, library, deconvPeaks, fragmentMassTolerance, minMatchedFragments);
	}

	// ── Private: per-candidate scoring ───────────────────────────────────────

	/// <summary>
	/// Scores one library entry against the sorted experimental peaks.
	/// Returns null when the candidate does not pass <paramref name="minMatchedFragments"/>.
	/// </summary>
	private static MslProteoformScoringResult? ScoreCandidate(
		MslLibraryEntry entry,
		DeconvolutedPeak[] sortedPeaks,
		double fragmentMassTolerance,
		int minMatchedFragments)
	{
		var matchedLibIntensities = new List<float>(entry.Fragments.Count);
		var matchedExpIntensities = new List<float>(entry.Fragments.Count);
		var matchedIons = new List<MatchedFragmentIon>(entry.Fragments.Count);

		foreach (MslFragmentIon libFrag in entry.Fragments)
		{
			// Guard: Charge = 0 produces neutral mass = 0.0, which is incorrect and
			// would silently fail to match any real experimental peak. MslWriter.ValidateEntries
			// flags Charge = 0 as an error, but that check is not called automatically.
			if (libFrag.Charge <= 0) continue;

			// Step 1: Convert library fragment m/z + charge to neutral mass.
			double fragNeutralMass = ((double)libFrag.Mz * libFrag.Charge)
									 - (libFrag.Charge * ProtonMass);

			// Step 2: Absolute tolerance in Da = ppm-fraction × mass.
			double absTol = fragmentMassTolerance * fragNeutralMass;

			// Step 3: Binary search for a matching experimental peak.
			int matchIdx = BinarySearchPeak(sortedPeaks, fragNeutralMass, absTol);
			if (matchIdx < 0)
				continue;

			float expIntensity = sortedPeaks[matchIdx].Intensity;
			float libIntensity = libFrag.Intensity;

			matchedLibIntensities.Add(libIntensity);
			matchedExpIntensities.Add(expIntensity);

			// Step 3 (annotation): build a MatchedFragmentIon for output.
			FragmentationTerminus terminus = libFrag.ProductType switch
			{
				ProductType.b => FragmentationTerminus.N,
				ProductType.c => FragmentationTerminus.N,
				ProductType.a => FragmentationTerminus.N,
				ProductType.y => FragmentationTerminus.C,
				ProductType.z => FragmentationTerminus.C,
				ProductType.w => FragmentationTerminus.C,
				_ => FragmentationTerminus.None
			};

			var product = new Product(
				libFrag.ProductType,
				terminus,
				neutralMass: 0.0,
				libFrag.FragmentNumber,
				libFrag.ResiduePosition,
				libFrag.NeutralLoss);

			matchedIons.Add(new MatchedFragmentIon(
				product,
				libFrag.Mz,
				libFrag.Intensity,
				libFrag.Charge));
		}

		int matchedCount = matchedLibIntensities.Count;

		// Step 6: Exclude candidates below the minimum match threshold.
		if (matchedCount < minMatchedFragments)
			return null;

		// Steps 4–5: Compute spectral angle then composite score.
		double spectralAngle = ComputeSpectralAngle(matchedLibIntensities, matchedExpIntensities);
		double compositeScore = spectralAngle * Math.Log2(1.0 + matchedCount);

		return new MslProteoformScoringResult
		{
			LibraryEntry = entry,
			CompositeScore = compositeScore,
			SpectralAngle = spectralAngle,
			MatchedFragmentCount = matchedCount,
			TotalLibraryFragments = entry.Fragments.Count,
			MatchedFragments = matchedIons.AsReadOnly()
		};
	}

	// ── Private: spectral angle ───────────────────────────────────────────────

	/// <summary>
	/// Computes the arccos-scaled spectral angle (SA) between two equal-length
	/// matched intensity vectors.
	///
	/// <para>Formula:  SA = 1 − (2/π) × arccos( dot(lib,exp) / (|lib| × |exp|) )</para>
	///
	/// <para>Range: 0–1 where 1.0 = perfect spectral match and 0.0 = orthogonal spectra.
	/// This matches the SA metric used throughout MetaMorpheus for bottom-up PSM scoring
	/// (see Bittremieux et al., J. Proteome Res. 2022), ensuring that a single threshold
	/// (e.g. SA > 0.7) applies consistently across both bottom-up and top-down workflows.</para>
	///
	/// <para>The cosine is clamped to [−1, 1] before arccos to guard against floating-point
	/// values microscopically outside the valid domain (e.g. 1.0000000000000002).</para>
	///
	/// <para>Returns 0.0 when either vector has zero magnitude.</para>
	/// </summary>
	private static double ComputeSpectralAngle(
		List<float> libIntensities,
		List<float> expIntensities)
	{
		double dot = 0.0;
		double normLib = 0.0;
		double normExp = 0.0;

		for (int i = 0; i < libIntensities.Count; i++)
		{
			double l = libIntensities[i];
			double e = expIntensities[i];
			dot += l * e;
			normLib += l * l;
			normExp += e * e;
		}

		double denom = Math.Sqrt(normLib) * Math.Sqrt(normExp);
		if (denom <= 0.0) return 0.0;

		// Clamp before arccos: floating-point arithmetic can produce values like
		// 1.0000000000000002 for identical vectors, which arccos would return NaN for.
		double cosine = Math.Clamp(dot / denom, -1.0, 1.0);

		return 1.0 - (2.0 / Math.PI) * Math.Acos(cosine);
	}

	// ── Private: binary search on sorted peaks ────────────────────────────────

	/// <summary>
	/// Locates the index of the <see cref="DeconvolutedPeak"/> whose
	/// <see cref="DeconvolutedPeak.NeutralMass"/> is within
	/// [<paramref name="targetMass"/> − <paramref name="absTol"/>,
	///  <paramref name="targetMass"/> + <paramref name="absTol"/>] and
	/// closest to <paramref name="targetMass"/>.
	/// Returns -1 when no peak is within tolerance.
	/// </summary>
	private static int BinarySearchPeak(
		DeconvolutedPeak[] sortedPeaks,
		double targetMass,
		double absTol)
	{
		double lo = targetMass - absTol;
		double hi = targetMass + absTol;

		// Lower-bound binary search for lo.
		int left = 0;
		int right = sortedPeaks.Length;

		while (left < right)
		{
			int mid = left + ((right - left) >> 1);
			if (sortedPeaks[mid].NeutralMass < lo)
				left = mid + 1;
			else
				right = mid;
		}

		// Scan forward to find the closest peak within the [lo, hi] window.
		int bestIdx = -1;
		double bestDist = double.MaxValue;

		for (int i = left; i < sortedPeaks.Length && sortedPeaks[i].NeutralMass <= hi; i++)
		{
			double dist = Math.Abs(sortedPeaks[i].NeutralMass - targetMass);
			if (dist < bestDist)
			{
				bestDist = dist;
				bestIdx = i;
			}
		}

		return bestIdx;
	}

	// ── Private: peak sorting ─────────────────────────────────────────────────

	/// <summary>
	/// Returns a copy of the experimental peak list sorted by neutral mass ascending.
	/// A copy is always made to provide a stable backing array for binary search across
	/// all candidates in the same Score call.
	/// </summary>
	private static DeconvolutedPeak[] SortPeaks(IReadOnlyList<DeconvolutedPeak> peaks)
	{
		var arr = new DeconvolutedPeak[peaks.Count];
		for (int i = 0; i < peaks.Count; i++)
			arr[i] = peaks[i];

		Array.Sort(arr, static (a, b) => a.NeutralMass.CompareTo(b.NeutralMass));
		return arr;
	}

	// ── Private: internal deconvolution ──────────────────────────────────────

	/// <summary>
	/// Deconvolutes an MS2 scan using mzLib's <see cref="ClassicDeconvolutionParameters"/>
	/// and returns the result as a list of <see cref="DeconvolutedPeak"/> values with
	/// intensities normalized so the most abundant peak = 1.0.
	/// </summary>
	private static List<DeconvolutedPeak> DeconvoluteScanInternal(MsDataScan scan)
	{
		// ClassicDeconvolutionParameters(minAssumedChargeState, maxAssumedChargeState,
		//                                deconvolutionTolerancePpm, intensityRatioLimit)
		var deconvParams = new ClassicDeconvolutionParameters(
			minCharge: 1,
			maxCharge: 60,
			deconPpm: 4,
			intensityRatio: 3);

		// Use Deconvoluter.Deconvolute(scan, parameters) — the correct mzLib entry point.
		var envelopes = Deconvoluter.Deconvolute(scan, deconvParams)
			.ToList();

		if (envelopes.Count == 0)
			return new List<DeconvolutedPeak>();

		double maxIntensity = envelopes.Max(e => e.TotalIntensity);
		if (maxIntensity <= 0.0)
			maxIntensity = 1.0;

		var peaks = new List<DeconvolutedPeak>(envelopes.Count);
		foreach (var env in envelopes)
		{
			peaks.Add(new DeconvolutedPeak(
				neutralMass: env.MonoisotopicMass,
				intensity: (float)(env.TotalIntensity / maxIntensity)));
		}

		return peaks;
	}
}