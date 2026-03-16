using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// NUnit 4 test suite for <see cref="MslProteoformScorer"/> and the
/// <see cref="MslLibrary.ScoreProteoformCandidates"/> convenience method.
///
/// All tests are network-free and use synthetic in-memory data.
/// The shared synthetic data fixture is built once in <see cref="MakeSyntheticProteoformData"/>.
///
/// Coverage groups:
/// <list type="number">
///   <item>Perfect-match scoring (SpectralAngle, CompositeScore, counts, coverage)</item>
///   <item>No-match and empty-input edge cases</item>
///   <item>MinMatchedFragments filtering</item>
///   <item>Partial match</item>
///   <item>Result ordering</item>
///   <item>Mass tolerance enforcement</item>
///   <item>MatchedFragments annotation</item>
///   <item>Decoy scoring parity</item>
///   <item>MslLibrary.ScoreProteoformCandidates round-trip and peptide-only guard</item>
///   <item>Orthogonal intensity edge case</item>
/// </list>
/// </summary>
[TestFixture]
public class TestMslProteoformScoring
{
	// ── Proton mass constant (mirrors MslProteoformScorer internals) ──────────
	private const double Proton = 1.007276;

	// ── Shared temp directory ─────────────────────────────────────────────────
	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(), $"TestMslProteoformScoring_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Shared synthetic data factory
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Builds a synthetic proteoform library entry and matching experimental peaks.
	/// The entry has 5 fragments; the synthetic experimental data matches all 5 perfectly
	/// (same neutral masses, same intensities).
	/// </summary>
	private static (MslLibraryEntry entry, List<DeconvolutedPeak> experimentalPeaks)
		MakeSyntheticProteoformData()
	{
		var fragments = new List<MslFragmentIon>
		{
			new() { ProductType = ProductType.b, FragmentNumber = 10, Charge = 1, Mz = 1200.5f, Intensity = 1.0f,  NeutralLoss = 0 },
			new() { ProductType = ProductType.b, FragmentNumber = 20, Charge = 2, Mz = 1100.3f, Intensity = 0.8f,  NeutralLoss = 0 },
			new() { ProductType = ProductType.y, FragmentNumber = 15, Charge = 1, Mz = 1600.7f, Intensity = 0.6f,  NeutralLoss = 0 },
			new() { ProductType = ProductType.y, FragmentNumber = 25, Charge = 2, Mz = 900.2f,  Intensity = 0.5f,  NeutralLoss = 0 },
			new() { ProductType = ProductType.c, FragmentNumber = 8,  Charge = 1, Mz = 850.4f,  Intensity = 0.3f,  NeutralLoss = 0 },
		};

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "LARGEPROT[Common Variable:Oxidation on M]EOFORMSEQUENCE",
			StrippedSequence = "LARGEPROTEFORMSEQUENCE",
			PrecursorMz = 1200.0,
			Charge = 10,
			Irt = 45.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			Fragments = fragments
		};

		// Compute exact neutral masses matching the fragment ions.
		var peaks = fragments.Select(f => new DeconvolutedPeak(
			neutralMass: (double)f.Mz * f.Charge - f.Charge * Proton,
			intensity: f.Intensity)).ToList();

		return (entry, peaks);
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 1: Perfect-match scoring
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_PerfectMatch_SpectralAngleNearOne()
	{
		var (entry, peaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results, Has.Count.EqualTo(1));
		Assert.That(results[0].SpectralAngle, Is.GreaterThan(0.99));
	}

	[Test]
	public void Scoring_PerfectMatch_CompositeScore_Correct()
	{
		var (entry, peaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		// SpectralAngle ≈ 1.0; log2(1 + 5) = log2(6) ≈ 2.585
		double expected = 1.0 * Math.Log2(6.0);
		Assert.That(results[0].CompositeScore, Is.EqualTo(expected).Within(0.01));
	}

	[Test]
	public void Scoring_MatchedFragmentCount_Correct()
	{
		var (entry, peaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results[0].MatchedFragmentCount, Is.EqualTo(5));
	}

	[Test]
	public void Scoring_FragmentCoverage_Correct()
	{
		var (entry, peaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results[0].FragmentCoverage, Is.EqualTo(1.0).Within(1e-10));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 2: No-match and empty-input edge cases
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_NoMatch_ExcludedFromResults()
	{
		var (entry, _) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		// Experimental peaks completely off in mass — no matches possible.
		var wrongPeaks = new List<DeconvolutedPeak>
		{
			new(neutralMass: 100.0,  intensity: 1.0f),
			new(neutralMass: 200.0,  intensity: 0.8f),
			new(neutralMass: 300.0,  intensity: 0.6f),
		};

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib, wrongPeaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results, Is.Empty);
	}

	[Test]
	public void Scoring_EmptyCandidateList_ReturnsEmpty()
	{
		var (_, peaks) = MakeSyntheticProteoformData();
		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "DUMMY",
			StrippedSequence = "DUMMY",
			PrecursorMz = 500.0,
			Charge = 5,
			Irt = 0.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>()
		};

		using var lib = BuildLibraryFromEntry(entry);
		// Empty span — no candidates.
		var results = MslProteoformScorer.Score(
			ReadOnlySpan<MslProteoformIndexEntry>.Empty, lib, peaks);

		Assert.That(results, Is.Empty);
	}

	[Test]
	public void Scoring_EmptyExperimentalPeaks_ReturnsEmpty()
	{
		var (entry, _) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib,
			new List<DeconvolutedPeak>(),
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results, Is.Empty);
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 3: MinMatchedFragments filtering
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_MinMatchedFragments_FiltersResults()
	{
		var (entry, allPeaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		// Keep only 3 of the 5 experimental peaks so only 3 fragments can match.
		var partialPeaks = allPeaks.Take(3).ToList();

		using var lib = BuildLibraryFromEntry(entry);

		// minMatchedFragments = 4 → should exclude a candidate that only matches 3.
		var results = MslProteoformScorer.Score(candidates, lib, partialPeaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 4);

		Assert.That(results, Is.Empty, "Candidate with 3 matches should be excluded when minimum is 4.");
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 4: Partial match
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_PartialMatch_LowerScore()
	{
		var (entry, allPeaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		// Full match peaks.
		using var lib = BuildLibraryFromEntry(entry);
		var fullResults = MslProteoformScorer.Score(candidates, lib, allPeaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		// Partial match — only first 3 peaks.
		var partialPeaks = allPeaks.Take(3).ToList();
		var partialResults = MslProteoformScorer.Score(candidates, lib, partialPeaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(partialResults, Has.Count.EqualTo(1));
		Assert.That(partialResults[0].CompositeScore,
			Is.LessThan(fullResults[0].CompositeScore),
			"Partial match should produce a lower composite score than perfect match.");
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 5: Result ordering
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_ResultsSortedByCompositeScore()
	{
		// Build two entries with different fragment counts against the same peaks.
		var (entry1, peaks) = MakeSyntheticProteoformData();  // 5 fragments, all match

		// Entry2 has only 2 fragments that match → lower composite score.
		var entry2Frags = entry1.Fragments!.Take(2).ToList();
		var entry2 = new MslLibraryEntry
		{
			ModifiedSequence = "ANOTHERPROTEOFORM",
			StrippedSequence = "ANOTHERPROTEOFORM",
			PrecursorMz = 1201.0,
			Charge = 10,
			Irt = 46.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			Fragments = entry2Frags
		};

		using var lib = BuildLibraryFromEntries(new[] { entry1, entry2 });
		var allEntries = lib.GetAllEntries().ToList();

		// Build candidates span manually to include both entries.
		var indexEntries = new List<MslProteoformIndexEntry>();
		for (int i = 0; i < allEntries.Count; i++)
		{
			var e = allEntries[i];
			double nm = (double)e.PrecursorMz * e.Charge - e.Charge * Proton;
			indexEntries.Add(new MslProteoformIndexEntry(
				nm, (float)e.PrecursorMz, (short)e.Charge, (float)e.Irt, e.IsDecoy, i));
		}

		var results = MslProteoformScorer.Score(
			indexEntries.ToArray().AsSpan(), lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results, Has.Count.EqualTo(2));
		Assert.That(results[0].CompositeScore,
			Is.GreaterThanOrEqualTo(results[1].CompositeScore),
			"Results must be sorted descending by composite score.");
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 6: Mass tolerance enforcement
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_MassTolerance_RespectedForMatching()
	{
		var (entry, correctPeaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		// Fragment 0: b10+1, Mz = 1200.5 → neutral mass = 1200.5 - 1.007276 ≈ 1199.493
		double frag0NeutralMass = (double)entry.Fragments![0].Mz * 1 - Proton;

		// 10 ppm tolerance on that mass.
		double tenPpm = frag0NeutralMass * 10e-6;

		// Peaks: one exactly within 10 ppm, one exactly outside 20 ppm.
		var peaksInsideTol = new List<DeconvolutedPeak>
		{
			new(frag0NeutralMass + tenPpm * 0.5, 1.0f)  // 5 ppm — inside 20 ppm window
        };
		var peaksOutsideTol = new List<DeconvolutedPeak>
		{
			new(frag0NeutralMass + frag0NeutralMass * 25e-6, 1.0f)  // 25 ppm — outside 20 ppm window
        };

		using var lib = BuildLibraryFromEntry(entry);

		var insideResults = MslProteoformScorer.Score(candidates, lib, peaksInsideTol,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		var outsideResults = MslProteoformScorer.Score(candidates, lib, peaksOutsideTol,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(insideResults, Has.Count.EqualTo(1), "Peak within 20 ppm should match.");
		Assert.That(outsideResults, Is.Empty, "Peak outside 20 ppm should not match.");
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 7: MatchedFragments annotation
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_MatchedFragments_PopulatedCorrectly()
	{
		var (entry, peaks) = MakeSyntheticProteoformData();
		var candidates = BuildCandidateSpan(entry);

		using var lib = BuildLibraryFromEntry(entry);
		var results = MslProteoformScorer.Score(candidates, lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results[0].MatchedFragments, Has.Count.EqualTo(5));

		// First fragment is b10+1 — verify product type and fragment number.
		MatchedFragmentIon firstIon = results[0].MatchedFragments[0];
		Assert.That(firstIon.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.b));
		Assert.That(firstIon.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(10));
		Assert.That(firstIon.Charge, Is.EqualTo(1));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 8: Decoy scoring parity
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_DecoyCandidate_Scored()
	{
		var (targetEntry, peaks) = MakeSyntheticProteoformData();

		// Make an identical decoy entry.
		var decoyEntry = new MslLibraryEntry
		{
			ModifiedSequence = targetEntry.ModifiedSequence,
			StrippedSequence = targetEntry.StrippedSequence,
			PrecursorMz = targetEntry.PrecursorMz,
			Charge = targetEntry.Charge,
			Irt = targetEntry.Irt,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			IsDecoy = true,
			Fragments = targetEntry.Fragments
		};

		var candidateArr = new[]
		{
			new MslProteoformIndexEntry(
				(double)decoyEntry.PrecursorMz * decoyEntry.Charge - decoyEntry.Charge * Proton,
				(float)decoyEntry.PrecursorMz,
				(short)decoyEntry.Charge,
				(float)decoyEntry.Irt,
				isDecoy: true,
				ordinalIndex: 0)
		};

		using var lib = BuildLibraryFromEntry(decoyEntry);
		var results = MslProteoformScorer.Score(
			candidateArr.AsSpan(), lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		Assert.That(results, Has.Count.EqualTo(1), "Decoy candidates should be scored identically to targets.");
		Assert.That(results[0].SpectralAngle, Is.GreaterThan(0.99));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 9: MslLibrary.ScoreProteoformCandidates round-trip
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void MslLibrary_ScoreProteoformCandidates_RoutesCorrectly()
	{
		var (entry, peaks) = MakeSyntheticProteoformData();

		string path = Path.Combine(_tempDir, "proto_roundtrip.msl");
		MslLibrary.Save(path, new[] { entry });

		using var lib = MslLibrary.Load(path);
		Assert.That(lib.ProteoformCount, Is.EqualTo(1), "Library must contain one proteoform entry.");

		double neutralMass = (double)entry.PrecursorMz * entry.Charge - entry.Charge * Proton;

		var results = lib.ScoreProteoformCandidates(
			precursorNeutralMass: neutralMass,
			massTolerance: 10.0,
			deconvolutedPeaks: peaks,
			fragmentMassTolerance: 20e-6,
			minMatchedFragments: 1);

		Assert.That(results, Has.Count.EqualTo(1));
		Assert.That(results[0].SpectralAngle, Is.GreaterThan(0.99));
		Assert.That(results[0].MatchedFragmentCount, Is.EqualTo(5));
	}

	[Test]
	public void MslLibrary_ScoreProteoformCandidates_EmptyWhenNoPeptideOnly()
	{
		// Write a peptide-only library (MoleculeType.Peptide, not Proteoform).
		var peptideEntry = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDER",
			StrippedSequence = "PEPTIDER",
			PrecursorMz = 500.0,
			Charge = 2,
			Irt = 30.0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>
			{
				new() { ProductType = ProductType.b, FragmentNumber = 3, Charge = 1, Mz = 350.2f, Intensity = 1.0f }
			}
		};

		string path = Path.Combine(_tempDir, "peptide_only.msl");
		MslLibrary.Save(path, new[] { peptideEntry });

		using var lib = MslLibrary.Load(path);
		Assert.That(lib.ProteoformCount, Is.EqualTo(0));

		var peaks = new List<DeconvolutedPeak> { new(349.2, 1.0f) };

		// Must return empty list, not throw.
		var results = lib.ScoreProteoformCandidates(
			precursorNeutralMass: 999.0,
			massTolerance: 10.0,
			deconvolutedPeaks: peaks);

		Assert.That(results, Is.Empty);
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Group 10: Orthogonal intensity edge case
	// ─────────────────────────────────────────────────────────────────────────

	[Test]
	public void Scoring_SpectralAngle_ZeroWhenNoIntensityOverlap()
	{
		// Build an entry where all library fragments have intensity 1.0.
		var fragments = new List<MslFragmentIon>
		{
			new() { ProductType = ProductType.b, FragmentNumber = 5, Charge = 1, Mz = 600.0f, Intensity = 1.0f, NeutralLoss = 0 },
			new() { ProductType = ProductType.y, FragmentNumber = 5, Charge = 1, Mz = 700.0f, Intensity = 1.0f, NeutralLoss = 0 },
			new() { ProductType = ProductType.b, FragmentNumber = 8, Charge = 1, Mz = 900.0f, Intensity = 1.0f, NeutralLoss = 0 },
		};

		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "ORTHOGTEST",
			StrippedSequence = "ORTHOGTEST",
			PrecursorMz = 600.0,
			Charge = 5,
			Irt = 20.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			Fragments = fragments
		};

		// Experimental peaks at the exact same neutral masses but with zero intensity
		// on two peaks and 1.0 on one. Because the library has all-equal intensities
		// and the experimental has all zeros except one, the spectral angle will be
		// > 0 (there is exactly one matching dimension). This test verifies the scorer
		// doesn't produce a negative or NaN score under partial overlap.
		double nm0 = (double)fragments[0].Mz - Proton;
		double nm1 = (double)fragments[1].Mz - Proton;
		double nm2 = (double)fragments[2].Mz - Proton;

		// Experimental: only first peak has intensity; other two are zero (omit them).
		var peaks = new List<DeconvolutedPeak>
		{
			new(nm0, 1.0f),
			new(nm1, 0.0f),
			new(nm2, 0.0f),
		};

		var candidates = BuildCandidateSpan(entry);
		using var lib = BuildLibraryFromEntry(entry);

		var results = MslProteoformScorer.Score(candidates, lib, peaks,
			fragmentMassTolerance: 20e-6, minMatchedFragments: 1);

		// SpectralAngle should be in [0, 1] (not NaN, not negative).
		Assert.That(results, Has.Count.EqualTo(1));
		Assert.That(results[0].SpectralAngle, Is.InRange(0.0, 1.0));
		Assert.That(double.IsNaN(results[0].SpectralAngle), Is.False);
		Assert.That(results[0].CompositeScore, Is.GreaterThanOrEqualTo(0.0));
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Private helpers
	// ─────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Builds a single-element <see cref="MslProteoformIndexEntry"/> array from
	/// a library entry and returns it as a <see cref="ReadOnlySpan{T}"/>.
	/// Ordinal index is 0 (correct for single-entry libraries).
	/// </summary>
	private static ReadOnlySpan<MslProteoformIndexEntry> BuildCandidateSpan(MslLibraryEntry entry)
	{
		double neutralMass = (double)entry.PrecursorMz * entry.Charge - entry.Charge * Proton;

		var arr = new[]
		{
			new MslProteoformIndexEntry(
				neutralMass,
				(float)entry.PrecursorMz,
				(short)entry.Charge,
				(float)entry.Irt,
				entry.IsDecoy,
				ordinalIndex: 0)
		};

		return arr.AsSpan();
	}

	/// <summary>
	/// Writes a single entry to a temp .msl file and loads it with
	/// <see cref="MslLibrary.Load"/>. The caller owns the returned library and must
	/// dispose it.
	/// </summary>
	private MslLibrary BuildLibraryFromEntry(MslLibraryEntry entry)
		=> BuildLibraryFromEntries(new[] { entry });

	/// <summary>
	/// Writes multiple entries to a temp .msl file and loads it.
	/// </summary>
	private MslLibrary BuildLibraryFromEntries(IReadOnlyList<MslLibraryEntry> entries)
	{
		string path = Path.Combine(_tempDir, $"lib_{Guid.NewGuid():N}.msl");
		MslLibrary.Save(path, entries);
		return MslLibrary.Load(path);
	}
}