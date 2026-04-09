using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// Tests targeting the findings from Prompt 3 — SpectralAngle vs. Cosine Similarity.
/// Written for the corrective Fix 6B: ComputeSpectralAngle now returns the true
/// arccos-scaled SA = 1 − (2/π) × arccos(cosine).
///
/// All existing tests in TestMslProteoformScoring that use identical/near-identical
/// vectors are unaffected because cosine == SA == 1.0 for those cases.
/// </summary>
[TestFixture]
public class TestMslPrompt3SpectralAngle
{
	private const double Proton = 1.007276;

	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(),
			$"TestMslPrompt3_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	// ── Helpers — same pattern as Testmslproteoformscoring.cs ────────────

	private MslLibrary BuildLibraryFromEntry(MslLibraryEntry entry)
	{
		string path = Path.Combine(_tempDir, $"lib_{Guid.NewGuid():N}.msl");
		MslLibrary.Save(path, new[] { entry });
		return MslLibrary.Load(path);
	}

	/// <summary>
	/// Constructs a ReadOnlySpan of one MslProteoformIndexEntry for the given
	/// entry, with OrdinalIndex = 0 (correct for single-entry libraries).
	/// Mirrors BuildCandidateSpan from Testmslproteoformscoring.cs.
	/// </summary>
	private static ReadOnlySpan<MslProteoformIndexEntry> CandidateSpan(MslLibraryEntry entry)
	{
		double neutralMass = (double)entry.PrecursorMz * entry.ChargeState
							 - entry.ChargeState * Proton;
		return new[]
		{
			new MslProteoformIndexEntry(
				neutralMass,
				(float)entry.PrecursorMz,
				(short)entry.ChargeState,
				(float)entry.RetentionTime,
				entry.IsDecoy,
				ordinalIndex: 0)
		}.AsSpan();
	}

	/// <summary>
	/// Builds a library entry with fragments whose neutral masses exactly match
	/// the experimental peaks, paired with the given intensity vectors.
	/// Each fragment gets a distinct m/z (300, 400, 500, …) with charge=1.
	/// </summary>
	private static (MslLibraryEntry entry, List<DeconvolutedPeak> peaks)
		MakePair(float[] libIntensities, float[] expIntensities)
	{
		var frags = new List<MslFragmentIon>();
		var peaks = new List<DeconvolutedPeak>();

		for (int i = 0; i < libIntensities.Length; i++)
		{
			float mz = 300.0f + i * 100.0f;
			double neutralMass = mz * 1 - 1 * Proton;

			frags.Add(new MslFragmentIon
			{
				ProductType = ProductType.b,
				FragmentNumber = i + 1,
				Charge = 1,
				Mz = mz,
				Intensity = libIntensities[i],
				NeutralLoss = 0.0,
				ResiduePosition = i + 1
			});

			peaks.Add(new DeconvolutedPeak(neutralMass, expIntensities[i]));
		}

		var entry = new MslLibraryEntry
		{
			FullSequence = "TESTPROTEOFORM",
			BaseSequence = "TESTPROTEOFORM",
			PrecursorMz = 1200.0,
			ChargeState = 10,
			RetentionTime = 45.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = frags
		};

		return (entry, peaks);
	}

	private MslProteoformScoringResult ScorePair(
		float[] libIntensities,
		float[] expIntensities,
		int minFragments = 1)
	{
		var (entry, peaks) = MakePair(libIntensities, expIntensities);

		using MslLibrary lib = BuildLibraryFromEntry(entry);

		var results = MslProteoformScorer.Score(
			CandidateSpan(entry), lib, peaks,
			fragmentMassTolerance: 100e-6,
			minMatchedFragments: minFragments);

		Assert.That(results, Is.Not.Empty,
			"Score must return at least one result for a matched entry.");

		return results[0];
	}

	// ═════════════════════════════════════════════════════════════════════
	// S1 — The discriminating test
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// THE DEFINITIVE TEST.
	///
	/// For lib=[1,0] and exp=[1,1]:
	///   cosine similarity = 1/√2 ≈ 0.7071
	///   true SA           = 1 − (2/π)×arccos(1/√2) = 0.5000 exactly
	///
	/// After Fix 6B, SpectralAngle must return 0.5, not 0.7071.
	/// </summary>
	[Test]
	public void SpectralAngle_LibOneZero_ExpOneOne_ReturnsTrueSA_NotCosine()
	{
		float[] lib = { 1.0f, 0.0f };
		float[] exp = { 1.0f, 1.0f };

		const double expectedSA = 0.5;
		double expectedCosine = 1.0 / Math.Sqrt(2.0);   // ≈ 0.7071

		var result = ScorePair(lib, exp);

		Assert.That(result.SpectralAngle,
			Is.EqualTo(expectedSA).Within(1e-9),
			$"After Fix 6B, SpectralAngle([1,0],[1,1]) must equal true SA = {expectedSA}. " +
			$"Got {result.SpectralAngle:F6}. If ~0.7071 the fix has not been applied.");

		Assert.That(Math.Abs(result.SpectralAngle - expectedCosine),
			Is.GreaterThan(0.1),
			$"SpectralAngle must not equal cosine similarity ({expectedCosine:F4}).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// S2 — Identical vectors: SA = 1.0 (unchanged by fix)
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void SpectralAngle_IdenticalVectors_IsOne()
	{
		float[] intensities = { 0.8f, 0.6f, 0.4f, 0.9f, 1.0f };
		var result = ScorePair(intensities, intensities);

		// float32 intensities accumulate small rounding errors in the dot product,
		// so the cosine is not exactly 1.0 (e.g. 0.9999999905...). After the arccos
		// transform, SA is correspondingly close to but not exactly 1.0.
		// 1e-6 is sufficient to confirm this is a float32 precision artefact, not
		// a computation error, and remains far from any meaningfully wrong value.
		Assert.That(result.SpectralAngle, Is.EqualTo(1.0).Within(1e-6),
			"Identical intensity vectors must yield SpectralAngle = 1.0 (within float32 precision).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// S3 — Orthogonal vectors: SA = 0.0 (unchanged by fix)
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void SpectralAngle_OrthogonalVectors_IsZero()
	{
		float[] lib = { 1.0f, 0.0f };
		float[] exp = { 0.0f, 1.0f };

		var result = ScorePair(lib, exp);

		Assert.That(result.SpectralAngle, Is.EqualTo(0.0).Within(1e-9),
			"Orthogonal intensity vectors must yield SpectralAngle = 0.0.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// S4 — Intermediate value: SA ≠ cosine
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// For lib=[1,0.5] and exp=[0.8,0.6]:
	///   cosine ≈ 0.9839 → true SA ≈ 0.8855
	/// After Fix 6B the result must be near 0.8855, not near 0.9839.
	/// </summary>
	[Test]
	public void SpectralAngle_PartialMatch_ReturnsTrueSA()
	{
		float[] lib = { 1.0f, 0.5f };
		float[] exp = { 0.8f, 0.6f };

		double dot = 1.0 * 0.8 + 0.5 * 0.6;
		double normLib = Math.Sqrt(1.0 * 1.0 + 0.5 * 0.5);
		double normExp = Math.Sqrt(0.8 * 0.8 + 0.6 * 0.6);
		double cosine = dot / (normLib * normExp);
		double expectedSA = 1.0 - (2.0 / Math.PI) * Math.Acos(cosine);

		var result = ScorePair(lib, exp);

		Assert.That(result.SpectralAngle,
			Is.EqualTo(expectedSA).Within(1e-4),
			$"SpectralAngle must equal true SA ≈ {expectedSA:F4}, not cosine ≈ {cosine:F4}.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// S5 — arccos clamp guard: result never NaN
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void SpectralAngle_IsNeverNaN_ClampGuardWorks()
	{
		float[] intensities = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
		var result = ScorePair(intensities, intensities);

		Assert.That(double.IsNaN(result.SpectralAngle), Is.False,
			"SpectralAngle must never be NaN — Math.Clamp guard prevents arccos(>1).");
		Assert.That(result.SpectralAngle, Is.InRange(0.0, 1.0),
			"SpectralAngle must always be in [0, 1].");
	}

	// ═════════════════════════════════════════════════════════════════════
	// S6 — Composite score uses true SA
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void CompositeScore_EqualsSpectralAngle_Times_Log2FragCount()
	{
		float[] lib = { 1.0f, 0.5f };
		float[] exp = { 0.8f, 0.6f };

		var result = ScorePair(lib, exp);

		double expected = result.SpectralAngle * Math.Log2(1.0 + result.MatchedFragmentCount);

		Assert.That(result.CompositeScore, Is.EqualTo(expected).Within(1e-9),
			"CompositeScore must equal SpectralAngle × log2(1 + MatchedFragmentCount).");
	}

	[Test]
	public void CompositeScore_PerfectMatch_OneFragment_IsOne()
	{
		var result = ScorePair(new float[] { 1.0f }, new float[] { 1.0f });

		Assert.That(result.SpectralAngle, Is.EqualTo(1.0).Within(1e-9));
		Assert.That(result.CompositeScore, Is.EqualTo(1.0).Within(1e-9),
			"SA=1.0 × log2(2)=1.0 = 1.0.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// S7 — Neutral mass conversion correctness
	// ═════════════════════════════════════════════════════════════════════

	[TestCase(1, TestName = "NeutralMassConversion_z1")]
	[TestCase(2, TestName = "NeutralMassConversion_z2")]
	[TestCase(3, TestName = "NeutralMassConversion_z3")]
	public void NeutralMassConversion_FindsMatchAtCorrectNeutralMass(int charge)
	{
		float mz = 400.0f;
		double expectedNeutralMass = (double)mz * charge - charge * Proton;

		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 1200.0,
			ChargeState = 10,
			RetentionTime = 45.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new()
				{
					ProductType = ProductType.b, FragmentNumber = 3,
					Charge = charge, Mz = mz, Intensity = 1.0f,
					NeutralLoss = 0.0, ResiduePosition = 3
				}
			}
		};

		using MslLibrary lib = BuildLibraryFromEntry(entry);

		var peak = new DeconvolutedPeak(expectedNeutralMass, 1.0f);

		var results = MslProteoformScorer.Score(
			CandidateSpan(entry), lib,
			new List<DeconvolutedPeak> { peak },
			fragmentMassTolerance: 5e-6,
			minMatchedFragments: 1);

		Assert.That(results, Is.Not.Empty,
			$"z={charge}: must find match at neutral_mass = mz×z − z×protonMass = {expectedNeutralMass:F4}.");
		Assert.That(results[0].MatchedFragmentCount, Is.EqualTo(1));
	}

	// ═════════════════════════════════════════════════════════════════════
	// S8 — ChargeState = 0 guard
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// After Fix 6B, a Charge=0 fragment is skipped by the guard in ScoreCandidate.
	/// The valid z=1 fragment in the same entry still produces a match.
	/// </summary>
	[Test]
	public void ScoreCandidate_ChargeZeroFragment_IsSkipped()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 1200.0,
			ChargeState = 10,
			RetentionTime = 45.0,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new()
				{
					ProductType = ProductType.b, FragmentNumber = 1,
					Charge = 1, Mz = 300.0f, Intensity = 1.0f,   // valid — will match
                    NeutralLoss = 0.0, ResiduePosition = 1
				},
				new()
				{
					ProductType = ProductType.b, FragmentNumber = 2,
					Charge = 0, Mz = 400.0f, Intensity = 0.8f,   // ChargeState=0 — must be skipped
                    NeutralLoss = 0.0, ResiduePosition = 2
				}
			}
		};

		using MslLibrary lib = BuildLibraryFromEntry(entry);

		double validNeutralMass = 300.0 * 1 - 1 * Proton;
		var peaks = new List<DeconvolutedPeak>
		{
			new DeconvolutedPeak(validNeutralMass, 1.0f)
		};

		var results = MslProteoformScorer.Score(
			CandidateSpan(entry), lib, peaks,
			fragmentMassTolerance: 10e-6,
			minMatchedFragments: 1);

		Assert.That(results, Is.Not.Empty, "The z=1 fragment must still match.");
		Assert.That(results[0].MatchedFragmentCount, Is.EqualTo(1),
			"Only 1 fragment must match — the z=0 fragment must be skipped by the guard.");
	}

	[Test]
	public void ValidateEntries_ChargeZeroFragment_ReportsError()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 500.0,
			ChargeState = 2,
			RetentionTime = 30.0,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Predicted,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = float.NaN,
			ElutionGroupId = 0,
			IsDecoy = false,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new()
				{
					ProductType = ProductType.b, FragmentNumber = 1,
					Charge = 0, Mz = 300.0f, Intensity = 1.0f,
					NeutralLoss = 0.0, ResiduePosition = 1
				}
			}
		};

		List<string> errors = MslWriter.ValidateEntries(
			new List<MslLibraryEntry> { entry });

		Assert.That(errors.Any(e => e.Contains("ChargeState")), Is.True,
			"ValidateEntries must report a ChargeState error for a fragment with ChargeState=0.");
	}
}