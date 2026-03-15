using MassSpectrometry;
using NUnit.Framework;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina;                                          // MslFragmentModelRouter + extension
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;  // Prosit2020IntensityHCD
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// NUnit 4 tests for <see cref="MslFragmentModelRouter"/> and the
/// <c>MslLibrary.PredictFragments()</c> extension method.
///
/// Tests that require a live Koina connection are marked
/// [Ignore("Requires Koina server")] and do not run in CI.
/// All routing/grouping/warning tests are network-free.
/// </summary>
[TestFixture]
public class TestMslFragmentModelRouter
{
	// ── Helper ────────────────────────────────────────────────────────────────

	private static MslLibraryEntry MakeEntry(
		string sequence, int charge, DissociationType dt, int nce, double irt = 0.0) =>
		new()
		{
			ModifiedSequence = sequence,
			StrippedSequence = sequence,
			PrecursorMz = 500.0,
			Charge = charge,
			Irt = irt,
			DissociationType = dt,
			Nce = nce,
			Fragments = new List<MslFragmentIon>(),
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			IonMobility = 0.0,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = float.NaN,
			ElutionGroupId = 0,
			IsDecoy = false
		};

	// ── ResolveModelType ──────────────────────────────────────────────────────

	[Test]
	public void Router_HcdWithNce28_ResolvesToProsit2020HCD()
	{
		Type? result = MslFragmentModelRouter.ResolveModelType(DissociationType.HCD, 28);
		Assert.That(result, Is.EqualTo(typeof(Prosit2020IntensityHCD)));
	}

	[Test]
	public void Router_UnknownDissociation_FallsBackToHcd28()
	{
		Type? result = MslFragmentModelRouter.ResolveModelType(DissociationType.Unknown, 0);
		Assert.That(result, Is.EqualTo(typeof(Prosit2020IntensityHCD)));
	}

	[Test]
	public void Router_UnsupportedDissociation_ReturnsNull()
	{
		Assert.That(MslFragmentModelRouter.ResolveModelType(DissociationType.ETD, 25), Is.Null);
	}

	// ── IsSupported ───────────────────────────────────────────────────────────

	[Test]
	public void Router_IsSupported_TrueForHcd()
	{
		Assert.That(MslFragmentModelRouter.IsSupported(DissociationType.HCD, 28), Is.True);
	}

	[Test]
	public void Router_IsSupported_FalseForUnregistered()
	{
		Assert.That(MslFragmentModelRouter.IsSupported(DissociationType.ETD, 25), Is.False);
	}

	// ── Grouping / normalisation ──────────────────────────────────────────────

	[Test]
	public void Router_GroupsByDissociationTypeAndNce()
	{
		// HCD-28 and HCD-35 are distinct groups but both route to the same model type
		Assert.That(MslFragmentModelRouter.IsSupported(DissociationType.HCD, 28), Is.True);
		Assert.That(MslFragmentModelRouter.IsSupported(DissociationType.HCD, 35), Is.True);
		Assert.That(
			MslFragmentModelRouter.ResolveModelType(DissociationType.HCD, 28),
			Is.EqualTo(MslFragmentModelRouter.ResolveModelType(DissociationType.HCD, 35)));
	}

	[Test]
	public void Router_ZeroNce_UsesDefault()
	{
		Assert.That(MslFragmentModelRouter.IsSupported(DissociationType.HCD, 0), Is.True);
		Assert.That(
			MslFragmentModelRouter.ResolveModelType(DissociationType.HCD, 0),
			Is.EqualTo(typeof(Prosit2020IntensityHCD)));
	}

	// ── Warning behaviour ─────────────────────────────────────────────────────

	[Test]
	public void Router_UnsupportedEntry_SkippedWithWarning()
	{
		// ETD is unregistered → rejected at factory-lookup, no HTTP call made
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE", charge: 2, dt: DissociationType.ETD, nce: 25)
		};

		var results = MslFragmentModelRouter.Predict(entries, out string? warnings);

		Assert.That(warnings, Is.Not.Null);
		Assert.That(results, Is.Empty);
	}

	[Test]
	public void Router_EmptyInput_NoWarnings()
	{
		var results = MslFragmentModelRouter.Predict(
			new List<MslLibraryEntry>(), out string? warnings);

		Assert.That(warnings, Is.Null);
		Assert.That(results, Is.Empty);
	}

	// ── Live Koina tests (ignored in CI) ─────────────────────────────────────

	[Test]
	[Ignore("Requires Koina server — run manually with network access")]
	public void Router_LivePrediction_SmallBatch()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE",    charge: 2, dt: DissociationType.HCD, nce: 28),
			MakeEntry("LESLIEK",    charge: 2, dt: DissociationType.HCD, nce: 28),
			MakeEntry("ELVISLIVER", charge: 3, dt: DissociationType.HCD, nce: 28)
		};

		var results = MslFragmentModelRouter.Predict(entries, out string? warnings);

		Assert.That(results, Has.Count.EqualTo(3));
		foreach (var spectrum in results.Values)
			Assert.That(spectrum.MatchedFragmentIons, Is.Not.Empty,
				$"Spectrum {spectrum.Name} should have fragment ions");
		Assert.That(warnings, Is.Null);
	}

	[Test]
	[Ignore("Requires Koina server — run manually with network access")]
	public void Router_LivePrediction_NceAffectsIntensities()
	{
		var e20 = new List<MslLibraryEntry> { MakeEntry("PEPTIDE", 2, DissociationType.HCD, 20) };
		var e35 = new List<MslLibraryEntry> { MakeEntry("PEPTIDE", 2, DissociationType.HCD, 35) };

		var r20 = MslFragmentModelRouter.Predict(e20, out _);
		var r35 = MslFragmentModelRouter.Predict(e35, out _);

		Assert.That(r20, Has.Count.EqualTo(1));
		Assert.That(r35, Has.Count.EqualTo(1));

		double sum20 = r20.Values.First().MatchedFragmentIons.Sum(f => f.Intensity);
		double sum35 = r35.Values.First().MatchedFragmentIons.Sum(f => f.Intensity);

		Assert.That(sum20, Is.Not.EqualTo(sum35).Within(1e-6),
			"NCE 20 and NCE 35 should produce different intensity distributions");
	}

	[Test]
	[Ignore("Requires Koina server — run manually with network access")]
	public void MslLibrary_PredictFragments_UpdatesFragments()
	{
		const string testLibraryPath = "TestData/predicted_no_fragments.msl";

		Assert.That(File.Exists(testLibraryPath), Is.True,
			$"Pre-condition: test file must exist at {testLibraryPath}");

		using var library = MslLibrary.Load(testLibraryPath);

		int emptyBefore = library.GetAllEntries()
			.Count(e => e.Source == MslFormat.SourceType.Predicted && e.Fragments.Count == 0);

		Assert.That(emptyBefore, Is.GreaterThan(0));

		// Extension method from MslLibraryPredictionExtensions
		int updated = library.PredictFragments();

		Assert.That(updated, Is.EqualTo(emptyBefore));
		Assert.That(
			library.GetAllEntries()
				.Count(e => e.Source == MslFormat.SourceType.Predicted && e.Fragments.Count == 0),
			Is.Zero);
	}
}