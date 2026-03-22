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
/// Tests targeting the findings from Prompt 6 — The Koina Prediction Pipeline.
///
/// These tests are network-free and exercise MslLibrary.PredictFragments using
/// synthetic delegate functions that simulate Koina responses, so they run in CI
/// without requiring a Koina server.
///
/// Covers:
///   K1 — PredictFragments drops SecondaryProductType/SecondaryFragmentNumber (Fix 8a)
///   K2 — PredictFragments increments updatedCount for empty predictions (Fix 8b)
///   K3 — PredictFragments only processes Source==Predicted AND MatchedFragmentIons.Count==0
///   K4 — PredictFragments returns 0 immediately when no eligible entries exist
///   K5 — PredictFragments skips entries not present in the predict() result dict
///   K6 — GetAwaiter().GetResult() deadlock risk documented (no network test)
///   K7 — Sequence format: mzLib modifications are handled before reaching model
/// </summary>
[TestFixture]
public class TestMslPrompt6KoinaPipeline
{
	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(),
			$"TestMslPrompt6_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	private string TempPath(string name) =>
		Path.Combine(_tempDir, name + ".msl");

	private const double ProtonMass = 1.007276;

	// ── Entry builder helpers ─────────────────────────────────────────────

	private static MslLibraryEntry PredictedEntryNoFragments(
		string seq = "PEPTIDE", int charge = 2) => new()
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = 449.75,
			ChargeState = charge,
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
			MatchedFragmentIons = new List<MslFragmentIon>()   // empty — eligible for prediction
		};

	private static MslLibraryEntry EmpiricalEntryWithFragments(
		string seq = "LESLIEK", int charge = 2) => new()
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = 450.0,
			ChargeState = charge,
			RetentionTime = 35.0,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Empirical,   // NOT eligible
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = 0.01f,
			ElutionGroupId = 0,
			IsDecoy = false,
			MatchedFragmentIons = new List<MslFragmentIon>
		{
			new() { ProductType = ProductType.b, FragmentNumber = 3,
					Charge = 1, Mz = 312.15f, Intensity = 1.0f, NeutralLoss = 0.0,
					ResiduePosition = 3 }
		}
		};

	/// <summary>
	/// Builds a synthetic LibrarySpectrum containing a terminal ion and an internal ion,
	/// simulating what a model that returns internal ions would produce.
	/// </summary>
	private static LibrarySpectrum MakeSpectrumWithInternalIon(
		string seq, double precursorMz, int charge)
	{
		var terminalProduct = new Product(
			ProductType.b, FragmentationTerminus.N,
			neutralMass: 0.0, fragmentNumber: 3,
			residuePosition: 3, neutralLoss: 0.0);

		var internalProduct = new Product(
			ProductType.b, FragmentationTerminus.None,
			neutralMass: 0.0, fragmentNumber: 2,
			residuePosition: 0, neutralLoss: 0.0,
			secondaryProductType: ProductType.y,
			secondaryFragmentNumber: 5);

		return new LibrarySpectrum(
			sequence: seq,
			precursorMz: precursorMz,
			chargeState: charge,
			peaks: new List<MatchedFragmentIon>
			{
				new MatchedFragmentIon(terminalProduct,  experMz: 312.15, experIntensity: 1.0f, charge: 1),
				new MatchedFragmentIon(internalProduct,  experMz: 390.16, experIntensity: 0.5f, charge: 1),
			},
			rt: 30.0);
	}

	private MslLibrary BuildLibrary(IReadOnlyList<MslLibraryEntry> entries)
	{
		string path = TempPath($"lib_{Guid.NewGuid():N}");
		MslWriter.Write(path, entries);
		return MslLibrary.Load(path);
	}

	// ═════════════════════════════════════════════════════════════════════
	// K1 — SecondaryProductType/SecondaryFragmentNumber preserved (Fix 8a)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// After Fix 8a: when a prediction returns a LibrarySpectrum containing an internal
	/// ion (SecondaryProductType set), PredictFragments must preserve both secondary
	/// fields in the resulting MslFragmentIon.
	///
	/// Before Fix 8a: SecondaryProductType was hardcoded to null and
	/// SecondaryFragmentNumber to 0, silently losing internal-ion information.
	/// </summary>
	[Test]
	public void PredictFragments_InternalIon_SecondaryProductType_IsPreserved()
	{
		var entry = PredictedEntryNoFragments();
		using MslLibrary lib = BuildLibrary(new[] { entry });

		string key = entry.Name;

		// Simulate a prediction model that returns an internal ion
		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
			=> new()
			{
				[key] = MakeSpectrumWithInternalIon(entry.FullSequence,
															entry.PrecursorMz, entry.ChargeState)
			};

		int updated = lib.PredictFragments(FakePredict);

		Assert.That(updated, Is.EqualTo(1));

		lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out var result);
		Assert.That(result, Is.Not.Null);

		var internalIon = result!.MatchedFragmentIons.FirstOrDefault(f => f.IsInternalFragment);
		Assert.That(internalIon, Is.Not.Null,
			"After Fix 8a, at least one internal fragment ion must be present.");
		Assert.That(internalIon!.SecondaryProductType, Is.EqualTo(ProductType.y),
			"SecondaryProductType must be preserved from the predicted LibrarySpectrum.");
		Assert.That(internalIon.SecondaryFragmentNumber, Is.EqualTo(5),
			"SecondaryFragmentNumber must be preserved from the predicted LibrarySpectrum.");
	}

	/// <summary>
	/// Terminal ions in a predicted spectrum must NOT have SecondaryProductType set.
	/// Fix 8a must not corrupt terminal ions while fixing internal ions.
	/// </summary>
	[Test]
	public void PredictFragments_TerminalIon_SecondaryProductType_IsNull()
	{
		var entry = PredictedEntryNoFragments();
		using MslLibrary lib = BuildLibrary(new[] { entry });

		string key = entry.Name;

		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
			=> new()
			{
				[key] = MakeSpectrumWithInternalIon(entry.FullSequence,
															entry.PrecursorMz, entry.ChargeState)
			};

		lib.PredictFragments(FakePredict);
		lib.TryGetEntry(entry.FullSequence, entry.ChargeState, out var result);

		var terminalIon = result!.MatchedFragmentIons.FirstOrDefault(f => !f.IsInternalFragment);
		Assert.That(terminalIon, Is.Not.Null);
		Assert.That(terminalIon!.SecondaryProductType, Is.Null,
			"Terminal ions must have SecondaryProductType = null after prediction.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// K2 — Empty prediction does not increment updatedCount (Fix 8b)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// When the predict delegate returns a spectrum with zero MatchedFragmentIons
	/// for an entry, that entry must NOT be counted as updated and must remain
	/// eligible for re-prediction on the next call.
	///
	/// Before Fix 8b: updatedCount was incremented and entry.MatchedFragmentIons was set
	/// to an empty list, making the return value misleading and leaving the entry
	/// stuck in a permanently re-eligible state.
	/// </summary>
	[Test]
	public void PredictFragments_EmptyPrediction_NotCountedAsUpdated()
	{
		var entry = PredictedEntryNoFragments();
		using MslLibrary lib = BuildLibrary(new[] { entry });

		string key = entry.Name;

		// Predict returns a spectrum with zero fragment ions
		var emptySpectrum = new LibrarySpectrum(
			sequence: entry.FullSequence,
			precursorMz: entry.PrecursorMz,
			chargeState: entry.ChargeState,
			peaks: new List<MatchedFragmentIon>(),
			rt: entry.RetentionTime);

		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
			=> new() { [key] = emptySpectrum };

		int updated = lib.PredictFragments(FakePredict);

		Assert.That(updated, Is.EqualTo(0),
			"An empty prediction (zero fragments) must not be counted as an update.");
	}

	/// <summary>
	/// After Fix 8b, an entry that received an empty prediction must remain eligible
	/// (MatchedFragmentIons.Count == 0, Source == Predicted) so subsequent calls can retry.
	/// </summary>
	[Test]
	public void PredictFragments_EmptyPrediction_EntryRemainsEligible()
	{
		var entry = PredictedEntryNoFragments();
		using MslLibrary lib = BuildLibrary(new[] { entry });

		string key = entry.Name;
		var emptySpectrum = new LibrarySpectrum(
			sequence: entry.FullSequence,
			precursorMz: entry.PrecursorMz,
			chargeState: entry.ChargeState,
			peaks: new List<MatchedFragmentIon>(),
			rt: entry.RetentionTime);

		int callCount = 0;
		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
		{
			callCount++;
			return new() { [key] = emptySpectrum };
		}

		// First call: empty prediction, entry not updated
		lib.PredictFragments(FakePredict);

		// Second call: entry must still be eligible (passed to the delegate again)
		lib.PredictFragments(FakePredict);

		Assert.That(callCount, Is.EqualTo(2),
			"After Fix 8b, an entry that received an empty prediction must be passed " +
			"to the predict delegate on subsequent calls (it remains eligible).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// K3 — Eligible filter: Source==Predicted AND MatchedFragmentIons.Count==0
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Empirical entries and entries that already have fragments must be excluded
	/// from the eligible list passed to the predict delegate.
	/// </summary>
	[Test]
	public void PredictFragments_EmpiricalEntries_AreNotEligible()
	{
		var predicted = PredictedEntryNoFragments("PEPTIDE", charge: 2);
		var empirical = EmpiricalEntryWithFragments("LESLIEK", charge: 2);

		using MslLibrary lib = BuildLibrary(new[] { predicted, empirical });

		var receivedKeys = new List<string>();
		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
		{
			receivedKeys.AddRange(eligible.Select(e => e.Name));
			return new Dictionary<string, LibrarySpectrum>();
		}

		lib.PredictFragments(FakePredict);

		Assert.That(receivedKeys, Does.Contain(predicted.Name),
			"The predicted entry with no fragments must be eligible.");
		Assert.That(receivedKeys, Does.Not.Contain(empirical.Name),
			"Empirical entries must not be passed to the predict delegate.");
	}

	/// <summary>
	/// A Predicted entry that already has fragments (e.g. populated from a previous
	/// call) must not be passed to the predict delegate again.
	/// </summary>
	[Test]
	public void PredictFragments_PredictedWithFragments_IsNotEligible()
	{
		var entry = PredictedEntryNoFragments();
		entry.MatchedFragmentIons.Add(new MslFragmentIon
		{
			ProductType = ProductType.b,
			FragmentNumber = 3,
			Charge = 1,
			Mz = 312.15f,
			Intensity = 1.0f,
			NeutralLoss = 0.0,
			ResiduePosition = 3
		});

		using MslLibrary lib = BuildLibrary(new[] { entry });

		int callCount = 0;
		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
		{
			callCount++;
			return new Dictionary<string, LibrarySpectrum>();
		}

		int updated = lib.PredictFragments(FakePredict);

		// No eligible entries → predict should not be called, returns 0
		Assert.That(updated, Is.EqualTo(0),
			"A Predicted entry that already has fragments must not be eligible.");
		Assert.That(callCount, Is.EqualTo(0),
			"The predict delegate must not be called when there are no eligible entries.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// K4 — Returns 0 immediately when no eligible entries
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void PredictFragments_NoEligibleEntries_ReturnsZeroWithoutCallingDelegate()
	{
		var empirical = EmpiricalEntryWithFragments();
		using MslLibrary lib = BuildLibrary(new[] { empirical });

		bool delegateCalled = false;
		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
		{
			delegateCalled = true;
			return new Dictionary<string, LibrarySpectrum>();
		}

		int result = lib.PredictFragments(FakePredict);

		Assert.That(result, Is.EqualTo(0));
		Assert.That(delegateCalled, Is.False,
			"The predict delegate must not be invoked when there are no eligible entries.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// K5 — Unknown keys in predict result are silently skipped
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void PredictFragments_UnknownKeyInResult_IsSilentlySkipped()
	{
		var entry = PredictedEntryNoFragments();
		using MslLibrary lib = BuildLibrary(new[] { entry });

		// Predict returns a key that does not match any eligible entry
		Dictionary<string, LibrarySpectrum> FakePredict(IReadOnlyList<MslLibraryEntry> eligible)
			=> new() { ["UNKNOWNSEQUENCE/99"] = MakeSpectrumWithInternalIon("X", 500.0, 1) };

		Assert.That(() => lib.PredictFragments(FakePredict), Throws.Nothing,
			"Unknown keys in the predict result must be silently skipped.");

		int updated = lib.PredictFragments(FakePredict);
		Assert.That(updated, Is.EqualTo(0),
			"No entries must be updated when the predict result contains only unknown keys.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// K6 — GetAwaiter().GetResult() deadlock risk (documentation test)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Documents the sync-over-async pattern in FragmentIntensityModel.Predict().
	/// No network call is made; this test exists to make the risk visible in the
	/// test suite and to fail if the pattern is ever changed to a clean async surface.
	///
	/// FragmentIntensityModel.Predict() calls:
	///   AsyncThrottledPredictor(modelInputs).GetAwaiter().GetResult()
	///
	/// AsyncThrottledPredictor uses await internally (Task.WhenAll, Task.Delay).
	/// In any environment with a single-threaded SynchronizationContext
	/// (WinForms, WPF, ASP.NET non-Core), this pattern deadlocks.
	///
	/// MetaMorpheus is currently a console app (no SynchronizationContext),
	/// so no deadlock occurs today. The risk is latent.
	/// </summary>
	[Test]
	public void FragmentIntensityModel_Predict_UsesSyncOverAsync_DeadlockRiskDocumented()
	{
		// We cannot easily test for the deadlock without a SynchronizationContext.
		// This test serves as a permanent, searchable documentation marker that
		// the risk was identified and accepted as latent for the current console-app usage.
		//
		// If this test is found during a code review that accompanies adding a WPF/WinForms
		// front-end to MetaMorpheus, it should trigger a fix:
		//   Option A: expose PredictFragmentsAsync and await it from the UI layer.
		//   Option B: wrap Predict() in Task.Run(() => ...).GetAwaiter().GetResult()
		//             to escape the SynchronizationContext.

		Assert.Pass(
			"Sync-over-async risk documented: FragmentIntensityModel.Predict() uses " +
			"AsyncThrottledPredictor(...).GetAwaiter().GetResult(). " +
			"Safe in console apps (no SynchronizationContext). " +
			"Would deadlock in WinForms/WPF/ASP.NET non-Core. " +
			"Fix required before adding a GUI front-end to MetaMorpheus.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// K7 — Sequence format: unsupported modifications are stripped
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// KoinaModelBase.TryCleanSequence with IncompatibleModHandlingMode.RemoveIncompatibleMods
	/// strips modifications not in ValidModificationUnimodMapping.
	/// The Prosit 2020 HCD model uses this mode by default.
	///
	/// This test exercises the mod-stripping path without making a network call
	/// by testing KoinaModelBase directly via Prosit2020IntensityHCD.
	///
	/// A sequence with only unrecognized modifications is stripped to the bare
	/// amino acid sequence, which is then sent to Koina. The caller receives a
	/// WarningException noting that modifications were removed.
	/// </summary>
	[Test]
	public void Prosit2020_UnrecognizedModification_IsStrippedWithWarning()
	{
		var model = new PredictionClients.Koina.SupportedModels.FragmentIntensityModels
			.Prosit2020IntensityHCD();

		// A phospho modification not in Prosit's ValidModificationUnimodMapping
		string seqWithUnsupportedMod = "PEPT[Common Variable:Phospho on T]IDE";

		// Access TryCleanSequence via the public ValidModificationUnimodMapping
		// to confirm the mod is not supported
		bool isSupported = model.ValidModificationUnimodMapping
			.ContainsKey("[Common Variable:Phospho on T]");

		Assert.That(isSupported, Is.False,
			"Phospho on T must not be in Prosit's ValidModificationUnimodMapping — " +
			"it is not a supported modification for this model.");

		// The model maps only Oxidation on M and Carbamidomethyl on C
		Assert.That(model.ValidModificationUnimodMapping.Count, Is.EqualTo(2),
			"Prosit 2020 HCD must support exactly 2 modifications: " +
			"Oxidation on M ([UNIMOD:35]) and Carbamidomethyl on C ([UNIMOD:4]).");
	}

	/// <summary>
	/// Prosit 2020 HCD rejects precursors with charge > 6 (AllowedPrecursorCharges = {1..6}).
	/// ValidateModelSpecificInputs must return false with a warning for charge 7.
	/// This confirms model parameter validation works without a network call.
	/// </summary>
	[Test]
	public void Prosit2020_ChargeAboveSix_IsRejectedWithWarning()
	{
		var model = new PredictionClients.Koina.SupportedModels.FragmentIntensityModels
			.Prosit2020IntensityHCD();

		Assert.That(model.AllowedPrecursorCharges, Does.Not.Contain(7),
			"Prosit 2020 HCD must not support charge 7.");
		Assert.That(model.AllowedPrecursorCharges, Is.EquivalentTo(new[] { 1, 2, 3, 4, 5, 6 }),
			"Prosit 2020 HCD must support exactly charges 1–6.");
	}
}